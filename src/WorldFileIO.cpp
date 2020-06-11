#include "World.h"

void World::writeContactMap(std::string filename)
{
  std::ofstream file;
  file.open(filename.c_str());
  file<<noOfParticles<<std::endl;
  for (int i = 0; i < noOfParticles; i++)
    {
      for (int j = 0 ; j < noOfParticles ; j++)
	{
	  file<<contactMap[i][j]<<",";
	}
      file<<std::endl;
    }
  file.close();
}

void World::writeSeparationIJ(std::string filename)
{
  std::ofstream file;
  file.open(filename.c_str());
  file<<"noOfParticles,"<<noOfParticles<<std::endl;
  file<<"entries,"<<noOfSepEntries<<std::endl;
  for (int i = 0; i < noOfParticles; i++)
    {
      file<<separationIJ[i]<<",";
    }
  file<<std::endl;
  for (int i = 0; i < noOfParticles; i++)
    {
      file<<separationSquaredIJ[i]<<",";
    }
  file<<std::endl;
  file.close();
}

void World::writeState(std::string filename)
{
  std::ofstream file;
  file.open(filename.c_str());
  file<<noOfParticles<<std::endl;
  file<<"Step "<<savedStep<<std::endl;
  file<<"Time "<<savedTime<<std::endl;
  file<<"Checkpoint "<<checkpointNumber<<std::endl;
  for(int i = 0; i < noOfParticles ; i++)
    file<<particles[i]->getAACode()<<" "<<savedPositions[i].x<<" "<<savedPositions[i].y<<" "<<savedPositions[i].z<<std::endl;
  file.close();

  std::ostringstream os;
  os << checkpointNumber;

  std::ifstream src(filename.c_str(), std::ios::binary);
  std::string copyFilename = filename + os.str();
  std::ofstream dst(copyFilename.c_str(), std::ios::binary);
  dst << src.rdbuf();

  stateSaved = true;
  checkpointNumber++;
}


//similar to the below code except it only writes the timestep and then the
//particle positions and appends to the same file
void World::writeStateBinaryFast(std::string filename)
{
  std::ofstream file;
  file.open(filename.c_str(), std::ios::out | std::ios::app | std::ios::binary);
  file.write((char *) &savedTime, sizeof(TYPE_FLOAT));
  //size of Vector3D is 24 (3 x sizeof(double)) so nothing extra
  //this is more than twice as fast as iterating:
  file.write((char *) &(savedPositions[0]), noOfParticles*sizeof(Vector3D));
  file.close();
}

// \todo buffer?
// pure C functions are much faster?
void World::writeStateBinary(std::string filename)
{
  std::ofstream file;
  file.open(filename.c_str(), std::ios::out | std::ios::binary | std::ios::trunc);
  file.write((char *) &noOfParticles, sizeof(int));
  file.write((char *) &savedStep, sizeof(int));
  file.write((char *) &savedTime, sizeof(TYPE_FLOAT));
  file.write((char *) &checkpointNumber, sizeof(int));
  for(int i = 0 ; i < noOfParticles ; i++)
    {
      char AAcode = particles[i]->getAACode();
      file.write(&AAcode, 1);
      file.write((char *) &(savedPositions[i]), sizeof(Vector3D));
    }
  file.close();

  // copies the checkpoint to a file with the same name and a checkpoint number
  // appended at the end, this is useful if you want to save multiple checkpoints
  // throughout a simulation
  // 
  std::ostringstream os;
  os << checkpointNumber;

  std::ifstream src(filename.c_str(), std::ios::binary);
  std::string copyFilename = filename + os.str();
  std::ofstream dst(copyFilename.c_str(), std::ios::binary);

  dst << src.rdbuf();

  stateSaved = true;
  checkpointNumber++;
}

void World::loadState(std::string filename)
{
  std::ifstream file;
  file.open(filename.c_str());
  if(!file.is_open())
    return;

  std::string line;
  std::getline(file, line);
  std::getline(file, line);
  std::istringstream fss(line);
  std::string val;
  std::getline(fss, val, ' ');
  std::getline(fss, val);
  step = atoi(val.c_str());
  
  std::istringstream gss(line);
  std::getline(gss, val, ' ');
  std::getline(gss, val);
  time = atof(val.c_str());

  std::istringstream hss(line);
  std::getline(hss, val, ' ');
  std::getline(hss, val);
  checkpointNumber = atoi(val.c_str());
  
  char AA;
  TYPE_FLOAT X, Y, Z;
  int j = 0;
  while(std::getline(file, line))
    {    
      std::istringstream ss(line);
      std::getline(ss, val, ' ');
      AA = val.c_str()[0];
      if( AA == particles[j]->getAACode() )
	std::cout << "Loading: " <<AA<<std::endl;
      else
	std::cout << "Warning, AA expected: "<<particles[j]->getAACode()<<", AA in checkpoint file: "<<AA <<std::endl;

      std::getline(ss, val, ' ');
      X = atof(val.c_str());
      std::getline(ss, val, ' ');
      Y = atof(val.c_str());
      std::getline(ss, val);
      Z = atof(val.c_str());
      particles[j]->r = Vector3D(X, Y, Z)+WorldSettings::offsetVector;
      previous[j] = particles[j]->r;
      j++;
    }
  
  file.close();
}

void World::loadStateBinary(std::string filename)
{
  std::ifstream file;
  file.open(filename.c_str(), std::ios::in | std::ios::binary);
  if(!file.is_open())
    return;
  int numParticles;
  
  file.read((char *) &numParticles, sizeof(int));
  file.read((char *) &step, sizeof(int));
  file.read((char *) &time, sizeof(TYPE_FLOAT));
  file.read((char *) &checkpointNumber, sizeof(int));
  
  char AA;
  Vector3D R;
  for( int i = 0 ; i < numParticles ; i++ )
    {
      file.read(&AA, 1);
      if( AA == particles[i]->getAACode() )
	std::cout << "Loading: " <<AA<<std::endl;
      else
	std::cout << "Warning, AA expected: "<<particles[i]->getAACode()<<", AA in checkpoint file: "<<AA <<std::endl;
      file.read((char *) &R, sizeof(Vector3D));
      particles[i]->r = R + WorldSettings::offsetVector;
      previous[i] = particles[i]->r;
    }
  file.close();
}
