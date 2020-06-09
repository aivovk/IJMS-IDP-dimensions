#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <csignal>

#include "Display.h"

#include "World.h"
#include "Vector3D.h"
#include "NormalDistribution.h"

TYPE_FLOAT ** outputArray; ///< buffer for file output
int numLines; ///< number of lines to write at once
int numEntries; ///< entries per line

void processSIGTERM(int sig)
{
  std::cout<< "SIGTERM RECEIVED!" << std::endl;
  WorldSettings::terminate = true;
}

void processSIGCheckpoint(int sig)
{
  std::cout<< "SIGCheckpoint RECEIVED!" << std::endl;
  WorldSettings::saveCheckpoint = true;
}

void writeOutput(std::ofstream * _file, int currentLine)
{
  std::stringstream outputBlock;
  for (int iLine = 0 ; iLine < currentLine ; iLine++){
    for (int iEntry = 0 ; iEntry < numEntries ; iEntry++){
      outputBlock<<outputArray[iLine][iEntry]<<",";
    }
    outputBlock<<std::endl;
  }
  (*_file) << outputBlock.str();
}

/// \todo check that references are correct
void writeOutputBinary(std::ofstream * _file, int currentLine){
  (*_file).write((char *) &(outputArray[0][0]), currentLine*numEntries*sizeof(TYPE_FLOAT));
}
		      

int main ( int argc, char** argv )
{ 
  signal(SIGTERM, processSIGTERM);
  signal(SIGINT, processSIGCheckpoint);

  int isFinished = 0;

  // check if a different configuration file has been specified
  if (argc == 1)
    WorldSettings::initialize();
  else
    for (int i = 1; i < argc; i+=2) /// \todo why inc by 2?
      WorldSettings::initialize(argv[i]);
  
  TYPE_FLOAT realDT = WorldSettings::dt;
  World world;
    
  if (WorldSettings::outSetting == WorldSettings::FILE)
    {
      // file write array
      int currentLine = 0;
      numEntries = 13;
      numLines = 20;
      outputArray = new TYPE_FLOAT*[numLines];
      for(int i = 0; i < numLines; ++i)
	{
	  outputArray[i] = new TYPE_FLOAT[numEntries];
	}
      
      /* open the file for output */
      std::stringstream filenamestream;
      filenamestream << WorldSettings::outFile << ".txt";
      std::string filename = filenamestream.str();
      std::ofstream file;
      file.open(filename.c_str(), std::fstream::out | std::fstream::app );
      //file.open(filename.c_str(), std::ios::out | std::ios::app | std::ios::binary);
      file<<"#seed: "<<WorldSettings::seed<<std::endl;
      
      std::cout<<"Writing to: "<<filename<<std::endl;
      
      while (world.getStep() < WorldSettings::steps && !isFinished)
        {
	  //world.resetDisplacementSquared();
	  realDT = world.simulate();

	  // make a snapshot of current configuration
	  if(WorldSettings::stepsBetweenSnapshot > 0 &&
	     world.getStep() % WorldSettings::stepsBetweenSnapshot ==0){
	    world.savePositions(); 
	    world.saveStateBinaryFast(WorldSettings::snapshotFile);
	  }
	  
	  if(world.getStep() % WorldSettings::stepsBetweenFileOutput == 0)
	    {
	      if (WorldSettings::stepsBetweenSeparation > 0)
		world.updateSeparationIJ(realDT);
	      
	      Vector3D com = world.centreOfMass() - WorldSettings::offsetVector;
	      Vector3D re = world.endToEndVector();
	      Vector3D eigen_rg = world.asphericity();
	 
	      // write the output array to file and reset the counter
	      if ( currentLine == numLines) {
		writeOutput(&file, currentLine);
		currentLine = 0;
	      }
	      
	      outputArray[currentLine][0] = world.getTime();
	      //hydrodynamic radius of the first partile
	      outputArray[currentLine][1] = world.displacementSquared(0)/realDT;
	      outputArray[currentLine][2] = world.rkInverse();
	      outputArray[currentLine][3] = com.x;
	      outputArray[currentLine][4] = com.y;
	      outputArray[currentLine][5] = com.z;
	      outputArray[currentLine][6] = re.x;
	      outputArray[currentLine][7] = re.y;
	      outputArray[currentLine][8] = re.z;
	      outputArray[currentLine][9] = eigen_rg.x;
	      outputArray[currentLine][10] = eigen_rg.y;
	      outputArray[currentLine][11] = eigen_rg.z;
	      outputArray[currentLine][12] = world.averageBondLength();

	      currentLine++;  
		  
	    }

	  if(WorldSettings::stepsBetweenSeparation > 0 &&
	     world.getStep() % WorldSettings::stepsBetweenSeparation ==0){
	    world.saveSeparationIJ(WorldSettings::separationFile);
	  }
	  	  
	  if(WorldSettings::terminate)
	    {
	      //world.saveContactMap(WorldSettings::contactMapFile);
	      world.saveSeparationIJ(WorldSettings::separationFile);
	      world.savePositions();
	      world.saveStateBinary(WorldSettings::checkpointFile);
	      isFinished = 1;
	    }
          if(WorldSettings::saveCheckpoint)
	    {
	      //world.saveContactMap(WorldSettings::contactMapFile);
	      world.savePositions();
	      world.saveStateBinary(WorldSettings::checkpointFile);
	      WorldSettings::saveCheckpoint = false;
	    }
	  
        }
      // write the incomplete output array to file
      writeOutput(&file, currentLine);
      file.close();

      for(int i = 0; i < numLines; ++i) {
	delete [] outputArray[i];
      }
      delete [] outputArray;
      
    } 
  else if (WorldSettings::outSetting == WorldSettings::SCREEN)
    {   
      #ifdef SFML

      Display displayWindow;
             
      while(world.getStep() < WorldSettings::steps && !isFinished)
        {
	  displayWindow.handleEvents();

	  world.simulate();
	  	  
	  displayWindow.draw(world);

	  if(WorldSettings::terminate)
	    {
	      world.savePositions();
	      world.saveStateBinary(WorldSettings::checkpointFile);
	      isFinished = 1;
	    }
	  if(WorldSettings::saveCheckpoint)
	    {
	      world.savePositions();
	      world.saveStateBinary(WorldSettings::checkpointFile);
	      WorldSettings::saveCheckpoint = false;
	    }
        }
      #endif
    } 

  /*  if the simulation finishes all steps, it did not receive a SIGTERM,
   *  and the terminate flag was not set, save a checkpoint of the state
   *  at the last file write (will not work with SCREEN)
   */
  if(!WorldSettings::terminate){
    world.saveSeparationIJ(WorldSettings::separationFile);
    world.savePositions();
    world.saveStateBinary(WorldSettings::checkpointFile);
  }

  delete WorldSettings::normalDistribution;
    
  std::cout<<"LJ Error Count: "<<WorldSettings::errorCountLJ<<
    " out of "<<WorldSettings::totalCountLJ<<std::endl;
  std::cout<<"FENE Error Count: "<<WorldSettings::errorCountFENE<<
    " out of "<<WorldSettings::totalCountFENE<<std::endl;
  return isFinished;
}
