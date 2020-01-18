#include "WorldSettings.h"


const char * WorldSettings::defaultSettingsFile = "WorldSettings.conf";

std::string WorldSettings::outFile;
std::string WorldSettings::snapshotFile;
std::string WorldSettings::separationFile;
int WorldSettings::stepsBetweenFileOutput;
int WorldSettings::stepsBetweenSnapshot;
int WorldSettings::stepsBetweenSeparation;
std::string WorldSettings::nupFile;
std::string WorldSettings::AAFile;
std::string WorldSettings::contactMapFile;
std::string WorldSettings::checkpointFile;
int WorldSettings::checkpointCounter;
NormalDistribution * WorldSettings::nd;
PeriodicBoundary WorldSettings::pbc;
AA *WorldSettings::aaProperties;
unsigned int WorldSettings::seed;
int WorldSettings::skipNeighbourUpdate;
TYPE_FLOAT WorldSettings::LJSize;
TYPE_FLOAT WorldSettings::bondLength;
TYPE_FLOAT WorldSettings::avgBondLength;
TYPE_FLOAT WorldSettings::debyeLength;
TYPE_FLOAT WorldSettings::feneLength;
TYPE_FLOAT WorldSettings::maxFENELengthSquared;
TYPE_FLOAT WorldSettings::offset;
Vector3D WorldSettings::offsetVector;
TYPE_FLOAT WorldSettings::coulombStrength;
TYPE_FLOAT WorldSettings::hydrophobicStrength;
TYPE_FLOAT WorldSettings::eLJ;
TYPE_FLOAT WorldSettings::dt;
TYPE_FLOAT WorldSettings::sqrt_dt;
int WorldSettings::steps;
int WorldSettings::eqSteps;
int WorldSettings::skip;
int WorldSettings::noOfMonomers;

TYPE_FLOAT WorldSettings::hydrophobicCutoff;
TYPE_FLOAT WorldSettings::neighbourListBuffer;

WorldSettings::TYPE_OUTPUT WorldSettings::outSetting;
WorldSettings::TYPE_FORCE_CHARGE WorldSettings::typeForceCharge;
WorldSettings::TYPE_FORCE_REPULSIVE WorldSettings::typeForceRepulsive;
WorldSettings::TYPE_FORCE_HYDROPHOBIC WorldSettings::typeForceHydrophobic;
WorldSettings::TYPE_FORCE_BOND WorldSettings::typeForceBond;
WorldSettings::TYPE_HYDRO WorldSettings::typeHydrodynamics;
WorldSettings::TYPE_FORCE_EXTERNAL WorldSettings::typeForceExternal;

WorldSettings::TYPE_INITIAL_CONDITION WorldSettings::initialCondition;

//flags
bool WorldSettings::terminate;
bool WorldSettings::saveCheckpoint;

//error totals
unsigned int WorldSettings::errorCountFENE;
unsigned int WorldSettings::errorCountLJ;
unsigned int WorldSettings::totalCountFENE;
unsigned int WorldSettings::totalCountLJ;

// error maxes
TYPE_FLOAT WorldSettings::MAX_FORCE;

WorldSettings::WorldSettings(){}

WorldSettings::~WorldSettings(){}

void WorldSettings::initialize()
{
  initialize(defaultSettingsFile);
}

void WorldSettings::initialize(const char * fileName)
{
  setDefaults();
  readSettingsFile(fileName);
  //if (neighbourListBuffer == 0)
  //skipNeighbourUpdate = 1;
  //else
  //skipNeighbourUpdate = (int) (neighbourListBuffer * bondLength / MAX_FORCE);
  nd = new NormalDistribution(NormalDistribution::GSL, seed);
  aaProperties = new AA(AAFile, avgBondLength);
  sqrt_dt = sqrt(WorldSettings::dt);
  maxFENELengthSquared = WorldSettings::feneLength 
    * WorldSettings::feneLength 
    * WorldSettings::bondLength
    * WorldSettings::bondLength;
  saveSettings(fileName); // to outFile.conf
}

void WorldSettings::setDefaults()
{
  //seed with current computer clock time in nanoseconds
  timespec tp;
  clock_gettime(CLOCK_REALTIME, &tp);
  seed = tp.tv_nsec;
  
  outFile = "outfile";
  stepsBetweenFileOutput = 10000;
  stepsBetweenSnapshot = 10000;
  stepsBetweenSeparation = 10000;
  outSetting = SCREEN;

  snapshotFile = "outfile.xyz";
  separationFile = "outfile.sep";
  
  contactMapFile = "contactmap.txt";
  checkpointFile = "checkpoint.txt";
  checkpointCounter = 1;

  LJSize = sqrt(3)/sqrt(2);
  bondLength = sqrt(3)/sqrt(2);
  avgBondLength = 1.32;
  debyeLength = 3;
  feneLength = 3;
  coulombStrength = 1.243;
  hydrophobicStrength = 0;
  eLJ = 1;
  
  MAX_FORCE = 0.5 * bondLength;
  
  typeForceRepulsive = REPULSIVE_NORMAL;
  typeForceHydrophobic = HYDROPHOBIC_NORMAL;
  typeForceBond = BOND_FENE;
  typeForceCharge = DEBYE;
  typeForceExternal = EXT_NONE;

  typeHydrodynamics = HYDRO_NONE;

  initialCondition = IC_LINE;

  dt = 0.001;
  steps = 10000000;
  eqSteps = 0;
  skip = 0;

  offset = 1500*bondLength + 0.04;
  offsetVector = Vector3D(offset, offset, offset);
  
  nupFile = "Nups/TestNup1.txt";
  AAFile = "Nups/Normal.aa";
  
  //set flags
  terminate = false;
  saveCheckpoint = false;

  //error counts
  errorCountFENE = 0;
  errorCountLJ = 0;
  totalCountFENE = 0;
  totalCountLJ = 0;
  
}

int WorldSettings::readSettingsFile(const char * fileName)
{
  std::ifstream file;
  file.open(fileName);

  if(!file.is_open())
    {
      std::cout<<"Unable to read settings file: "<<fileName<<std::endl;
      return -1;
    }

  std::string line;

  while(std::getline(file, line))
    {
      std::istringstream ss(line);
      std::string property;

      if(std::getline(ss, property, '='))
        {
	  std::string value;
	  if(std::getline(ss, value,'#'))
            {
	      if(addProperty(property.c_str(), value.c_str()) < 0)
                {
		  std::cout<<"Invalid property/value combination: "<<property<<"="<<value<<std::endl;
                }
            }
        }
    }

  file.close();

  return 0;
}

int WorldSettings::addProperty(const char * property, const char * value)
{
  if (strcmp(property, "nupFile") == 0)
    WorldSettings::nupFile = value;
  else if (strcmp(property, "AAFile") == 0)
    WorldSettings::AAFile = value;
  else if (strcmp(property, "avgBondLength") == 0)
    WorldSettings::avgBondLength = atof(value);
  else if (strcmp(property, "skipNeighbourUpdate") == 0)
    WorldSettings::skipNeighbourUpdate = atoi(value);
  else if (strcmp(property, "hydrophobicCutoff") == 0)
    WorldSettings::hydrophobicCutoff = atof(value);
  else if (strcmp(property, "neighbourListBuffer") == 0)
    WorldSettings::neighbourListBuffer = atof(value);
  else if (strcmp(property, "debyeLength") == 0)
    WorldSettings::debyeLength = atof(value);
  else if (strcmp(property, "feneLength") == 0)
    WorldSettings::feneLength = atof(value);
  else if (strcmp(property, "offset") == 0){
    WorldSettings::offset = atof(value);
    offsetVector = Vector3D(offset, offset, offset);
  }
  else if (strcmp(property, "coulombStrength") == 0)
    WorldSettings::coulombStrength = atof(value);
  else if (strcmp(property, "hydrophobicStrength") == 0)
    WorldSettings::hydrophobicStrength = atof(value);
  else if (strcmp(property, "eLJ") == 0)
    WorldSettings::eLJ = atof(value);
  else if (strcmp(property, "graftingDistance") == 0) {}
    //graftingDistance = atof(value);
  else if (strcmp(property, "numberOfNanoparticles") == 0) {}
    //numberOfNanoparticles = atoi(value);
  else if (strcmp(property, "MAX_FORCE") == 0)
    {
      MAX_FORCE = atof(value) * bondLength;
    }
  else if (strcmp(property, "typeForceCharge") == 0)
    {
      if (strcmp(value, "TYPE_FORCE_CHARGE_COULOMB") == 0)
	{
	  typeForceCharge = COULOMB;
	  //funcForceHydrophobic = &forceLJHydrophobic;
	}
      else if (strcmp(value, "TYPE_FORCE_CHARGE_DEBYE") == 0)
	{
	  typeForceCharge = DEBYE;
	  //funcForceHydrophobic = &forceLJHydrophobicSameRange;
	}
    }
  else if (strcmp(property, "typeForceRepulsive") == 0)
    {
      if (strcmp(value, "TYPE_FORCE_REPULSIVE_NORMAL") == 0)
	{
	  typeForceRepulsive = REPULSIVE_NORMAL;
	}
      if (strcmp(value, "TYPE_FORCE_REPULSIVE_126") == 0)
	{
	  typeForceRepulsive = REPULSIVE_LJ126;
	}
    }
  else if (strcmp(property, "typeForceHydrophobic") == 0)
    {
      if (strcmp(value, "TYPE_FORCE_HYDROPHOBIC_NORMAL") == 0)
	{
	  typeForceHydrophobic = HYDROPHOBIC_NORMAL;
	}
      if (strcmp(value, "TYPE_FORCE_HYDROPHOBIC_126") == 0)
	{
	  typeForceHydrophobic = HYDROPHOBIC_LJ126;
	}
      else if (strcmp(value, "TYPE_FORCE_HYDROPHOBIC_SAME_RANGE") == 0)
	{
	  typeForceHydrophobic = HYDROPHOBIC_SAME_RANGE;
	}
    }
  else if (strcmp(property, "typeForceBond") == 0)
    {
      if (strcmp(value, "TYPE_FORCE_BOND_HOOKE") == 0)
	{
	  typeForceBond = BOND_HOOKE;
	  //funcForceBond = &forceSpringHooke;
	}
      else if (strcmp(value, "TYPE_FORCE_BOND_EXP") == 0)
	{
	  typeForceBond = BOND_EXP;
	}
      else if (strcmp(value, "TYPE_FORCE_BOND_FENE") == 0)
	{
	  typeForceBond = BOND_FENE;
	  //funcForceBond = &forceSpringFENE;
	}
    }
  else if (strcmp(property, "typeHydrodynamics") == 0)
    {
      if (strcmp(value, "TYPE_HYDRO_NONE") == 0)
	{
	  typeHydrodynamics = HYDRO_NONE;
	}
      else if (strcmp(value, "TYPE_HYDRO_RPY") == 0)
	{
	  typeHydrodynamics = HYDRO_RPY;
	}
    }

  else if (strcmp(property, "initialCondition") == 0)
    {
      if (strcmp(value, "IC_LINE") == 0)
	{
	  initialCondition = IC_LINE;
	}
      else if (strcmp(value, "IC_RW") == 0)
	{
	  initialCondition = IC_RW;
	}
      else if (strcmp(value, "IC_SAW") == 0)
	{
	  initialCondition = IC_SAW;
	}
    }
  
  else if (strcmp(property, "typeForceExternal") == 0)
    {
      if (strcmp(value, "TYPE_FORCE_EXTERNAL_BLANK") == 0)
	{
	  typeForceExternal = EXT_NONE;
	}
      else if (strcmp(value, "TYPE_FORCE_EXTERNAL_ZWALL") == 0)
	{
	  typeForceExternal = EXT_ZWALL;
	}
    }
  else if (strcmp(property, "dt") == 0)
    WorldSettings::dt = atof(value);
  else if (strcmp(property, "steps") == 0)
    WorldSettings::steps = atoi(value);
  else if (strcmp(property, "eqSteps") == 0)
    WorldSettings::eqSteps = atoi(value);
  else if (strcmp(property, "skip") == 0)
    WorldSettings::skip = atoi(value);
  else if (strcmp(property, "seed") == 0)
    seed = atoi(value);
  else if (strcmp(property, "outFile") == 0)
    outFile = value;
  else if (strcmp(property, "snapshotFile") == 0)
    snapshotFile = value;
  else if (strcmp(property, "stepsBetweenFileOutput") == 0)
    stepsBetweenFileOutput = atoi(value);
  else if (strcmp(property, "stepsBetweenSnapshot") == 0)
    stepsBetweenSnapshot = atoi(value);
  else if (strcmp(property, "stepsBetweenSeparation") == 0)
    stepsBetweenSeparation = atoi(value);
  else if (strcmp(property, "contactMapFile") == 0)
    contactMapFile = value;
  else if (strcmp(property, "checkpointFile") == 0)
    checkpointFile = value;
  else if (strcmp(property, "outSetting") == 0)
    {
      if(strcmp(value, "FILE") == 0)
	outSetting = FILE;
      else if (strcmp(value, "SCREEN")==0)
	outSetting = SCREEN;
    }
  else
    {
      std::cout<<"Unknown property: "<<property<<std::endl;
      return 0;
    }
  return 0;
}

/// \todo does not include all settings
void WorldSettings::saveSettings(const char * fileName)
{
  std::stringstream filenamestream;
  filenamestream << WorldSettings::outFile << ".conf";
  std::string filename = filenamestream.str();
  std::ofstream file;
  file.open(filename.c_str());
  std::cout<<"Saving settings to: "<<filename<<std::endl;

  file<<"#Settings loaded from: "<<fileName<<std::endl;
  file<<"outFile="<<outFile<<std::endl;
  if(outSetting == SCREEN)
    file<<"outSetting=SCREEN"<<std::endl;
  else if(outSetting == FILE)
    file<<"outSetting=FILE"<<std::endl;
  file<<"stepsBetweenFileOutput="<<stepsBetweenFileOutput<<std::endl;
  file<<"stepsBetweenSnapshot="<<stepsBetweenSnapshot<<std::endl;
  file<<"contactMapFile="<<contactMapFile<<std::endl;
  file<<"checkpointFile="<<checkpointFile<<std::endl;
  file<<"snapshotFile="<<snapshotFile<<std::endl;

  file<<"nupFile="<<nupFile<<std::endl;
  file<<"bondLength="<<bondLength<<std::endl;
  file<<"avgBondLength="<<avgBondLength<<std::endl;
  file<<"debyeLength="<<debyeLength<<std::endl;
  file<<"feneLength="<<feneLength<<std::endl;
  file<<"coulombStrength="<<coulombStrength<<std::endl;
  file<<"hydrophobicStrength="<<hydrophobicStrength<<std::endl;
  file<<"eLJ="<<eLJ<<std::endl;
  //file<<"graftingDistance="<<graftingDistance<<std::endl;
  //file<<"numberOfNanoparticles="<<numberOfNanoparticles<<std::endl;
  file<<"skipNeighbourUpdate="<<skipNeighbourUpdate<<std::endl;

  if (typeForceCharge == COULOMB)
    file<<"typeForceCharge=TYPE_FORCE_CHARGE_COULOMB"<<std::endl;
  else if (typeForceCharge == DEBYE)
    file<<"typeForceCharge=TYPE_FORCE_CHARGE_DEBYE"<<std::endl;
  file<<std::endl;
  if (typeForceHydrophobic == HYDROPHOBIC_NORMAL)
    file<<"typeForceHydrophobic=TYPE_FORCE_HYDROPHOBIC_NORMAL"<<std::endl;
  else if (typeForceHydrophobic ==HYDROPHOBIC_SAME_RANGE)
    file<<"typeForceHydrophobic=TYPE_FORCE_HYDROPHOBIC_SAME_RANGE"<<std::endl;
  file<<std::endl;
  if (typeForceBond == BOND_HOOKE)
    file<<"typeForceBond=TYPE_FORCE_BOND_HOOKE"<<std::endl;
  else if (typeForceBond == BOND_FENE)
    file<<"typeForceBond=TYPE_FORCE_BOND_FENE"<<std::endl;
  file<<std::endl;
  file<<"dt="<<dt<<std::endl;
  file<<"steps="<<steps<<std::endl;
  file<<"eqSteps="<<eqSteps<<std::endl;
  file<<"skip="<<skip<<std::endl;
  file<<"offset="<<offset<<std::endl;
  file<<"seed="<<seed<<std::endl;

}
