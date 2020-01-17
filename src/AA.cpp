#include "AA.h"

TYPE_FLOAT const AA::realBondLength = 3.792;
TYPE_FLOAT const AA::LJFactor = 1.0;

// characters not used: B, J, O, U, X, Z
AA::AA(std::string AAFile, TYPE_FLOAT abl){
  avgBondLength = abl;
  std::cout<<"Loading Amino Acid properties: "<<AAFile<<std::endl;
  std::ifstream file;
  file.open(AAFile.c_str(), std::ios::in);
  if(file.is_open())
    {
      // first line defines number of particles
      std::string line;
      std::getline(file, line);
      std::istringstream fss(line);
      std::string val;
      std::getline(fss, val, ',');
      std::getline(fss, val, ',');
      numberOfAA = atoi(val.c_str());
      
      // second line defines the smallest LJ diameter in real units
      std::getline(file, line);
      fss.str(line);
      std::getline(fss, val, ',');
      std::getline(fss, val, ',');
      minAASize = atof(val.c_str());

      // list of positively charged particles
      std::getline(file, line);
      fss.str(line);
      std::getline(fss, val, ',');
      std::getline(fss, val, ',');
      positive = val;

      // list of negatively charged particles
      std::getline(file, line);
      fss.str(line);
      std::getline(fss, val, ',');
      std::getline(fss, val, ',');
      negative = val;

      // list of nonpolar/hydrophobic/cohesive particles
      std::getline(file, line);
      fss.str(line);
      std::getline(fss, val, ',');
      std::getline(fss, val, ',');
      nonpolar = val;

      // column headings
      std::getline(file, line);

      // 
      std::cout<<"Number of AA: "<<numberOfAA<<std::endl;
      std::cout<<"Positive: "<<positive<<std::endl;
      std::cout<<"Negative: "<<negative<<std::endl;
      std::cout<<"Hydrophobic: "<<nonpolar<<std::endl;

      // one line per particle
      for (int i = 0; i < numberOfAA; i++)
	{
	  // 1 char name
	  std::getline(file, line);
	  std::istringstream fss2(line);
	  std::getline(fss2, val, ',');
	  aaList.push_back(val[0]);

	  // string name
	  std::getline(fss2, val, ',');
	  aaNames.push_back(val);

	  // LJ diameter
	  std::getline(fss2, val, ',');
	  aaSetBdiameter.push_back(atof(val.c_str()));

	  // hydrophobic/cohesive strength
	  std::getline(fss2, val, ',');
	  aaHydrophobicStrength.push_back(atof(val.c_str()));

	  // hydrodynamic diameter
	  std::getline(fss2, val, ',');
	  aaSetHdiameter.push_back(atof(val.c_str()));

	  // charge
	  std::getline(fss2, val, ',');
	  aaCharge.push_back(atof(val.c_str()));
	  
	  // convert hydrodynamic and LJ diameters to simulation units
	  // and fill property arrays and maps
	  aaHDiameters[aaList[i]] = aaSetHdiameter[i] * avgBondLength / realBondLength;
	  aaHDiameter.push_back(aaSetHdiameter[i] * avgBondLength / realBondLength);
	  aaLJDiameters[aaList[i]] = LJFactor * aaSetBdiameter[i] * avgBondLength / realBondLength;
	  aaLJDiameter.push_back(LJFactor * aaSetBdiameter[i] * avgBondLength / realBondLength);

	  aaHStrengths[aaList[i]] = aaHydrophobicStrength[i];
	  aaChargeStrengths[aaList[i]] = aaCharge[i];
	  
	  aaMap[aaList[i]]=aaNames[i];
	  
	  std::cout<<"AA: "<<aaList[i]<<" ";
	  std::cout<<"Name: "<<aaNames[i]<<" ";
	  std::cout<<"H Size: "<<aaHDiameter[i]<<" ";
	  std::cout<<"LJ Size: "<<aaLJDiameter[i]<<" ";
	  std::cout<<"HStrength: "<<aaHydrophobicStrength[i]<<" ";
	  std::cout<<"Charge: "<<aaCharge[i]<<std::endl;
	}
      
      file.close();
      
      
    }
  else //default case
    {
      std::cout<<"Failed to open: "<<AAFile<<std::endl;
      
    }
}
  
TYPE_FLOAT AA::getLJRadius(char aaCode){
  return 0.5 * aaLJDiameters[aaCode];
}

TYPE_FLOAT AA::getHRadius(char aaCode){
  return 0.5 * aaHDiameters[aaCode];
}

/// \todo explain why? note: needs to be larger than 0.5*WorldSettings::LJSize
TYPE_FLOAT AA::maxLJRadius(){
  TYPE_FLOAT maxLJRad = 0.0;
  for (int i = 0 ; i < numberOfAA ; i ++)
    if (aaLJDiameter[i] > maxLJRad)
      maxLJRad = 0.5 * aaLJDiameter[i];
  return maxLJRad;
}
TYPE_FLOAT AA::getCharge(char aaCode){
  if(positive.find(aaCode, 0) != std::string::npos)
    {
      return aaChargeStrengths[aaCode];
    }
  if(negative.find(aaCode, 0) != std::string::npos)
    {
      return aaChargeStrengths[aaCode];
    }
  return 0;
}
TYPE_FLOAT AA::getHydrophobicity(char aaCode){
  if(nonpolar.find(aaCode, 0) != std::string::npos)
    {
      return aaHStrengths[aaCode];
    }
  return 0;
}
