#ifndef AA_H
#define AA_H

#include <map>
#include <vector>
#include <string>
#include <iostream>
#include <cstring>
#include <stdio.h>
#include <fstream>
#include <string>
#include <sstream>
#include <stdlib.h> 
#include "Vector3D.h"

/**
  AA Class
  \brief Describes the available particles and their properties.

  Read in Amino Acid (or general Particle) properties, and convert to simulation
  units, in order to be used by the Particle class.

  Notes:
  nonpolar refers to hydrophobic/cohesive
 */
class AA{
 public:
  /**
    \param AAFile - (CSV File)
    \param abl - average bond length in simulation units

    AAFile structure (entries are INT, CHAR, FLOAT, or STR):
    START
    numberOfAA,INT
    minAAsize,FLOAT
    positive,STR
    negative,STR
    nonpolar,STR
    aaList,aaNames,aaSetBradii,aaHydrophobicStrength,aaSetHdiameter,aaCharge
    CHAR,STRING,FLOAT,FLOAT,FLOAT,FLOAT
    ...
    ... (numberOfAA entries)
    END

    The list of headings is not entirely correct:
    Columns 1, 2 are particle names in single character and full name
    Column 3 is the LJ diameter in real units
    Column 4 is the hydrophobic strength in simulation units [kT]
    Column 5 is the hydrodynamic diameter in real units
    Column 6 is the charge in units of elementary charge

    Column 1 are the particle names that are used to specify a monomer sequence

    Distances from the AAFile, are converted into simulations units [S] so that
    that the bond length in real units [U] is realBondLength (= 3.792):
    r [S] = r [U] * abl [S] / 3.792 [U]
   */
  AA(std::string AAFile, TYPE_FLOAT abl);

  /// 3.792 [Angstroms]
  static const TYPE_FLOAT realBondLength; //Angstrom Units

  /// correction factor to LJ sizes for testing
  static const TYPE_FLOAT LJFactor; 

  /// number of different particles defined in AA File
  int numberOfAA;
  
  /// array of 1 letter particle names
  std::vector<char> aaList; ///< the property arrays are ordered in the same way
  
  TYPE_FLOAT getHRadius(char aaCode); ///< hydrodynamic radius
  TYPE_FLOAT getLJRadius(char aaCode);
  TYPE_FLOAT getCharge(char aaCode);
  TYPE_FLOAT getHydrophobicity(char aaCode);
  TYPE_FLOAT maxLJRadius(); ///< LJ radius of the largest particle

  //char to size map
  /// LJ diameters
  std::map<char, TYPE_FLOAT> aaLJDiameters; ///< SIM units
  /// Hydrodynamic diameters
  std::map<char, TYPE_FLOAT> aaHDiameters; ///< SIM units

  //int to size map
  /// LJ diameters
  std::vector<TYPE_FLOAT> aaLJDiameter; ///< SIM units
  /// Hydrodynamic diameters
  std::vector<TYPE_FLOAT> aaHDiameter; ///< SIM units

  // used to construct the size arrays
  TYPE_FLOAT minAASize;
  /// LJ diameters
  std::vector<TYPE_FLOAT> aaSetBdiameter; ///< Ang units
  /// Hydrodynamic diameters
  std::vector<TYPE_FLOAT> aaSetHdiameter; ///< Ang units

  /// strengths used for hydrophobic/cohesive attraction
  std::map<char, TYPE_FLOAT> aaHStrengths;
  std::vector<TYPE_FLOAT> aaHydrophobicStrength;

  //charge
  std::map<char, TYPE_FLOAT> aaChargeStrengths;
  std::vector<TYPE_FLOAT> aaCharge;

  /// string particle names
  std::map<char, std::string> aaMap;
  std::vector<std::string> aaNames;

  /** list of 1 character names that have the given property */
  std::string positive;
  std::string negative;
  std::string nonpolar; ///< hydrophobic/cohesive
 
 private:
  TYPE_FLOAT avgBondLength; //SIM Units
};

#endif
