#ifndef WORLDSETTINGS_H
#define WORLDSETTINGS_H

#include <iostream>
#include <cstring>
#include <stdio.h>
#include <fstream>
#include <string>
#include <sstream>

#include "NormalDistribution.h"
#include "PeriodicBoundary.h"
#include "Vector3D.h"
#include "AA.h"

/**
 * WorldSettings Class
 *
 * \brief Store global properties of a simulation run.
 *
 * Store global properties that remain unchanged during a simulation
 * run. Settings are loaded from a specified file at the start of a simulation.
 * Default settings file: "./WorldSettings.conf"
 *
 * Each line must contain the following syntax:
 *
 * PROPERTY=VALUE
 * PROPERTY=VALUE #Comment
 * #Comments with = signs will be interpreted of as commands but not recognized
 * #since no PROPERTY names have # signs
 *
 * If a PROPERTY name is unrecognized or its VALUE is invalid a message will be
 * printed and the PROPERTY will be set to its default value.
 *
 * PROPERTY             DEFAULT DESCRIPTION
 *
 * nupFile                      file containing sequence of 1-letter AA names
 *                              must be in Nups directory
 *
 * AAFile                       file containing AA (Particle) properties for the
 *                              sequence given in nupFile
 *
 * skipNeighbourUpdate  1       number of steps between neighbour array updates
 *
 * bondLength           1.22... expected Rouse model bond length sqrt(1.5)
 *                              used as a unit of length
 *
 * debyeLength          3       [in simulation units]
 * 
 * feneLength           3       l_max of FENE force [in units of bondLength]
 *
 * offset               1500    all x,y,z coordinates are offset by this amount
 *                              must be large enough so that coordinates cannot
 *                              reach a negative value
 *                              [in units of bondLength]
 *
 * coulombStrength      1.243   constant multiplier of coulomb force term
 *                              depends on actual bond length within the
 *                              simulation:
 *                              =(4.124/3.792)*(actual bond length/bondLength)
 *                              =(4.124/3.792)*(1.4/sqrt(1.5))
 *
 * hydrophobicStrength  1.6     constant multiplier of hydrophobic force term
 *
 * dt                   0.001   time step
 *
 * steps                10^7    number of steps the simulation will run for
 *
 * eqSteps              0       number of steps to run for before turning on
 *                              HI and hydrophobic forces
 *
 * outSetting                   FILE or SCREEN
 *
 * outFile                      filename to write simulation results.
 *
 * OTHER        DESCRIPTION
 *
 * sqrt_dt      square root of dt
 *
 * nd           pointer to the NormalDistribution object for RNG
 *
 * terminate    signals the program to save the state (checkpoint)
 *
 */
class WorldSettings
{
 public:
  WorldSettings();
  ~WorldSettings();

  static void initialize();
  static void initialize(const char * fileName);

  // load the settings
  static int readSettingsFile(const char * fileName);
  
  // settings are rewritten to outFile.conf in case the configuration file
  // fileName (usually defaultSettingsFile) is incomplete or is changed in the
  // future
  static void saveSettings(const char * fileName);

  static const char * defaultSettingsFile;

  // file describing AA (Particle) properties that will be used
  static std::string AAFile;

  // file containing monomer sequence
  static std::string nupFile;

  enum TYPE_OUTPUT {FILE, SCREEN};
  static TYPE_OUTPUT outSetting;

  enum TYPE_FORCE_CHARGE {COULOMB, DEBYE};
  static TYPE_FORCE_CHARGE typeForceCharge;
  
  enum TYPE_FORCE_REPULSIVE {REPULSIVE_NORMAL, REPULSIVE_LJ126, REPULSIVE_EXP};
  static TYPE_FORCE_REPULSIVE typeForceRepulsive;

  enum TYPE_FORCE_HYDROPHOBIC {HYDROPHOBIC_NORMAL,
			       HYDROPHOBIC_LJ126,
			       HYDROPHOBIC_SAME_RANGE};
  static TYPE_FORCE_HYDROPHOBIC typeForceHydrophobic;

  enum TYPE_FORCE_BOND {BOND_HOOKE, BOND_FENE, BOND_EXP};
  static TYPE_FORCE_BOND typeForceBond;

  
  enum TYPE_HYDRO {HYDRO_NONE, HYDRO_RPY};
  static TYPE_HYDRO typeHydrodynamics;
  
  enum TYPE_FORCE_EXTERNAL {EXT_ZWALL, EXT_NONE};
  static TYPE_FORCE_EXTERNAL typeForceExternal;

  enum TYPE_INITIAL_CONDITION {IC_LINE, IC_RW, IC_SAW};
  static TYPE_INITIAL_CONDITION initialCondition;
  
  // output file to save results
  static std::string outFile;
  static int stepsBetweenFileOutput;

  // output file to save snapshots
  static std::string snapshotFile;
  static int stepsBetweenSnapshot;

  // output file for interparticle distances
  static std::string separationFile;
  static int stepsBetweenSeparation;
	
  static std::string contactMapFile;

  /// checkpoint filename to save/load from
  static std::string checkpointFile;
  static int checkpointCounter;
  
  static NormalDistribution * normalDistribution;
  static PeriodicBoundary pbc;
  static unsigned int seed;
  static AA * aaProperties;

  // simulation details
  static TYPE_FLOAT LJSize; ///< \todo explain
  static TYPE_FLOAT bondLength; ///< $\sqrt{1.5}$
  static TYPE_FLOAT avgBondLength; ///< \todo explan
  static TYPE_FLOAT debyeLength; ///< \todo Units?
  static TYPE_FLOAT feneLength; ///< max bond length \todo Units?
  static TYPE_FLOAT maxFENELengthSquared; ///< \todo explain
  static TYPE_FLOAT coulombStrength; ///< strength of charged forces \todo Units?
  static TYPE_FLOAT hydrophobicStrength; ///< strength of cohesive forces [kT]
  static TYPE_FLOAT eLJ; ///< strength of repulsive forces [kT]

  static TYPE_FLOAT hydrophobicCutoff; ///< \todo units
  static TYPE_FLOAT neighbourListBuffer; ///< \todo units

  /// multiple of steps to perform neighbour list updates, should be
  /// neighbourListBuffer/MAX_FORCE \todo Units
  static int skipNeighbourUpdate;

  /// time step
  static TYPE_FLOAT dt;
  static TYPE_FLOAT sqrt_dt;
  /// total simulation steps to run
  static int steps;   
  /// number of steps to run without hydrodynamic interactions and hydrophobic
  /// forces
  static int eqSteps;
  /// number of steps to run before file IO? \todo not used?
  static int skip;
  /// amount to offset initial coordinates by (cannot be <0)
  static TYPE_FLOAT offset;
  static Vector3D offsetVector;

  // flags
  static bool terminate;
  static bool saveCheckpoint;

  // counters if Particle positions try to change more than MAX
  static unsigned int errorCountFENE;
  static unsigned int errorCountLJ;
  static unsigned int totalCountFENE;
  static unsigned int totalCountLJ;

  /// error maxes \todo units
  static TYPE_FLOAT MAX_FORCE;

  /// \todo explain
  static int noOfMonomers;

 protected:
 private:
  static void setDefaults();
  static int addProperty(const char * property, const char * value);
};

#endif // WORLDSETTINGS_H
