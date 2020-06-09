#ifndef POLYMER_H
#define POLYMER_H

//#include <cmath>
#include <iostream>

#include <GL/gl.h>
#include <GL/glu.h>
//#include <omp.h>

#include "NormalDistribution.h"
#include "CubeSpace.h"
#include "Vector3D.h"
#include "Force.h"

/**
   Polymer Class
   \brief Polymer is an array of Particle objects connected by spring forces
 */
class Polymer
{
 public:

  Polymer();


  /**
     Polymer constructor
     \param noOfMonomers
     \param sequence list of 1 char AA (Particle type) names
     \param random if initial condition is a random walk \todo not used
     \param fixedStart if the first monomer (i==0) is fixed (constant position)
     \param start position of the first monomer (i==0)
     \param orientation direction of monomers (for non-random initial condition)
     \param particles pointer to array of monomer positions (to fill)
   */
  Polymer(int noOfMonomers,
	  std::string sequence,
	  bool random,
	  bool fixedStart,
	  Vector3D start,
	  Vector3D orientation,
	  std::vector<Particle *> * particles);

  ~Polymer();

  int getSize(){return noOfMonomers;};

  /**
     self-avoiding check (used to determine initial conditions)
     Note: uses bond length distance instead of particle diameter
     \param r distance between two monomers
     \return true if the particles overlap
  */
  bool distanceCheck(Vector3D r);

  /// calculate the change in position due to bond forces
  void simulate(TYPE_FLOAT t);
  
  /// vector from first to last monomer
  Vector3D endToEndVector();

  TYPE_FLOAT endToEndDistanceSquared();


  /// average of monomer positions (assuming same masses)
  Vector3D centreOfMass();

  TYPE_FLOAT radiusOfGyrationSquared();

  TYPE_FLOAT averageBondLengthSquared();

  /// draw white lines between adjacent monomers
  void draw(TYPE_FLOAT scale);

 protected:
 private:
  int noOfMonomers;

  // positions of the monomers
  std::vector<Particle> chain;
};

#endif // POLYMER_H
