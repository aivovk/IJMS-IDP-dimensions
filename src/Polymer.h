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
     \return true if the particles don't overlap
  */
  bool isSelfAvoiding(Vector3D r) const;

  /// calculate the change in position due to bond forces
  void simulate(TYPE_FLOAT t);
  
  /// vector from first to last monomer
  Vector3D endToEndVector() const;

  TYPE_FLOAT endToEndDistanceSquared() const;


  /// average of monomer positions (assuming same masses)
  Vector3D centreOfMass() const;

  TYPE_FLOAT radiusOfGyrationSquared() const;

  TYPE_FLOAT averageBondLengthSquared() const;

  /// draw cylinders lines between adjacent monomers
  void draw() const;

 protected:
 private:
  int noOfMonomers;

  // positions of the monomers
  std::vector<Particle> monomers;
};

#endif // POLYMER_H
