#ifndef PARTICLE_H
#define PARTICLE_H

#include "Vector3D.h"
#include "WorldSettings.h"

#define BLANK_PARTICLE -1
#define AMINO_ACID 1
#define NANO_PARTICLE 2

/**
  Particle Class
  \brief A particle's properties and position.
 */
class Particle{
 public:
  Particle();
  Particle(Vector3D r_,
	   int type_,
	   char typeChar_,
	   Particle * nextParticleOnChain_,
	   bool fixed_);
  
  Vector3D r; ///< position
  Vector3D dr; ///< change in position in the current time step

  inline TYPE_FLOAT getHRadius()
  {
    return radiusH;
  };
  inline TYPE_FLOAT getLJRadius()
  {
    return radiusLJ;
  };    
  inline int getType()
  {
    return type;
  };
  inline char getAACode()
  {
    return typeChar;
  };
  inline Particle * getNextParticleOnChain()
  {
    return nextParticleOnChain;
  };
  inline TYPE_FLOAT getCharge()
  {
    return charge;
  };
  inline TYPE_FLOAT getHydrophobicity()
  {
    return hydrophobicity;
  };
  inline bool isFixed()
  {
    return fixed;
  };
  
 private:
  int type; ///< \todo not used?
  char typeChar; ///< 1 char code to reference properties in AA class
  TYPE_FLOAT radiusH; ///< hydrodynamic radius
  TYPE_FLOAT radiusLJ; ///< repulsive/cohesive LJ radius

  Particle * nextParticleOnChain; ///< \todo not used?

  TYPE_FLOAT charge;
  TYPE_FLOAT hydrophobicity;
  bool fixed; ///< fixed particles do not move
};


#endif
