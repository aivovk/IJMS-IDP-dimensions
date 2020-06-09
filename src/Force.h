#ifndef FORCE_H
#define FORCE_H

#include "WorldSettings.h"
#include "Vector3D.h"

/**
  \brief Interparticle (bond, repulsive, cohesive, charged), random, and
  external forces

  the selector functions choose one of the options based on WorldSettings
  
  r is the vector between the two particles
  
  All quantities use SIMULATION units except:
  size has BL units = sqrt(1.5) * SIM units
 */

//selector
Vector3D forceBond(const Vector3D &r);
//options
Vector3D forceBondHooke(const Vector3D &r);
Vector3D forceBondFENE(const Vector3D &r);
Vector3D forceBondExp(const Vector3D &r);

//random forces
Vector3D noiseTerm(TYPE_FLOAT stddev);
/**
  Returns a random number from a normal distribution with variance equal to the
  time step
 */
Vector3D noiseTerm();

//selector
Vector3D forceRepulsive(const Vector3D &r, TYPE_FLOAT size);
//options
Vector3D forceExpRepulsive(const Vector3D &r, TYPE_FLOAT size);
Vector3D forceLJ86Repulsive(const Vector3D &r, TYPE_FLOAT size);
Vector3D forceLJ126Repulsive(const Vector3D &r, TYPE_FLOAT size);
Vector3D forceLJ96Repulsive(const Vector3D &r, TYPE_FLOAT size);

//selector
Vector3D forceCohesive(const Vector3D &r, TYPE_FLOAT size);
//options
Vector3D forceLJ86Attractive(const Vector3D &r, TYPE_FLOAT size);
Vector3D forceLJ126Attractive(const Vector3D &r, TYPE_FLOAT size);
Vector3D forceLJ86AttractiveSameRange(const Vector3D &r, TYPE_FLOAT size);

//selector
Vector3D forceCharge(const Vector3D &r, TYPE_FLOAT charge);
//options
Vector3D forceChargeCoulomb(const Vector3D &r, TYPE_FLOAT charge);
Vector3D forceChargeDebye(const Vector3D &r, TYPE_FLOAT charge);

//external forces
Vector3D forceExternal(const Vector3D &r, TYPE_FLOAT size);
Vector3D forceExternalZWall(const Vector3D &r, TYPE_FLOAT size);

#endif // FORCE_H
