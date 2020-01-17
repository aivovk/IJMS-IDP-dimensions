#ifndef FORCE_H
#define FORCE_H

#include "WorldSettings.h"
#include "Vector3D.h"

/**
  \brief Interparticle (bond, repulsive, cohesive, charged), random, and
  external forces

  the selector functions chooses one of the options based on WorldSettings
  
  \todo size has BL units = sqrt(1.5) * SIM units ?
 */

//selector
Vector3D forceBond(const Vector3D &r);
//options
Vector3D forceSpringHooke(const Vector3D &r);
Vector3D forceSpringFENE(const Vector3D &r);
Vector3D forceSpringExp(const Vector3D &r);

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
Vector3D forceLennardJones(const Vector3D &r, TYPE_FLOAT size);
Vector3D forceLJRepulsive(const Vector3D &r, TYPE_FLOAT size);
Vector3D forceLJ126Repulsive(const Vector3D &r, TYPE_FLOAT size);
Vector3D forceLJ96Repulsive(const Vector3D &r, TYPE_FLOAT size);

//selector
Vector3D forceHydrophobic(const Vector3D &r, TYPE_FLOAT size);
//options
Vector3D forceLJHydrophobic(const Vector3D &r, TYPE_FLOAT size);
Vector3D forceLJ126Hydrophobic(const Vector3D &r, TYPE_FLOAT size);
Vector3D forceLJHydrophobicSameRange(const Vector3D &r, TYPE_FLOAT size);

//selector
Vector3D forceCharge(const Vector3D &r, TYPE_FLOAT charge);
//options
Vector3D forceCoulomb(const Vector3D &r, TYPE_FLOAT charge);
Vector3D forceDebye(const Vector3D &r, TYPE_FLOAT charge);

//external forces
Vector3D forceExternal(const Vector3D &r, TYPE_FLOAT size);
Vector3D forceExternalZWall(const Vector3D &r, TYPE_FLOAT size);

#endif // FORCE_H
