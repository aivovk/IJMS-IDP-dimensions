#include "Particle.h"

Particle::Particle()
{
  r = Vector3D();
  dr = Vector3D();
  type = BLANK_PARTICLE;
  typeChar = '0';
  radiusH = 1;
  radiusLJ = 1;
  nextParticleOnChain = NULL;
  charge = 0;
  hydrophobicity = 0;
  fixed = false;
}

Particle::Particle(Vector3D r_,
		   int type_,
		   char typeChar_,
		   Particle * nextParticleOnChain_,
		   bool fixed_){
  r = r_;
  dr = Vector3D();
  type = type_;
  typeChar = typeChar_;
  radiusH = WorldSettings::aaProperties->getHRadius(typeChar);
  radiusLJ = WorldSettings::aaProperties->getLJRadius(typeChar);
  nextParticleOnChain = nextParticleOnChain_;
  charge = WorldSettings::aaProperties->getCharge(typeChar);
  hydrophobicity = WorldSettings::aaProperties->getHydrophobicity(typeChar);
  fixed = fixed_;
}

