#ifndef CUBESPACE_H
#define CUBESPACE_H

#include <iostream>
#include <vector>
#include <set>

#include "Particle.h"
#include "PeriodicBoundary.h"
#include "WorldSettings.h"

/**
   Class CubeSpace
   \brief (AKA Cell List)
 
   Represents infinite space using a periodic box (cube) with sidelength
   length*cubeSize.  The cube is subdivided into length^3 smaller cubes each
   with sidelength cubeSize.

   \member particles points to the list of particle coordinates in the World
   class. Each particle is placed into its corresponding cube based on it's
   coordinates. This way all particles which may possibly interact with it using
   a cutoff LJ potential can easily be found by computing the indices of the
   surrounding cubes (complexity O(1)).
 
   Note: cubeSize must be greater or equal to the maximum LJ
   interaction distance (cutoff)

   \todo Check Units
 */
class CubeSpace
{
 public:
  CubeSpace(){};
  CubeSpace(int length, TYPE_FLOAT cubeSize, int noOfParticles);
  
  /// place all particles into the cubes
  void update(const std::vector<Particle *>& particles);

  /// remove all particles from the cubes, must be called before particle
  /// positions are changed and next update
  void clear(const std::vector<Particle *>& particles);
  
  /// compute a list of indices of neighbouring particles for each particle
  std::vector< std::set<int> > * getNeighbours(const std::vector<Particle *>& particles);
  
  
 protected:
 private:
  static const int stencil[14][3];

  int xLength;
  int yLength;
  int zLength;
  TYPE_FLOAT cubeSize; ///<
  int noOfParticles;

  /// array of cubes containing lists of particles within each
  std::vector< std::vector< std::vector< std::vector<int> > > > cubes;

  /// neighbour list for each particle
  std::vector< std::set< int > > neighbours;

  /// get list of particles within a cube pased on position
  std::vector<int>& getCube(const Vector3D& rParticle);

};

#endif // CUBESPACE_H
