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
  CubeSpace();
  CubeSpace(int length, std::vector<Particle *> * particles, TYPE_FLOAT cubeSize_);
  ~CubeSpace();
  
  /// place all particles into the cubes
  void update();

  /// remove all particles from the cubes
  void clear();
  
  /// compute a list of indices of neighbouring particles for each particle
  std::vector< std::vector<int> > * getNeighbours();
  
  
 protected:
 private:

  static const int stencil[14][3];
  
  // helper variables used in getNeighbours
  std::vector<int>::iterator iter;
  std::vector<int> * cube;
  int ix, iy, iz;
  int bx, by, bz;

  /// temporary pointer to list of particles in a cube used in clear and update
  std::vector<int> ** cube_;

  int noOfParticles;
  int length; ///< number of smaller cubes in 1 dimension
  int xLength;
  int yLength;
  int zLength;
  TYPE_FLOAT cubeSize; ///<

  /// list of particle coordinates
  std::vector<Particle *> * particles;

  /// array of cubes containing lists of particles within each
  std::vector< std::vector< std::vector< std::vector<int> * > > > cubes;

  /// neighbour list for each particle
  std::vector< std::vector< int > > neighbours;

  /// get list of particles within a cube pased on position
  std::vector<int> ** getCube(Vector3D rParticle);

  /// get list of particles within a cube based on a single particle index
  std::vector<int> ** getCube(int iParticle);
};

#endif // CUBESPACE_H
