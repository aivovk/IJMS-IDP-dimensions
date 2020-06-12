#include "CubeSpace.h"

// coordinates (relative to current cube) of cubes to check for neighbours
// only need to check half of the 26 neighbours and the current cube
const int CubeSpace::stencil[14][3] = {{0,0,0},
				       {1,0,0},
				       {-1,1,0},
				       {0,1,0},
				       {1,1,0},
				       {-1,-1,1},
				       {0,-1,1},
				       {1,-1,1},
				       {-1,0,1},
				       {0,0,1},
				       {1,0,1},
				       {-1,1,1},
				       {0,1,1},
				       {1,1,1}};

CubeSpace::CubeSpace(int length, TYPE_FLOAT cubeSize, int noOfParticles)
  : xLength(length),
    yLength(length),
    zLength(length),
    cubeSize(cubeSize),
    noOfParticles(noOfParticles)
{    
    /// \todo not expected to work with periodic boundary conditions
    if(WorldSettings::pbc.isPeriodicX)
      xLength = WorldSettings::pbc.xSize;
    if(WorldSettings::pbc.isPeriodicY)
      yLength = WorldSettings::pbc.ySize;
    if(WorldSettings::pbc.isPeriodicZ)
      zLength = WorldSettings::pbc.zSize;

    neighbours.resize(noOfParticles);

    cubes.resize(xLength);
    for (int i = 0; i < xLength; i++)
      {
        cubes[i].resize(yLength);
	
        for (int j = 0; j < yLength; j++)
        {
            cubes[i][j].resize(zLength);
        }
    }
}

// remove all particles from the cubes
void CubeSpace::clear(const std::vector<Particle *>& particles)
{
    for (int i = 0; i < noOfParticles; i++)
    {
      getCube(particles[i]->r).clear();
    }
}

// places particles into the cubes
void CubeSpace::update(const std::vector<Particle *>& particles)
{
    for (int i = 0 ; i < noOfParticles ; i++)
    {
      getCube(particles[i]->r).push_back(i);
    }
}

// compute a list of neighbouring particles for each particle
std::vector< std::set< int > > * CubeSpace::getNeighbours(const std::vector<Particle *>& particles)
{
  int ix, iy, iz, bx, by, bz;
  for (int iParticle = 0 ; iParticle < noOfParticles ; iParticle++)
    {
      neighbours[iParticle].clear();
      
      //loop through surrounding cubes
      ix = particles[iParticle]->r.x/cubeSize;
      iy = particles[iParticle]->r.y/cubeSize;
      iz = particles[iParticle]->r.z/cubeSize;
      
      for (int iStencil = 0; iStencil < 14; iStencil++)
	{
	  bx = ix + stencil[iStencil][0];
	  by = iy + stencil[iStencil][1];
	  bz = iz + stencil[iStencil][2];
	  
	  std::vector<int>& cube = cubes[bx%xLength][by%yLength][bz%zLength];
	  for (std::vector<int>::iterator iter = cube.begin(); iter!= cube.end(); iter++)
	    {
	      // for particles in the same cube
	      // only particles with higher indices are included
	      // for particles in other cubes, all indices are included
	      // since only checking half of neighbouring cubes
	      /// \todo why do duplicates occur?
	      if (iStencil != 0 || *iter > iParticle)
		{
		  neighbours[iParticle].insert(*iter);
		}
	    }
	}
    }
  return &neighbours;
}

std::vector<int>& CubeSpace::getCube(const Vector3D& rParticle)
{
  int ix = rParticle.x/cubeSize;
  ix %= xLength;
  int iy = rParticle.y/cubeSize;
  iy %= yLength;
  int iz = rParticle.z/cubeSize;
  iz %= zLength;
  return cubes[ix][iy][iz];
}
