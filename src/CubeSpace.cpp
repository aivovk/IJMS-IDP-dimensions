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

CubeSpace::CubeSpace()
{

}

CubeSpace::CubeSpace(int length, std::vector<Particle *> * particles, TYPE_FLOAT cubeSize_)
{

    this->length = length;
    cubeSize = cubeSize_;
    this->particles = particles;

    xLength = length;
    yLength = length;
    zLength = length;
    
    /// \todo not expected to work with periodic boundary conditions
    if(WorldSettings::pbc.isPeriodicX)
      xLength = WorldSettings::pbc.xSize;
    if(WorldSettings::pbc.isPeriodicY)
      yLength = WorldSettings::pbc.ySize;
    if(WorldSettings::pbc.isPeriodicZ)
      zLength = WorldSettings::pbc.zSize;

    noOfParticles = particles->size();
    neighbours.resize(noOfParticles);

    // allocate memory for the cubes
    cubes.resize(xLength);
    for (int i = 0; i < xLength; i++)
    {
        cubes[i].resize(yLength);

        for (int j = 0; j < yLength; j++)
        {
            cubes[i][j].resize(zLength);
            for (int k = 0; k < zLength; k++)
            {
                cubes[i][j][k] = NULL;
            }
        }
    }
}

CubeSpace::~CubeSpace()
{
    for(int i = 0; i < xLength; i++)
        for(int j = 0; j < yLength; j++)
            for(int k = 0; k < zLength; k++)
            {
                if(cubes[i][j][k] != NULL)
                    delete cubes[i][j][k];
            }
    cubes.clear();
}

// remove all particles from the cubes
void CubeSpace::clear()
{
    for (int i = 0; i < noOfParticles; i++)
    {
        cube_ = getCube(i);
        if(*cube_ != NULL)
        {
	  (*cube_)->clear();
        }

    }
}

// places particles into the cubes
void CubeSpace::update()
{
    for (int i = 0 ; i < noOfParticles ; i++)
    {
        cube_ = getCube(i);
	
	// create a new set if the cube is empty
        if((*cube_) == NULL)
        {
	  (*cube_) = new std::vector<int>;
        }
	
	// add the particle to the set
        (*cube_)->push_back(i);
    }
}

// compute a list of neighbouring particles for each particle
std::vector< std::vector< int > > * CubeSpace::getNeighbours()
{
  for (int iParticle = 0 ; iParticle < noOfParticles ; iParticle++)
    {
      neighbours[iParticle].clear();
      
      //loop through surrounding cubes
      ix = (*particles)[iParticle]->r.x/cubeSize;
      iy = (*particles)[iParticle]->r.y/cubeSize;
      iz = (*particles)[iParticle]->r.z/cubeSize;
      
      for (int iStencil = 0; iStencil < 14; iStencil++)
	{
	  bx = ix + stencil[iStencil][0];
	  by = iy + stencil[iStencil][1];
	  bz = iz + stencil[iStencil][2];
	  
	  cube = cubes[bx%xLength][by%yLength][bz%zLength];
	  if (cube != NULL)
	    for (iter = cube->begin(); iter!= cube->end(); iter++)
	      {
		// for particles in the same cube
		if (iStencil == 0)
		  {
		    // only particles with higher indices are included
		    if(*iter > iParticle)
		      {
			// check neighbour list for duplicates
			/// \todo replace with set?
			bool found = false;
			for (unsigned int i = 0 ; !found && i < neighbours[iParticle].size() ; i++)
			  if (neighbours[iParticle][i] == *iter)
			    found = true;
			if (!found)
			  neighbours[iParticle].push_back(*iter);
		      }
		  }
		else // for particles in other cubes, all indices are included
		     // since only checking half of neighbouring cubes
		  { /// \todo why do duplicates occur?
		    bool found = false;
		    for (unsigned int i = 0 ; !found && i < neighbours[iParticle].size() ; i++)
		      if (neighbours[iParticle][i] == *iter)
			found = true;
		    if (!found)
		      neighbours[iParticle].push_back(*iter);
		  }
	      }
	}
    }
  return &neighbours;
}

std::vector<int> ** CubeSpace::getCube(int iParticle)
{
  ix = (*particles)[iParticle]->r.x/cubeSize;
  ix %= xLength;
  iy = (*particles)[iParticle]->r.y/cubeSize;
  iy %= yLength;
  iz = (*particles)[iParticle]->r.z/cubeSize;
  iz %= zLength;
  return &cubes[ix][iy][iz];
}

std::vector<int> ** CubeSpace::getCube(Vector3D rParticle)
{
  ix = rParticle.x/cubeSize;
  ix %= xLength;
  iy = rParticle.y/cubeSize;
  iy %= yLength;
  iz = rParticle.z/cubeSize;
  iz %= zLength;
  return &cubes[ix][iy][iz];
}
