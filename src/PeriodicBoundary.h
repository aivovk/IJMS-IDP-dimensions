#ifndef PERIODICBOUNDARY_H
#define PERIODICBOUNDARY_H

#include <iostream>
#include "Vector3D.h"

/**
   \todo 2020/01/14 Hasn't been updated so this class only works with isPeriodic
   set to false
 */
class PeriodicBoundary
{
 public:
  PeriodicBoundary();
  PeriodicBoundary(bool _isPeriodicX,
		   bool _isPeriodicY,
		   bool _isPeriodicZ,
		   TYPE_FLOAT _xMin,
		   TYPE_FLOAT _xMax,
		   TYPE_FLOAT _yMin,
		   TYPE_FLOAT _yMax,
		   TYPE_FLOAT _zMin,
		   TYPE_FLOAT _zMax);

  ~PeriodicBoundary();

  /*
    Convert position
   */
  inline Vector3D convertR(Vector3D r, TYPE_FLOAT bl)
  {
    if(isPeriodicX)
      {
	while ( r.x > xSize * bl / 2 )
	  r.x -= xSize * bl;
	while ( r.x < -xSize * bl / 2)
	  r.x += xSize * bl;
      }
    if(isPeriodicY)
      {
	while ( r.y > ySize * bl / 2 )
	  r.y -= ySize * bl;
	while ( r.y < -ySize * bl / 2)
	  r.y += ySize * bl;
      }
    if(isPeriodicZ)
    {
      while ( r.z > zSize * bl / 2 )
	r.z -= zSize * bl;
      while ( r.z < -zSize * bl / 2)
	r.z += zSize * bl;
    }
    return r;
  };

  // periodic boundary conditions
  bool isPeriodicX;
  bool isPeriodicY;
  bool isPeriodicZ;
  
  // boundaries of box
  TYPE_FLOAT xMin;
  TYPE_FLOAT xMax;
  TYPE_FLOAT yMin;
  TYPE_FLOAT yMax;
  TYPE_FLOAT zMin;
  TYPE_FLOAT zMax;
  
  // length of box
  TYPE_FLOAT xSize;
  TYPE_FLOAT ySize;
  TYPE_FLOAT zSize;
 protected:
 private:

};

#endif // PERIODICBOUNDARY_H
