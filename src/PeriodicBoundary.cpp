#include "PeriodicBoundary.h"

PeriodicBoundary::PeriodicBoundary()
{

}

PeriodicBoundary::PeriodicBoundary(bool _isPeriodicX,
        bool _isPeriodicY,
        bool _isPeriodicZ,
        TYPE_FLOAT _xMin,
        TYPE_FLOAT _xMax,
        TYPE_FLOAT _yMin,
        TYPE_FLOAT _yMax,
        TYPE_FLOAT _zMin,
        TYPE_FLOAT _zMax)
{
    isPeriodicX = _isPeriodicX;
    isPeriodicY = _isPeriodicY;
    isPeriodicZ = _isPeriodicZ;
    xMin = _xMin;
    xMax = _xMax;
    yMin = _yMin;
    yMax = _yMax;
    zMin = _zMin;
    zMax = _zMax;

    xSize = xMax - xMin;
    ySize = yMax - yMin;
    zSize = zMax - zMin;
    //std::cout<<xMin<<", "<<xMax<<", "<<yMin<<", "<<yMax<<", "<<zMin<<", "<<zMax<<", "<<xSize<<", "<<ySize<<", "<<zSize<<std::endl;
}

PeriodicBoundary::~PeriodicBoundary()
{
    //dtor
}


