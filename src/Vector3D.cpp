#include "Vector3D.h"

std::ostream & operator<<(std::ostream & out, const Vector3 & r)
{
    out<<r.x<<","<<r.y<<","<<r.z;
    return out;
}
