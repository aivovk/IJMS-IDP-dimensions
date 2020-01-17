#include "Matrix3.h"

std::ostream & operator<<(std::ostream & out, const Matrix3 & m)
{
  out<<"[ "<<m.xx<<" , "<<m.xy<<" , "<<m.xz<<"\n";
  out<<"  "<<m.yx<<" , "<<m.yy<<" , "<<m.yz<<"\n";
  out<<"  "<<m.zx<<" , "<<m.zy<<" , "<<m.zz<<" ]";
  return out;
}
