#ifndef MATRIX3_H
#define MATRIX3_H

#include <iostream>
#include <stdexcept>
#include "Vector3D.h"

/**
  3x3 Matrix
 */
class Matrix3
{
 public:
  // entries
  TYPE_FLOAT xx, xy, xz, yx, yy, yz, zx, zy, zz;

  inline Matrix3():xx(0),
			       xy(0),
			       xz(0),
			       yx(0),
			       yy(0),
			       yz(0),
			       zx(0),
			       zy(0),
			       zz(0) {} ;
  // diagonal matrix
  inline Matrix3(TYPE_FLOAT c):xx(c),
			       xy(0),
			       xz(0),
			       yx(0),
			       yy(c),
			       yz(0),
			       zx(0),
			       zy(0),
			       zz(c) {} ;
  inline Matrix3(TYPE_FLOAT xx_,
		 TYPE_FLOAT xy_,
		 TYPE_FLOAT xz_,
		 TYPE_FLOAT yx_,
		 TYPE_FLOAT yy_,
		 TYPE_FLOAT yz_,
		 TYPE_FLOAT zx_,
		 TYPE_FLOAT zy_,
		 TYPE_FLOAT zz_):xx(xx_),
				 xy(xy_),
				 xz(xz_),
				 yx(yx_),
				 yy(yy_),
				 yz(yz_),
				 zx(zx_),
				 zy(zy_),
				 zz(zz_) {};
  inline Matrix3(const Vector3D & r1, const Vector3D & r2):xx(r1.x*r2.x),
					   xy(r1.x*r2.y),
					   xz(r1.x*r2.z),
					   yx(r1.y*r2.x),
					   yy(r1.y*r2.y),
					   yz(r1.y*r2.z),
					   zx(r1.z*r2.x),
					   zy(r1.z*r2.y),
					   zz(r1.z*r2.z) {};

   inline Matrix3(const Matrix3 & m):xx(m.xx),
				    xy(m.xy),
				    xz(m.xz),
				    yx(m.yx),
				    yy(m.yy),
				    yz(m.yz),
				    zx(m.zx),
				    zy(m.zy),
				    zz(m.zz) {};

  inline Matrix3& operator=(const Matrix3 & m)
    {
      xx = m.xx;
      xy = m.xy;
      xz = m.xz;
      yx = m.yx;
      yy = m.yy;
      yz = m.yz;
      zx = m.zx;
      zy = m.zy;
      zz = m.zz;
      return *this;
    };
  TYPE_FLOAT trace()
  {
    return xx+yy+zz;
  };

  /// \todo explain
  TYPE_FLOAT m ()
  {
    return xx*yy + yy*zz + zz*xx - xy*yx - yz*zy - xz*zx;
  };
  TYPE_FLOAT get(int i, int j)
  {
    if (i<0 || j <0 || i>2 || j>2)
      throw std::invalid_argument( "matrix index out of bounds" );
    if (i==0)
      {
	if (j== 0)
	  return xx;
	if (j== 1)
	  return xy;
	if (j== 2)
	  return xz;
      }
    if (i==1)
      {
	if (j== 0)
	  return yx;
	if (j== 1)
	  return yy;
	if (j== 2)
	  return yz;
      }
    if (i==2)
      {
	if (j== 0)
	  return zx;
	if (j== 1)
	  return zy;
	if (j== 2)
	  return zz;
      }
    return nan("");
  };
  void set(int i, int j, TYPE_FLOAT d)
  {
    if (i<0 || j <0 || i>2 || j>2)
      throw std::invalid_argument( "matrix index out of bounds" );
    if (i==0)
      {
	if (j== 0)
	  xx = d;
	if (j== 1)
	  xy = d;
	if (j== 2)
	  xz = d;
      }
    if (i==1)
      {
	if (j== 0)
	  yx = d;
	if (j== 1)
	  yy = d;
	if (j== 2)
	  yz = d;
      }
    if (i==2)
      {
	if (j== 0)
	  zx = d;
	if (j== 1)
	  zy = d;
	if (j== 2)
	  zz = d;
      }
  };

 protected:
 private:

};

inline Matrix3 operator+(const Matrix3 & m1, const Matrix3 & m2)
{
  return Matrix3(m1.xx + m2.xx,
		 m1.xy + m2.xy,
		 m1.xz + m2.xz,
		 m1.yx + m2.yx,
		 m1.yy + m2.yy,
		 m1.yz + m2.yz,
		 m1.zx + m2.zx,
		 m1.zy + m2.zy,
		 m1.zz + m2.zz);
};

inline Matrix3 operator-(const Matrix3 & m1, const Matrix3 & m2)
{
  return Matrix3(m1.xx - m2.xx,
		 m1.xy - m2.xy,
		 m1.xz - m2.xz,
		 m1.yx - m2.yx,
		 m1.yy - m2.yy,
		 m1.yz - m2.yz,
		 m1.zx - m2.zx,
		 m1.zy - m2.zy,
		 m1.zz - m2.zz);
};

inline Matrix3 operator*(const Matrix3 & m, const TYPE_FLOAT c)
{
  return Matrix3(m.xx * c,
		 m.xy * c,
		 m.xz * c,
		 m.yx * c,
		 m.yy * c,
		 m.yz * c,
		 m.zx * c,
		 m.zy * c,
		 m.zz * c);
};

inline Matrix3 operator*(const TYPE_FLOAT c, const Matrix3 & m)
{
  return Matrix3(m.xx * c,
		 m.xy * c,
		 m.xz * c,
		 m.yx * c,
		 m.yy * c,
		 m.yz * c,
		 m.zx * c,
		 m.zy * c,
		 m.zz * c);
};
inline Vector3D operator*(const Matrix3 & m, const Vector3D & r)
{
    return Vector3(m.xx * r.x + m.xy * r.y + m.xz * r.z,
                   m.yx * r.x + m.yy * r.y + m.yz * r.z,
                   m.zx * r.x + m.zy * r.y + m.zz * r.z);
};

std::ostream & operator<<(std::ostream & out, const Matrix3 & m);

#endif // MATRIX3_H
