#ifndef VECTOR3D_H
#define VECTOR3D_H

#include <cmath>
#include <iostream>
#include <stdexcept>

class Vector3;

///floating point precision used throughout
typedef double TYPE_FLOAT;

typedef Vector3 Vector3D;

/**
   3 component vector
 */
class Vector3
{
 public:
  // coordinates
  TYPE_FLOAT x, y, z;

  // default is (0, 0, 0)
  inline Vector3():x(0), y(0), z(0) {};
  inline Vector3(TYPE_FLOAT x_, TYPE_FLOAT y_, TYPE_FLOAT z_):x(x_), y(y_), z(z_) {};
  inline Vector3(const Vector3 & r):x(r.x), y(r.y), z(r.z) {};

  inline Vector3& operator=(const Vector3 & r)
    {
      x = r.x;
      y = r.y;
      z = r.z;
      return *this;
    };

  inline Vector3 operator-() const
  {
    return Vector3(-x, -y, -z);
  };

  inline Vector3& operator+=(const Vector3 & r)
  {
    x += r.x;
    y += r.y;
    z += r.z;
    return *this;
  };

  inline Vector3& operator-=(const Vector3 & r)
  {
    x -= r.x;
    y -= r.y;
    z -= r.z;
    return *this;
  };

  inline TYPE_FLOAT magnitude() const
  {
    return sqrt(x*x + y*y + z*z);
    //or sqrt(this*this);
  };

  inline TYPE_FLOAT magnitudeSquared() const
  {
    return x*x + y*y + z*z;
    // or this * this;
  };
  TYPE_FLOAT get(int i) const
  {
    if (i==0)
      return x;
    if (i==1)
      return y;
    if (i==2)
      return z;
    if (i<0 ||  i>2)
      throw std::invalid_argument( "vector index out of bounds" );
    return nan("");
  };
 protected:
 private:

};

std::ostream & operator<<(std::ostream & out, const Vector3 & r);

/// individually multiply x, y, and z of two vectors together
/// (not dot product)
inline Vector3 mult(const Vector3 & r1, const Vector3 & r2)
{
  return Vector3(r1.x * r2.x, r1.y * r2.y, r1.z * r2.z);
};

inline Vector3 operator+(const Vector3 & r1, const Vector3 & r2)
{
  return Vector3(r1.x + r2.x, r1.y + r2.y, r1.z + r2.z);
};

inline Vector3 operator-(const Vector3 & r1, const Vector3 & r2)
{
  return Vector3(r1.x - r2.x, r1.y - r2.y, r1.z - r2.z);
};


// dot product
inline TYPE_FLOAT operator*(const Vector3 & r1, const Vector3 & r2)
{
  return r1.x * r2.x + r1.y * r2.y + r1.z * r2.z;
};

// multiplication and division by a constant
inline Vector3 operator*(const Vector3 & r, const TYPE_FLOAT c)
{
    return Vector3(r.x * c, r.y * c, r.z * c);
};

inline Vector3 operator*(const TYPE_FLOAT c, const Vector3 & r)
{
    return Vector3(r.x * c, r.y * c, r.z * c);
};

inline Vector3 operator/(const Vector3 & r, const TYPE_FLOAT c)
{
    return Vector3(r.x / c, r.y / c, r.z / c);
};

inline bool operator==(const Vector3& a, const Vector3& b){
  return a.x == b.x && a.y == b.y && a.z == b.z;
}
#endif // VECTOR3D_H
