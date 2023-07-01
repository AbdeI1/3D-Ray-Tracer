#ifndef VECTOR3_H_
#define VECTOR3_H_

#define _USE_MATH_DEFINES

#include <cmath>
#include <iostream>

#define EPS 0.0000001

// a simple class used for vector math, most function are self explanatory
class Vector3 {
  public:
    double x, y, z;
    Vector3() : Vector3(0) {}
    Vector3(double a) : Vector3(a, a, a) {}
    Vector3(double a, double b, double c) : x(a), y(b), z(c) {}
    bool operator==(const Vector3& v) const;
    double& operator[](int i);
    Vector3 operator+(const Vector3& v) const;
    Vector3 operator-(const Vector3& v) const;
    Vector3 operator*(double s) const;
    Vector3 operator/(double s) const;
    double operator*(const Vector3& v) const; //dot product
    Vector3 cross(const Vector3& v) const;
    double magnitude() const;
    void normalize();
    Vector3 unit() const; // normal vector in same direction
    double dist(const Vector3& v) const;
    friend std::ostream& operator<<(std::ostream& strm, const Vector3& v);
};

struct Ray {
  Vector3 origin;
  Vector3 direction;
};

#endif
