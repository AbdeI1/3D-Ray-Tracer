#include "vector3.h"

#define EPS 0.0000001

bool Vector3::operator==(const Vector3& v) const {
  return dist(v) < EPS;
}

double& Vector3::operator[](int i) {
  if(i < 0 || i >= 3) throw std::out_of_range("i not in range for vector");
  switch(i) {
    case 0: return x;
    case 1: return y;
    case 2: return z;
  }
  return x;
}

Vector3 Vector3::operator+(const Vector3& v) const {
  return Vector3(x+v.x, y+v.y, z+v.z);
}

Vector3 Vector3::operator-(const Vector3& v) const {
  return (*this)+(v*-1);
}

Vector3 Vector3::operator*(double s) const {
  return Vector3(x*s, y*s, z*s);
}

Vector3 Vector3::operator/(double s) const {
  return (*this)*(1/s);
}

double Vector3::operator*(const Vector3& v) const {
  return x*v.x + y*v.y + z*v.z;
}

Vector3 Vector3::cross(const Vector3& v) const {
  return Vector3(y*v.z-z*v.y, z*v.x-x*v.z, x*v.y-y*v.x);
}

double Vector3::magnitude() const {
  return sqrt(x*x + y*y + z*z);
}

void Vector3::normalize() {
  (*this) = unit();
}

Vector3 Vector3::unit() const {
  if(magnitude() < EPS) return (*this);
  return (*this)/magnitude();
}

double Vector3::dist(const Vector3& v) const {
  return ((*this)-v).magnitude();
}

std::ostream& operator<<(std::ostream& strm, const Vector3& v) {
  strm << "<" << v.x << ", " << v.y << ", " << v.z << ">";
  return strm;
}
