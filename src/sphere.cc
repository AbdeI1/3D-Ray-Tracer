#include "sphere.h"

std::vector<double> Sphere::rayIntersection(Ray r) {
  std::vector<double> intersections;
  // solves quadratic formula for all points on ray which are
  // length radius away from position
  // formula was given in lecture slides
  double A = pow(r.direction.magnitude(), 2);
  double B = 2*(r.direction*(r.origin - position));
  double C = pow((r.origin - position).magnitude(), 2) - pow(radius, 2);
  double disc = B*B - 4*A*C;
  if(disc == 0) {
    intersections.push_back(-B/(2*A));
  } else if (disc > 0) {
    intersections.push_back((-B + sqrt(disc))/(2*A));
    intersections.push_back((-B - sqrt(disc))/(2*A));
  }
  return intersections;
}

Vector3 Sphere::getTangent(Vector3 p) {
  Vector3 N = getNormal(p);
  if(N.magnitude() < EPS) return Vector3(0);
  Vector3 T = Vector3(-N.y, N.x, 0);
  return T.unit();
}

Vector3 Sphere::getUV(Vector3 p) {
  Vector3 n = getNormal(p);
  double phi = acos(n.z);
  double v = phi / M_PI;
  double theta = atan2(n.y, n.x);
  double u = theta < 0 ? (theta + 2*M_PI) / (2*M_PI) : theta / (2*M_PI);
  return Vector3(u, v, 0);
}
