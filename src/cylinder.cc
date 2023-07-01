#include "cylinder.h"

Cylinder::Cylinder(Vector3 c, Vector3 d, double r, double l, Material* m) : Object(c, d, m), radius(r), length(l) {
  Vector3 z = direction.unit();
  Vector3 x = Vector3(0, 0, -1);
  if(x.cross(z).magnitude() < EPS) x = Vector3(0, 1, 0);
  Vector3 y = z.cross(x).unit();
  x = y.cross(z).unit();
  hDir = x;
}

Vector3 Cylinder::getNormal(Vector3 p) {
  Vector3 v1 = p - position;
  Vector3 v2 = direction * (direction.unit() * v1);
  Vector3 p2 = position + v2;
  return (p - p2).unit();
}

Vector3 Cylinder::getTangent(Vector3 p) {
  Vector3 N = getNormal(p);
  if(N.magnitude() < EPS) return Vector3(0);
  Vector3 T = Vector3(-N.y, N.x, 0);
  return T.unit();
}

std::vector<double> Cylinder::rayIntersection(Ray r) {
  std::vector<double> intersections;
  // solves quadratic formula for all points on ray which are 
  // length radius away from the line defined by position and direction
  // derivation is provided in the submitted writeup
  Vector3 dn = direction.unit();
  double A = pow(r.direction.magnitude(), 2) - pow(dn*r.direction, 2);
  double B = 2.*((r.direction * (r.origin - position)) - (dn * r.direction)*(dn * (r.origin - position)));
  double C = pow((r.origin - position).magnitude(), 2) - pow(dn * (r.origin - position), 2) - pow(radius, 2);
  double disc = B*B - 4*A*C;
  if(disc >= 0) {
    for(double t : {(-B + sqrt(disc))/(2*A), (-B - sqrt(disc))/(2*A)}) {
      Vector3 p = r.origin + r.direction*t;
      double d = sqrt(abs(pow(p.dist(position), 2) - pow(radius, 2)));
      // it's only an intersection with the cylinder if the projection of the point onto the line
      // is at most length/2 distance away from position
      if(d < length/2) intersections.push_back(t);
    }
  }
  return intersections;
}

Vector3 Cylinder::getUV(Vector3 p) {
  Vector3 n = getNormal(p);
  Vector3 d = hDir.cross(n).unit();
  double theta = acos(hDir * n);
  if(((d * -1) - direction).magnitude() < EPS)
    theta = 2*M_PI - theta;
  double u = theta / (2.*M_PI);
  double v = (((position - p) * direction.unit()) / length) + 0.5;
  return Vector3(u, v, 0);
}
