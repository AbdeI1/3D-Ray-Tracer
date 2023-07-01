#include "triangle.h"

Triangle::Triangle(const Vector3& a, const Vector3& b, const Vector3& c) : v1(a), v2(b), v3(c) {
  N = (v2 - v1).cross(v3 - v1);
  planeD = -(N * v1);
  cacheP = a;
  cacheB = Vector3(1, 0, 0);
}

Vector3 Triangle::getBarycentricCoordinates(Vector3 p) {
  if(p == cacheP) return cacheB;
  Vector3 e1 = v2 - v1;
  Vector3 e2 = v3 - v1;
  Vector3 ep = p - v1;
  double  d11 = e1 * e1;
  double  d12 = e1 * e2;
  double  d22 = e2 * e2;
  double det = d11 * d22 - d12 * d12;
  if(fabs(det) < EPS)
    return Vector3(-1);
  double beta = (d22 * (e1 * ep) - d12 * (e2 * ep)) / det;
  double gamma = (d11 * (e2 * ep) - d12 * (e1 * ep)) / det;
  double alpha = 1 - (beta + gamma);
  Vector3 coords(alpha, beta, gamma);
  cacheP = p;
  cacheB = coords;
  return coords;
}

Vector3 Triangle::getNormal(Vector3 p) {
  if(n1) {
    Vector3 baryCoords = getBarycentricCoordinates(p);
    Vector3 n = *n1 * baryCoords[0] + *n2 * baryCoords[1] + *n3 * baryCoords[2];
    return n.unit();
  }
  return N.unit();
}

void Triangle::computeT() {
  if(t1 == nullptr) return;
  double du1 = (*t2)[0] - (*t1)[0];
  double dv1 = (*t2)[1] - (*t1)[1];
  double du2 = (*t3)[0] - (*t1)[0];
  double dv2 = (*t3)[1] - (*t1)[1];
  double denom = du2 * dv1 - du1 * dv2;
  if(fabs(denom) < EPS) return;
  double det = 1 / denom;
  T = (v2 - v1) * (det * -dv2) - (v3 - v1) * (det * dv1);
}

Vector3 Triangle::getTangent(Vector3 p) {
  return T.unit();
}

std::vector<double> Triangle::rayIntersection(Ray r) {
  double denom = r.direction * N;
  if(fabs(denom) < EPS) return {};
  double t = -(N * r.origin + planeD) / denom;
  Vector3 p = r.origin + r.direction*t;
  Vector3 baryCoords = getBarycentricCoordinates(p);
  for(int i = 0; i < 3; i++)
    if(baryCoords[i] <= 0 || baryCoords[i] >= 1) 
      return {};
  return {t};
}

Vector3 Triangle::getUV(Vector3 p) {
  if(t1) {
    Vector3 baryCoords = getBarycentricCoordinates(p);
    Vector3 uv = *t1 * baryCoords[0] + *t2 * baryCoords[1] + *t3 * baryCoords[2];
    return uv;
  }
  return Vector3(0, 0, -1);
}

std::ostream& operator<<(std::ostream& strm, const Triangle& t) {
  strm << (Object)t << "\n";
  strm << "v1: " << t.v1 << ", v2: " << t.v2 << ", v3: " << t.v3;
  if(t.t1)
    strm << "\nt1: " << *t.t1 << ", t2: " << *t.t2 << ", t3: " << *t.t3;
  if(t.n1)
    strm << "\nn1: " << *t.n1 << ", n2: " << *t.n2 << ", n3: " << *t.n3;
  return strm;
}
