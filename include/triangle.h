#ifndef TRIANGLE_H_
#define TRIANGLE_H_

#include "object.h"

class Triangle : public Object {
  public:
    const Vector3& v1, v2, v3; 
    Vector3* n1 = nullptr, *n2 = nullptr, *n3 = nullptr;
    Vector3* t1 = nullptr, *t2 = nullptr, *t3 = nullptr;
    Vector3 T, N;
    double planeD;
    Triangle(const Vector3& a, const Vector3& b, const Vector3& c);
    void computeT();
    Vector3 cacheP;
    Vector3 cacheB;
    Vector3 getBarycentricCoordinates(Vector3 p);
    Vector3 getTangent(Vector3 p);
    Vector3 getNormal(Vector3 p);
    std::vector<double> rayIntersection(Ray r);
    Vector3 getUV(Vector3 p);
    friend std::ostream& operator<<(std::ostream& strm, const Triangle& t);
};

#endif
