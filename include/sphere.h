#ifndef SPHERE_H_
#define SPHERE_H_

#include "object.h"

class Sphere : public Object {
  public:
    double radius;
    Sphere(Vector3 c = Vector3(0), double r = 1, Material* m = &Material::DEFAULT) : Object(c, Vector3(0, 0, 1), m), radius(r) {}
    Vector3 getTangent(Vector3 p);
    std::vector<double> rayIntersection(Ray r);
    Vector3 getUV(Vector3 p);
};

#endif
