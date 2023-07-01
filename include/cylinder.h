#ifndef CYLINDER_H_
#define CYLINDER_H_

#include "object.h"

class Cylinder : public Object {
  public:
    double radius;
    double length;
    Vector3 hDir;
    Cylinder(Vector3 c = Vector3(0), Vector3 d = Vector3(0, 1, 0), double r = 1, double l = 1, Material* m = &Material::DEFAULT);
    Vector3 getNormal(Vector3 p);
    Vector3 getTangent(Vector3 p);
    std::vector<double> rayIntersection(Ray r);
    Vector3 getUV(Vector3 p);
};

#endif
