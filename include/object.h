#ifndef OBJECT_H_
#define OBJECT_H_

#include "vector3.h"
#include "material.h"
#include "image.h"
#include <vector>

// a class representing an object in the scene
class Object {
  public:
    Vector3 position;
    Vector3 direction;
    Material* material; // pointer to material of object
    Image* texture = nullptr; // pointer to texture (null if untextured)
    Image* bump = nullptr; // pointer to bump (null if no bumb)
    std::vector<Object*> children;
    Object(Vector3 p = Vector3(0),  Vector3 d = Vector3(0, 0, 1), Material* m = &Material::DEFAULT) : position(p), direction(d), material(m) {}
    virtual ~Object() {}
    virtual Vector3 getNormal(Vector3 p);
    virtual Vector3 getTangent(Vector3 p);
    virtual Vector3 getUV(Vector3 p) { return Vector3(0.5, 0.5, -1); }
    virtual std::vector<double> rayIntersection(Ray r) { return {}; } // computes all intersection with a given ray, returns array of times, t, when ray hits object
    virtual std::vector<std::pair<Object*, double>> rRayIntersection(Ray r);
    friend std::ostream& operator<<(std::ostream& strm, const Object& o);
};

#endif
