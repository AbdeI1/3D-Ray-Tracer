#ifndef LIGHT_H_
#define LIGHT_H_

#include "object.h"

class Light : public Object {
  public:
    Color color;
    double c1 = 1, c2 = 0, c3 = 0;
    virtual Vector3 getL(Vector3 p) = 0;
    virtual double fatt(Vector3 p);
};

class PointLight : public Light {
  public:
    PointLight(Vector3 v);
    Vector3 getL(Vector3 p) { 
      return p - position;
    }
};

class DirectionalLight : public Light {
  public:
    DirectionalLight(Vector3 v);
    Vector3 getL(Vector3 p) { return direction; }
    double fatt(Vector3 p) { return 1; }
};

#endif
