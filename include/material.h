#ifndef MATERIAL_H_
#define MATERIAL_H_

#include "color.h"

// class for representing a material
class Material {
  public:
    static Material DEFAULT;
    Color ambient = Color(1);
    Color diffuse = Color(1);
    Color specular = Color(1);
    double kA = 0.2, kD = 0.8, kS = 0.3;
    double n = 2;
    double alpha = 0, eta = 0;
    Material() : Material(DEFAULT) {}
    Material(Color c) : Material(c, c, c) {}
    Material(Color c1, Color c2, Color c3) : ambient(c1), diffuse(c2), specular(c3) {}
    friend std::ostream& operator<<(std::ostream& strm, const Material& m);
};

#endif
