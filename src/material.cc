#include "material.h"

Material Material::DEFAULT(Color(1.));

std::ostream& operator<<(std::ostream& strm, const Material& m) {
  strm << "Ambient: " << m.kA << " " << m.ambient << "\n";
  strm << "Diffuse: " << m.kD << " " << m.diffuse << "\n";
  strm << "Specular: " << m.kS << " " << m.specular << " ^" << m.n;
  return strm;
}
