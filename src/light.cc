#include "light.h"

double Light::fatt(Vector3 p) {
  // light attenuation
  double d = position.dist(p);
  double denom = c1 + c2*d + c3*d*d;
  if(denom < EPS) return 1;
  return 1 / denom;
}

PointLight::PointLight(Vector3 v) {
  position = v;
}

DirectionalLight::DirectionalLight(Vector3 v) {
  position = v*-1000000; // light position approximated as super far away
  direction = v;
}
