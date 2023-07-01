#include "object.h"

Vector3 Object::getNormal(Vector3 p) {
  return (p - position).unit();
}

Vector3 Object::getTangent(Vector3 p) {
  return Vector3(0);
}

std::vector<std::pair<Object*, double>> Object::rRayIntersection(Ray r) {
  std::vector<std::pair<Object*, double>> intersections;
  auto v = rayIntersection(r);
  if(v.size() > 0) {
    for(auto t : v)
      intersections.push_back({this, t});
    for(auto o : children) {
      auto ci = o->rRayIntersection(r);
      intersections.insert(intersections.end(), ci.begin(), ci.end());
    }
  }
  return intersections;
}

std::ostream& operator<<(std::ostream& strm, const Object& o) {
  strm << "position: " << o.position << "\n";
  strm << "direction: " << o.direction << "\n";
  strm << "material: {\n";
  strm << *o.material << "}";
  return strm;
}
