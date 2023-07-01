#ifndef CAMERA_H_
#define CAMERA_H_

#include "object.h"

#define NaN sqrt(-1)

enum cameraMode {
  perspective,
  parallel
};

// class for camera, inherits from object since camera needs position and direction
class Camera : public Object {
  public:
    Vector3 up; // up direction of camera
    // fov values, only one will be used, the other should be set to NaN if not used
    double hfov; // horizontal fov
    double vfov; // vertical fov
    cameraMode mode = perspective; // projection mod of camera
    Camera(Vector3 p = Vector3(0), Vector3 d = Vector3(0, 0, 1), Vector3 u = Vector3(0, 1, 0), double hf = NaN, double vf = NaN) : Object(p, d), up(u), hfov(hf), vfov(vf) {}
};

#endif
