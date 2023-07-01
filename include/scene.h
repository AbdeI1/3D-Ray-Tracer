#ifndef SCENE_H_
#define SCENE_H_

#define MAX_RECURSION_DEPTH 10

#include <set>
#include "camera.h"
#include "image.h"
#include "light.h"

struct DepthCue {
  bool enabled = false;
  Color color;
  double aMin;
  double aMax;
  double distNear;
  double distFar;
};

struct Mesh {
  std::vector<Vector3> vertices;
  std::vector<Vector3> normals;
  std::vector<Vector3> textureCoords;
};


// class representng entire scene
class Scene {
  public:
    Camera camera; // the main camera to render from
    Color backgroundColor; // default color if no intersection is found
    double eta; // refraction index of surrounding environment
    DepthCue depthCue;
    Mesh mesh;
    std::vector<Object*> objects; // array of all the objects in the scene
    std::vector<Material*> materials; // array of all material in the scene (kept for memory management)
    std::vector<Image*> textures;
    std::vector<Image*> bumps;
    std::vector<Light*> lights;
    Scene();
    ~Scene();
    Image* render(int width, int height); // renders scene from perspective of camera and return pointer to Image created
    void readFromFile(std::ifstream&, std::vector<int>&, ppmFormat&);
  private:
    int rDepth = 0;
    std::vector<Object*> inside = {nullptr};
    std::pair<Object*, double> traceRay(Ray r, double eps = 0, std::set<Object*> ignore = std::set<Object*>()); // traces a ray, checks for intersections with each object
    double shadowCoefficient(Vector3 p, Light* l);
    Color shadeRay(Object* o, Ray r, double t); // shades ray with material
    Color traceAndShade(Ray r, double eps = 0, std::set<Object*> ignore = std::set<Object*>());
};

#endif
