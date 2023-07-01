#include "scene.h"

#include <sstream>

#include "sphere.h"
#include "cylinder.h"
#include "triangle.h"

# define PI 3.14159265358979323846

Scene::Scene() {
  mesh.vertices.push_back(Vector3(0));
  mesh.normals.push_back(Vector3(0));
  mesh.textureCoords.push_back(Vector3(0));
}

Scene::~Scene() {
  for(Light* l : lights) delete l;
  for(Material* m : materials) delete m;
  for(Image* t : textures) delete t;
  for(Image* b : bumps) delete b;
  for(Object* p : objects) delete p;
}

Image* Scene::render(int imgWidth, int imgHeight) {

  Image* img = new Image(imgWidth, imgHeight);
  double D = 10; // hardcoded distance to viewing plane

  // defines the w, u, v axes
  Vector3 w = camera.direction;
  Vector3 u = w.cross(camera.up).unit();
  Vector3 v = u.cross(w).unit();

  // if no camera fov is specified, sets it to 120
  if(std::isnan(camera.hfov) && std::isnan(camera.vfov)) camera.hfov = 120;

  double width, height;
  if(std::isnan(camera.vfov)) { // vertical fov is not specified
    width = 2*D*tan(camera.hfov*PI/360); // computes width with horizontal fov
    height = width * ((double)img->height()/img->width());
  } else {
    height = 2*D*tan(camera.vfov*PI/360); // computes height with vertical fov
    width = height * ((double)img->width()/img->height());
  }

  // defins the viewing plane boundries
  Vector3 ul, ur, ll, lr;
  ul = camera.position +  w*D - u*(width/2) + v*(height/2);
  ur = camera.position +  w*D + u*(width/2) + v*(height/2);
  ll = camera.position +  w*D - u*(width/2) - v*(height/2);
  lr = camera.position +  w*D + u*(width/2) - v*(height/2);

  Vector3 dh = (ur - ul)/(img->width()-1);
  Vector3 dv = (ll - ul)/(img->height()-1);

  // iterate through each pixel in image
  for(int i = 0; i < img->height(); i++) {
    for(int j = 0; j < img->width(); j++) {
      Vector3 p = ul + dv*i + dh*j;
      Ray r;
      switch(camera.mode) {
        case perspective:
          r.origin = camera.position; // computes perspective ray
          r.direction = p-r.origin;
          break;
        case parallel:
          r.origin = p - w*D; // computes parallel ray
          r.direction = w;
          break;
      }
      (*img)[i][j] = traceAndShade(r);
    }
  }

  return img;
}

std::pair<Object*, double> Scene::traceRay(Ray r, double eps, std::set<Object*> ignore) {
  double minT = -1;
  Object* minO = NULL;
  for(Object* p : objects) { // loops through each object in scene
    if(ignore.find(p) != ignore.end()) continue;
    for(double t : p->rayIntersection(r)) { // find all intersections with ray
      if(t*r.direction.magnitude() < eps) continue; // if intersection is before the ray origin, it is ignored
      if(minT == -1 || t < minT) { // keeps track of earliest intersection
        minT = t;
        minO = p;
      }
    }
  }
  return {minO, minT};
}

double Scene::shadowCoefficient(Vector3 p, Light* l) {
  double shadowAvg = 0;
  int total = 0;
  // jittering light around for soft shadoes
  for(int x = 0; x <= 0; x++) {
    for(int y = 0; y <= 0; y++) {
      for(int z = 0; z <= 0; z++) {
        double shadow = 1;
        Ray shadowRay;
        shadowRay.origin = p;
        Vector3 shadowLight = l->position;
        shadowLight = shadowLight + Vector3(x, y, z)*0.1;
        shadowRay.direction = (shadowLight - p);
        auto trace = traceRay(shadowRay, 0.01);
        std::set<Object*> Os = {};
        while(trace.second > 0 && trace.second < 1) {
          double a = trace.first->material->alpha;
          Os.insert(trace.first);
          shadow *= (1-a);
          shadowRay.origin = shadowRay.origin + shadowRay.direction*trace.second;
          shadowRay.direction = shadowLight - shadowRay.origin;
          trace = traceRay(shadowRay, 0.01, Os);
        }
        shadowAvg += shadow;
        total += 1;
      }
    }
  }
  shadowAvg /= total;
  return shadowAvg;
}

Color Scene::shadeRay(Object* o, Ray r, double t) {
  Vector3 p = r.origin + r.direction*t;
  Vector3 N = o->getNormal(p);
  Color c(0);
  double eta_i = inside[inside.size()-1] ? inside[inside.size()-1]->material->eta : eta;
  if(o == inside[inside.size()-1]) inside.pop_back();
  else inside.push_back(o);
  double eta_t = inside[inside.size()-1] ? inside[inside.size()-1]->material->eta : eta;
  Vector3 V = r.direction.unit()*-1;
  Color Oa = o->material->ambient;
  Color Od = o->material->diffuse;
  Color Os = o->material->specular;
  if(o->texture) {
    Vector3 uv = o->getUV(p);
    if(uv[2] != -1)
      Oa = Od = o->texture->bilinearInterpolate(uv[0], uv[1]);
  }
  if(o->bump) {
    Vector3 uv = o->getUV(p);
    if(uv[2] != -1) {
      Color c = o->bump->bilinearInterpolate(uv[0], uv[1]);
      Vector3 n;
      for(int i = 0; i < 3; i++)
        n[i] = c[i]*2 - 1;
      // TBN stuff
      Vector3 T = o->getTangent(p);
      Vector3 B = N.cross(T).unit();
      T = B.cross(N).unit();
      N = T*n.x + B*n.y + N*n.z;
    }
  }
  for(int i = 0; i < 3; i++) {
    // compute ambient light
    double aL = o->material->kA*Oa[i];
    if(aL < 0) aL = 0;
    c[i] += aL;
  }
  for(Light* l : lights) {
    double shadowAvg = shadowCoefficient(p, l);
    Vector3 L = l->getL(p).unit()*-1;
    Vector3 H = (L + V).unit();
    double att = l->fatt(p);
    for(int i = 0; i < 3; i++) {
      // compute diffuse light
      double dL = shadowAvg*att*l->color[i]*o->material->kD*Od[i]*(N * L);
      if(dL < 0) dL = 0;
      c[i] += dL;
      double x = N * H;
      if(x < 0) x = 0;
      // compute specular light
      double sL = shadowAvg*att*l->color[i]*o->material->kS*Os[i]*pow(x, o->material->n);
      if(sL < 0) sL = 0;
      c[i] += sL;
    }
  }
  // ray tracing
  if(o != inside[inside.size()-1]) N = N * -1;
  double F0 = pow((eta_t - eta_i)/(eta_t + eta_i), 2);
  Vector3 I = (r.direction*-1).unit();
  double a = I * N;
  Vector3 R = N*2*a - I;
  double Fr = F0 + (1 - F0)*pow(1 - a, 5);
  if (Fr < 0) Fr = 0;
  if (Fr > 1) Fr = 1;
  if(o->material->kS > 0 && rDepth < MAX_RECURSION_DEPTH) {
    Ray reflect;
    reflect.origin = p;
    reflect.direction = R;
    if(o == inside[inside.size()-1]) inside.pop_back();
    else inside.push_back(o);
    Color Rc = traceAndShade(reflect, 0.01);
    if(o == inside[inside.size()-1]) inside.pop_back();
    else inside.push_back(o);
    for(int i = 0; i < 3; i++)
      c[i] += Fr*Rc[i];
  }
  double Tr = (1 - Fr)*(1 - o->material->alpha);
  if (Tr < 0) Tr = 0;
  if (Tr > 1) Tr = 1;
  if(rDepth < MAX_RECURSION_DEPTH && Tr > 0 && eta_t > 0) {
    double refractionRatio = eta_i/eta_t;
    Vector3 T = (N*-1)*sqrt(1. - pow(refractionRatio, 2)*(1. - pow(a, 2))) + (N*a - I)*refractionRatio;
    if(!std::isnan(T.x)) {
      Ray transparent;
      transparent.origin = p;
      transparent.direction = T;
      Color Tc = traceAndShade(transparent, 0.01);
      for(int i = 0; i < 3; i++)
        c[i] += Tr*Tc[i];
    }
  }
  if(o == inside[inside.size()-1]) inside.pop_back();
  else inside.push_back(o);
  // depth cueing
  if(depthCue.enabled) {
    double d = p.dist(camera.position);
    double alpha;
    if(d <= depthCue.distNear) alpha = depthCue.aMax;
    else if (d >= depthCue.distFar) alpha = depthCue.aMin;
    else alpha = depthCue.aMin + (depthCue.aMax - depthCue.aMin)*((depthCue.distFar - d)/(depthCue.distFar - depthCue.distNear));
    for(int i = 0; i < 3; i++) c[i] = c[i]*alpha + depthCue.color[i]*(1 - alpha);
  }
  return c;
}

Color Scene::traceAndShade(Ray r, double eps, std::set<Object*> ignore) {
  rDepth++;
  auto trace = traceRay(r, eps, ignore);
  Color c = trace.first ? shadeRay(trace.first, r, trace.second) : backgroundColor;
  rDepth--;
  return c;
}

void Scene::readFromFile(std::ifstream& fin, std::vector<int>& dim, ppmFormat& format) {
  std::string line;
  int lineNumber = 1;
  while(std::getline(fin, line)) { // reads the file line by line
    std::stringstream ss(line);
    std::string w;
    ss >> w; // reads the line word by word
    // checks to see if first word matches an expected paramter
    if (w == "eye") {
      for(int i = 0; i < 3; i++) {
        // ensure that the next word exists
        if(!(ss >> w)) {
          std::cerr << "Error on line " << lineNumber <<  ": not enough arguments specified for the eye parameter" << std::endl;
          std::cerr << "expected specification:" << std::endl;
          std::cerr << "  eye [eyex] [eyey] [eyez]" << std::endl;
          throw new std::invalid_argument("input file not formatted correctly");
        }
        // validate the number and store it in pos
        try {
          camera.position[i] = stod(w); // attempts to convert the string to a number
        } catch(std::invalid_argument e) {
          // the string did not represent an double
          std::cerr << "Error on line " << lineNumber <<  ": invalid argument specified for eye (must be a double)" << std::endl;
          throw new std::invalid_argument("input file not formatted correctly");
        } catch(std::out_of_range e) {
          // the double was out of range for a double. Is this even possible? I don't know, but better safe than sorry
          std::cerr << "Error on line " << lineNumber <<  ": invalid argument specified for eye (must be a double)" << std::endl;
          throw new std::invalid_argument("input file not formatted correctly");
        }
      }
    } else if (w == "viewdir") {
      for(int i = 0; i < 3; i++) {
        // ensure that the next word exists
        if(!(ss >> w)) {
          std::cerr << "Error on line " << lineNumber <<  ": not enough arguments specified for the viewdir parameter" << std::endl;
          std::cerr << "expected specification:" << std::endl;
          std::cerr << "  viewdir [vdirx] [vdiry] [vdirz]" << std::endl;
          throw new std::invalid_argument("input file not formatted correctly");
        }
        // validate the number and store it in dir
        try {
          camera.direction[i] = stod(w); // attempts to convert the string to a number
        } catch(std::invalid_argument e) {
          // the string did not represent an double
          std::cerr << "Error on line " << lineNumber <<  ": invalid argument specified for viewdir (must be a double)" << std::endl;
          throw new std::invalid_argument("input file not formatted correctly");
        } catch(std::out_of_range e) {
          // the double was out of range for a double. Is this even possible? I don't know, but better safe than sorry
          std::cerr << "Error on line " << lineNumber <<  ": invalid argument specified for viewdir (must be a double)" << std::endl;
          throw new std::invalid_argument("input file not formatted correctly");
        }
      }
      if(camera.direction.magnitude() < EPS) {
        // vector has magnitude 0, does not specify a valid direction
        std::cerr << "Error on line " << lineNumber <<  ": vdir has magnitude 0, it does not specify a valid direction" << std::endl;
        throw new std::invalid_argument("input file not formatted correctly");
      }
      camera.direction.normalize();
    } else if (w == "updir") {
      for(int i = 0; i < 3; i++) {
        // ensure that the next word exists
        if(!(ss >> w)) {
          std::cerr << "Error on line " << lineNumber <<  ": not enough arguments specified for the updir parameter" << std::endl;
          std::cerr << "expected specification:" << std::endl;
          std::cerr << "  updir [upx] [upy] [upz]" << std::endl;
          throw new std::invalid_argument("input file not formatted correctly");
        }
        // validate the number and store it in dir
        try {
          camera.up[i] = stod(w); // attempts to convert the string to a number
        } catch(std::invalid_argument e) {
          // the string did not represent an double
          std::cerr << "Error on line " << lineNumber <<  ": invalid argument specified for updir (must be a double)" << std::endl;
          throw new std::invalid_argument("input file not formatted correctly");
        } catch(std::out_of_range e) {
          // the double was out of range for a double. Is this even possible? I don't know, but better safe than sorry
          std::cerr << "Error on line " << lineNumber <<  ": invalid argument specified for updir (must be a double)" << std::endl;
          throw new std::invalid_argument("input file not formatted correctly");
        }
        
      }
      if(camera.up.magnitude() < EPS) {
        // vector has magnitude 0, does not specify a valid direction
        std::cerr << "Error on line " << lineNumber <<  ": up has magnitude 0, it does not specify a valid direction" << std::endl;
        throw new std::invalid_argument("input file not formatted correctly");
      }
      camera.up.normalize();
    } else if (w == "hfov") {
      if(!(ss >> w)) {
        std::cerr << "Error on line " << lineNumber <<  ": not enough arguments specified for the hfov parameter" << std::endl;
        std::cerr << "expected specification:" << std::endl;
        std::cerr << "  hfov [fov]" << std::endl;
        throw new std::invalid_argument("input file not formatted correctly");
      }
      try {
        camera.hfov = stod(w); // attempts to convert the string to a number
      } catch(std::invalid_argument e) {
        // the string did not represent an double
        std::cerr << "Error on line " << lineNumber <<  ": invalid argument specified for hfov (must be a double)" << std::endl;
        throw new std::invalid_argument("input file not formatted correctly");
      } catch(std::out_of_range e) {
        // the double was out of range for a double. Is this even possible? I don't know, but better safe than sorry
        std::cerr << "Error on line " << lineNumber <<  ": invalid argument specified for hfov (must be a double)" << std::endl;
        throw new std::invalid_argument("input file not formatted correctly");
      }
      if(camera.hfov <= -180 || camera.hfov >= 180) {
        std::cerr << "Error on line " << lineNumber <<  ": fov must be between -180 and 180 exclusive" << std::endl;
        throw new std::invalid_argument("input file not formatted correctly");
      }
      camera.vfov = sqrt(-1);
    } else if (w == "vfov") {
      if(!(ss >> w)) {
        std::cerr << "Error on line " << lineNumber <<  ": not enough arguments specified for the hfov parameter" << std::endl;
        std::cerr << "expected specification:" << std::endl;
        std::cerr << "  vfov [fov]" << std::endl;
        throw new std::invalid_argument("input file not formatted correctly");
      }
      try {
        camera.vfov = stod(w); // attempts to convert the string to a number
      } catch(std::invalid_argument e) {
        // the string did not represent an double
        std::cerr << "Error on line " << lineNumber <<  ": invalid argument specified for vfov (must be a double)" << std::endl;
        throw new std::invalid_argument("input file not formatted correctly");
      } catch(std::out_of_range e) {
        // the double was out of range for a double. Is this even possible? I don't know, but better safe than sorry
        std::cerr << "Error on line " << lineNumber <<  ": invalid argument specified for vfov (must be a double)" << std::endl;
        throw new std::invalid_argument("input file not formatted correctly");
      }
      if(camera.vfov <= -180 || camera.vfov >= 180) {
        std::cerr << "Error on line " << lineNumber <<  ": fov must be between -180 and 180 exclusive" << std::endl;
        throw new std::invalid_argument("input file not formatted correctly");
      }
      camera.hfov = sqrt(-1);
    } else if(w == "imsize") {
      for(int i = 0; i < dim.size(); i++) {
        // ensure that the next word exists
        if(!(ss >> w)) {
          std::cerr << "Error on line " << lineNumber <<  ": not enough arguments specified for the imsize parameter" << std::endl;
          std::cerr << "expected specification:" << std::endl;
          std::cerr << "  imsize [width] [height]" << std::endl;
          throw new std::invalid_argument("input file not formatted correctly");
        }
        // validate the number and store it in dim
        try {
          dim[i] = stoi(w); // attempts to convert the string to a number
        } catch(std::invalid_argument e) {
          // the string did not represent an integer
          std::cerr << "Error on line " << lineNumber <<  ": invalid argument specified for imsize (must be a positive integer)" << std::endl;
          throw new std::invalid_argument("input file not formatted correctly");
        } catch(std::out_of_range e) {
          // the integer was out of range for an int
          std::cerr << "Error on line " << lineNumber <<  ": imsize argument is too big (must be a positive integer)" << std::endl;
          throw new std::invalid_argument("input file not formatted correctly");
        }
        if(dim[i] < 0) {
          // the integer was negative
          std::cerr << "Error on line " << lineNumber <<  ": invalid argument specified for imsize (must be a positive integer)" << std::endl;
          throw new std::invalid_argument("input file not formatted correctly");
        }
      }
    } else if (w == "bkgcolor") {
      double input[4];
      for(int i = 0; i < 4; i++) {
        // ensure that the next word exists
        if(!(ss >> w)) {
          std::cerr << "Error on line " << lineNumber << ": not enough arguments specified for the bkgcolor parameter" << std::endl;
          std::cerr << "expected specification:" << std::endl;
          std::cerr << "  bkgcolor [r] [g] [b] [eta]" << std::endl;
          throw new std::invalid_argument("input file not formatted correctly");
        }
        // validate the number and store it in bkgcolor
        double x;
        try {
          x = stod(w); // attempts to convert the string to a number
        } catch(std::invalid_argument e) {
          // the string did not represent an double
          std::cerr << "Error on line " << lineNumber <<  ": invalid argument specified for bkgcolor (must be a double)" << std::endl;
          throw new std::invalid_argument("input file not formatted correctly");
        } catch(std::out_of_range e) {
          // the integer was out of range for an int
          std::cerr << "Error on line " << lineNumber <<  ": invalid argument specified for bkgcolor (must be a double)" << std::endl;
          throw new std::invalid_argument("input file not formatted correctly");
        }
        input[i] = x;
      }
      backgroundColor = Color(input[0], input[1], input[2]);
      eta = input[3];
    } else if (w == "mtlcolor") {
      Material* m = new Material();
      double input[12];
      for(int i = 0; i < 12; i++) {
        // ensure that the next word exists
        if(!(ss >> w)) {
          std::cerr << "Error on line " << lineNumber << ": not enough arguments specified for the mtlcolor parameter" << std::endl;
          std::cerr << "expected specification:" << std::endl;
          std::cerr << "  mtlcolor [Odr] [Odg] [Odb] [Osr] [Osg] [Osb] [ka] [kd] [ks] [n] [alpha] [eta]" << std::endl;
          throw new std::invalid_argument("input file not formatted correctly");
        }
        // validate the number and store it in mtlcolor
        double x;
        try {
          x = stod(w); // attempts to convert the string to a number
        } catch(std::invalid_argument e) {
          // the string did not represent an double
          std::cerr << "Error on line " << lineNumber <<  ": invalid argument specified for mtlcolor (must be a double)" << std::endl;
          throw new std::invalid_argument("input file not formatted correctly");
        } catch(std::out_of_range e) {
          // the integer was out of range for an int
          std::cerr << "Error on line " << lineNumber <<  ": invalid argument specified for mtlcolor (must be a double)" << std::endl;
          throw new std::invalid_argument("input file not formatted correctly");
        }
        input[i] = x;
      }
      m->ambient = Color(input[0], input[1], input[2]);
      m->diffuse = Color(input[0], input[1], input[2]);
      m->specular = Color(input[3], input[4], input[5]);
      m->kA = input[6];
      m->kD = input[7];
      m->kS = input[8];
      m->n = input[9];
      m->alpha = input[10];
      m->eta = input[11];
      materials.push_back(m);
    } else if (w == "sphere") {
      double input[4];
      for(int i = 0; i < 4; i++) {
        // ensure that the next word exists
        if(!(ss >> w)) {
          std::cerr << "Error on line " << lineNumber << ": not enough arguments specified for the sphere parameter" << std::endl;
          std::cerr << "expected specification:" << std::endl;
          std::cerr << "  sphere [cx] [cy] [cz] [r]" << std::endl;
          throw new std::invalid_argument("input file not formatted correctly");
        }
        // validate the number and store it in input
        double x;
        try {
          x = stod(w); // attempts to convert the string to a number
        } catch(std::invalid_argument e) {
          // the string did not represent an double
          std::cerr << "Error on line " << lineNumber <<  ": invalid argument specified for sphere (must be a double)" << std::endl;
          throw new std::invalid_argument("input file not formatted correctly");
        } catch(std::out_of_range e) {
          // the string did not represent an double
          std::cerr << "Error on line " << lineNumber <<  ": invalid argument specified for sphere (must be a double)" << std::endl;
          throw new std::invalid_argument("input file not formatted correctly");
        }
        input[i] = x;
      }
      Sphere* s = new Sphere(Vector3(input[0], input[1], input[2]), input[3], materials[materials.size()-1]);
      if(textures.size() > 0) s->texture = textures[textures.size()-1];
      if(bumps.size() > 0) s->bump = bumps[bumps.size()-1];
      objects.push_back(s);
    } else if (w == "projection") {
      if(!(ss >> w)) {
        std::cerr << "Error on line " << lineNumber <<  ": not enough arguments specified for the projection parameter" << std::endl;
        std::cerr << "expected specification:" << std::endl;
        std::cerr << "  projection [perspective|parallel]" << std::endl;
        throw new std::invalid_argument("input file not formatted correctly");
      }
      // validate that the projection is either perspective or parallel
      if(w == "perspective") camera.mode = perspective;
      else if (w == "parallel") camera.mode = parallel;
      else {
        std::cerr << "Error on line " << lineNumber <<  ": invalid argument specified for projection (must be either \"perspective\" or \"parallel\")" << std::endl;
        throw new std::invalid_argument("input file not formatted correctly");
      }
    } else if (w == "cylinder") {
      double input[8];
      for(int i = 0; i < 8; i++) {
        // ensure that the next word exists
        if(!(ss >> w)) {
          std::cerr << "Error on line " << lineNumber << ": not enough arguments specified for the cylinder parameter" << std::endl;
          std::cerr << "expected specification:" << std::endl;
          std::cerr << "  cylinder [cx] [cy] [cz] [dx] [dy] [dz] [radius] [length]" << std::endl;
          throw new std::invalid_argument("input file not formatted correctly");
        }
        // validate the number and store it in input
        double x;
        try {
          x = stod(w); // attempts to convert the string to a number
        } catch(std::invalid_argument e) {
          // the string did not represent an double
          std::cerr << "Error on line " << lineNumber <<  ": invalid argument specified for cylinder (must be a double)" << std::endl;
          throw new std::invalid_argument("input file not formatted correctly");
        } catch(std::out_of_range e) {
          // the string did not represent an double
          std::cerr << "Error on line " << lineNumber <<  ": invalid argument specified for cylinder (must be a double)" << std::endl;
          throw new std::invalid_argument("input file not formatted correctly");
        }
        input[i] = x;
      }
      Cylinder* c = new Cylinder(Vector3(input[0], input[1], input[2]), Vector3(input[3], input[4], input[5]), input[6], input[7], materials[materials.size()-1]);
      if(textures.size() > 0) c->texture = textures[textures.size()-1];
      if(bumps.size() > 0) c->bump = bumps[bumps.size()-1];
      objects.push_back(c);
    } else if (w == "format") {
      if(!(ss >> w)) {
        std::cerr << "Error on line " << lineNumber <<  ": not enough arguments specified for the format parameter" << std::endl;
        std::cerr << "expected specification:" << std::endl;
        std::cerr << "  format [P3|P6]" << std::endl;
        throw new std::invalid_argument("input file not formatted correctly");
      }
      // validate that the format is either P3 or P6
      if(w == "P3") format = P3;
      else if (w == "P6") format = P6;
      else {
        std::cerr << "Error on line " << lineNumber <<  ": invalid argument specified for format (must be either \"P3\" or \"P6\")" << std::endl;
        throw new std::invalid_argument("input file not formatted correctly");
      }
    } else if (w == "light") {
      double input[6];
      for(int i = 0; i < 3; i++) {
        // ensure that the next word exists
        if(!(ss >> w)) {
          std::cerr << "Error on line " << lineNumber << ": not enough arguments specified for the light parameter" << std::endl;
          std::cerr << "expected specification:" << std::endl;
          std::cerr << "  light [x] [y] [z] [w] [r] [g] [b]" << std::endl;
          throw new std::invalid_argument("input file not formatted correctly");
        }
        // validate the number and store it in input
        double x;
        try {
          x = stod(w); // attempts to convert the string to a number
        } catch(std::invalid_argument e) {
          // the string did not represent an double
          std::cerr << "Error on line " << lineNumber <<  ": invalid argument specified for light (must be a double)" << std::endl;
          throw new std::invalid_argument("input file not formatted correctly");
        } catch(std::out_of_range e) {
          // the string did not represent an double
          std::cerr << "Error on line " << lineNumber <<  ": invalid argument specified for light (must be a double)" << std::endl;
          throw new std::invalid_argument("input file not formatted correctly");
        }
        input[i] = x;
      }
      if(!(ss >> w)) {
        std::cerr << "Error on line " << lineNumber << ": not enough arguments specified for the light parameter" << std::endl;
        std::cerr << "expected specification:" << std::endl;
        std::cerr << "  light [x] [y] [z] [w] [r] [g] [b]" << std::endl;
        throw new std::invalid_argument("input file not formatted correctly");
      }
      // validate the number and store it in input
      int t;
      try {
        t = stoi(w); // attempts to convert the string to a number
      } catch(std::invalid_argument e) {
        // the string did not represent an double
        std::cerr << "Error on line " << lineNumber <<  ": invalid argument specified for light (must be an int)" << std::endl;
        throw new std::invalid_argument("input file not formatted correctly");
      } catch(std::out_of_range e) {
        // the string did not represent an double
        std::cerr << "Error on line " << lineNumber <<  ": invalid argument specified for light (must be an int)" << std::endl;
        throw new std::invalid_argument("input file not formatted correctly");
      }
      if(t != 0 && t != 1) {
        std::cerr << "Error on line " << lineNumber <<  ": invalid argument specified for light (w must be either 0 or 1)" << std::endl;
        throw new std::invalid_argument("input file not formatted correctly");
      }
      for(int i = 3; i < 6; i++) {
        // ensure that the next word exists
        if(!(ss >> w)) {
          std::cerr << "Error on line " << lineNumber << ": not enough arguments specified for the light parameter" << std::endl;
          std::cerr << "expected specification:" << std::endl;
          std::cerr << "  light [x] [y] [z] [w] [r] [g] [b]" << std::endl;
          throw new std::invalid_argument("input file not formatted correctly");
        }
        // validate the number and store it in input
        double x;
        try {
          x = stod(w); // attempts to convert the string to a number
        } catch(std::invalid_argument e) {
          // the string did not represent an double
          std::cerr << "Error on line " << lineNumber <<  ": invalid argument specified for light (must be a double)" << std::endl;
          throw new std::invalid_argument("input file not formatted correctly");
        } catch(std::out_of_range e) {
          // the string did not represent an double
          std::cerr << "Error on line " << lineNumber <<  ": invalid argument specified for light (must be a double)" << std::endl;
          throw new std::invalid_argument("input file not formatted correctly");
        }
        input[i] = x;
      }
      Light* l;
      Vector3 v(input[0], input[1], input[2]);
      if(t == 1) l = new PointLight(v);
      else if (t == 0) l = new DirectionalLight(v);
      l->color = Color(input[3], input[4], input[5]);
      lights.push_back(l);
    } else if (w == "attlight") {
      double input[9];
      for(int i = 0; i < 3; i++) {
        // ensure that the next word exists
        if(!(ss >> w)) {
          std::cerr << "Error on line " << lineNumber << ": not enough arguments specified for the attlight parameter" << std::endl;
          std::cerr << "expected specification:" << std::endl;
          std::cerr << "  attlight [x] [y] [z] [w] [r] [g] [b] [c1] [c2] [c3]" << std::endl;
          throw new std::invalid_argument("input file not formatted correctly");
        }
        // validate the number and store it in input
        double x;
        try {
          x = stod(w); // attempts to convert the string to a number
        } catch(std::invalid_argument e) {
          // the string did not represent an double
          std::cerr << "Error on line " << lineNumber <<  ": invalid argument specified for attlight (must be a double)" << std::endl;
          throw new std::invalid_argument("input file not formatted correctly");
        } catch(std::out_of_range e) {
          // the string did not represent an double
          std::cerr << "Error on line " << lineNumber <<  ": invalid argument specified for attlight (must be a double)" << std::endl;
          throw new std::invalid_argument("input file not formatted correctly");
        }
        input[i] = x;
      }
      if(!(ss >> w)) {
        std::cerr << "Error on line " << lineNumber << ": not enough arguments specified for the attlight parameter" << std::endl;
        std::cerr << "expected specification:" << std::endl;
        std::cerr << "  attlight [x] [y] [z] [w] [r] [g] [b] [c1] [c2] [c3]" << std::endl;
        throw new std::invalid_argument("input file not formatted correctly");
      }
      // validate the number and store it in input
      int t;
      try {
        t = stoi(w); // attempts to convert the string to a number
      } catch(std::invalid_argument e) {
        // the string did not represent an double
        std::cerr << "Error on line " << lineNumber <<  ": invalid argument specified for attlight (must be an int)" << std::endl;
        throw new std::invalid_argument("input file not formatted correctly");
      } catch(std::out_of_range e) {
        // the string did not represent an double
        std::cerr << "Error on line " << lineNumber <<  ": invalid argument specified for attlight (must be an int)" << std::endl;
        throw new std::invalid_argument("input file not formatted correctly");
      }
      if(t != 0 && t != 1) {
        std::cerr << "Error on line " << lineNumber <<  ": invalid argument specified for attlight (w must be either 0 or 1)" << std::endl;
        throw new std::invalid_argument("input file not formatted correctly");
      }
      for(int i = 3; i < 9; i++) {
        // ensure that the next word exists
        if(!(ss >> w)) {
          std::cerr << "Error on line " << lineNumber << ": not enough arguments specified for the attlight parameter" << std::endl;
          std::cerr << "expected specification:" << std::endl;
          std::cerr << "  attlight [x] [y] [z] [w] [r] [g] [b] [c1] [c2] [c3]" << std::endl;
          throw new std::invalid_argument("input file not formatted correctly");
        }
        // validate the number and store it in input
        double x;
        try {
          x = stod(w); // attempts to convert the string to a number
        } catch(std::invalid_argument e) {
          // the string did not represent an double
          std::cerr << "Error on line " << lineNumber <<  ": invalid argument specified for light (must be a double)" << std::endl;
          throw new std::invalid_argument("input file not formatted correctly");
        } catch(std::out_of_range e) {
          // the string did not represent an double
          std::cerr << "Error on line " << lineNumber <<  ": invalid argument specified for light (must be a double)" << std::endl;
          throw new std::invalid_argument("input file not formatted correctly");
        }
        input[i] = x;
      }
      Light* l;
      Vector3 v(input[0], input[1], input[2]);
      if(t == 1) l = new PointLight(v);
      else if (t == 0) l = new DirectionalLight(v);
      l->color = Color(input[3], input[4], input[5]);
      l->c1 = input[6];
      l->c2 = input[7];
      l->c3 = input[8];
      lights.push_back(l);
    } else if (w == "depthcueing") {
      double input[7];
      for(int i = 0; i < 7; i++) {
        // ensure that the next word exists
        if(!(ss >> w)) {
          std::cerr << "Error on line " << lineNumber << ": not enough arguments specified for the depthcueing parameter" << std::endl;
          std::cerr << "expected specification:" << std::endl;
          std::cerr << "  depthcueing [dcr] [dcg] [dcb] [amax] [amin] [distmax] [distmin]" << std::endl;
          throw new std::invalid_argument("input file not formatted correctly");
        }
        // validate the number and store it in input
        double x;
        try {
          x = stod(w); // attempts to convert the string to a number
        } catch(std::invalid_argument e) {
          // the string did not represent an double
          std::cerr << "Error on line " << lineNumber <<  ": invalid argument specified for depthcueing (must be a double)" << std::endl;
          throw new std::invalid_argument("input file not formatted correctly");
        } catch(std::out_of_range e) {
          // the string did not represent an double
          std::cerr << "Error on line " << lineNumber <<  ": invalid argument specified for depthcueing (must be a double)" << std::endl;
          throw new std::invalid_argument("input file not formatted correctly");
        }
        input[i] = x;
      }
      depthCue.enabled = true;
      depthCue.color = Color(input[0], input[1], input[2]);
      depthCue.aMax = input[3];
      depthCue.aMin = input[4];
      depthCue.distFar = input[5];
      depthCue.distNear = input[6];
    } else if (w == "v") {
      double input[3];
      for(int i = 0; i < 3; i++) {
        // ensure that the next word exists
        if(!(ss >> w)) {
          std::cerr << "Error on line " << lineNumber << ": not enough arguments specified for the vertex parameter" << std::endl;
          std::cerr << "expected specification:" << std::endl;
          std::cerr << "  v [x] [y] [z]" << std::endl;
          throw new std::invalid_argument("input file not formatted correctly");
        }
        // validate the number and store it in input
        double x;
        try {
          x = stod(w); // attempts to convert the string to a number
        } catch(std::invalid_argument e) {
          // the string did not represent an double
          std::cerr << "Error on line " << lineNumber <<  ": invalid argument specified for vertex (must be a double)" << std::endl;
          throw new std::invalid_argument("input file not formatted correctly");
        } catch(std::out_of_range e) {
          // the string did not represent an double
          std::cerr << "Error on line " << lineNumber <<  ": invalid argument specified for vertex (must be a double)" << std::endl;
          throw new std::invalid_argument("input file not formatted correctly");
        }
        input[i] = x;
      }
      Vector3 v(input[0], input[1], input[2]);
      mesh.vertices.push_back(v);
    } else if (w == "vn") {
      double input[3];
      for(int i = 0; i < 3; i++) {
        // ensure that the next word exists
        if(!(ss >> w)) {
          std::cerr << "Error on line " << lineNumber << ": not enough arguments specified for the vertex normals parameter" << std::endl;
          std::cerr << "expected specification:" << std::endl;
          std::cerr << "  vn [x] [y] [z]" << std::endl;
          throw new std::invalid_argument("input file not formatted correctly");
        }
        // validate the number and store it in input
        double x;
        try {
          x = stod(w); // attempts to convert the string to a number
        } catch(std::invalid_argument e) {
          // the string did not represent an double
          std::cerr << "Error on line " << lineNumber <<  ": invalid argument specified for vertex normal (must be a double)" << std::endl;
          throw new std::invalid_argument("input file not formatted correctly");
        } catch(std::out_of_range e) {
          // the string did not represent an double
          std::cerr << "Error on line " << lineNumber <<  ": invalid argument specified for vertex normal (must be a double)" << std::endl;
          throw new std::invalid_argument("input file not formatted correctly");
        }
        input[i] = x;
      }
      Vector3 n(input[0], input[1], input[2]);
      mesh.normals.push_back(n.unit());
    } else if (w == "vt") {
      double input[2];
      for(int i = 0; i < 2; i++) {
        // ensure that the next word exists
        if(!(ss >> w)) {
          std::cerr << "Error on line " << lineNumber << ": not enough arguments specified for the texture Coordinates parameter" << std::endl;
          std::cerr << "expected specification:" << std::endl;
          std::cerr << "  vt [u] [v]" << std::endl;
          throw new std::invalid_argument("input file not formatted correctly");
        }
        // validate the number and store it in input
        double x;
        try {
          x = stod(w); // attempts to convert the string to a number
        } catch(std::invalid_argument e) {
          // the string did not represent an double
          std::cerr << "Error on line " << lineNumber <<  ": invalid argument specified for texture Coordinates  (must be a double)" << std::endl;
          throw new std::invalid_argument("input file not formatted correctly");
        } catch(std::out_of_range e) {
          // the string did not represent an double
          std::cerr << "Error on line " << lineNumber <<  ": invalid argument specified for texture Coordinates  (must be a double)" << std::endl;
          throw new std::invalid_argument("input file not formatted correctly");
        }
        input[i] -= (long)input[i];
        if(input[i] < 0) input[i] += 1;
        input[i] = x;
      }
      Vector3 t(input[0], input[1], 0);
      mesh.textureCoords.push_back(t);
    } else if (w == "f") {
      std::string line;
      for(int i = 0; i < 3; i++) {
        if(!(ss >> w)) {
          std::cerr << "Error on line " << lineNumber << ": not enough arguments specified for triangle face parameter" << std::endl;
          std::cerr << "expected specification:" << std::endl;
          std::cerr << "  f [v1{/vt1/vn1}] [v2{/vt2/vn2}] [v3{/vt3/vn3}]" << std::endl;
          throw new std::invalid_argument("input file not formatted correctly");
        }
        line += w + " ";
      }
      int i1, i2, i3, i4, i5, i6, i7, i8, i9;
      Triangle* t;
      try{
        if (sscanf(line.c_str(), "%d/%d/%d %d/%d/%d %d/%d/%d", &i1, &i2, &i3, &i4, &i5, &i6, &i7, &i8, &i9) == 9) {
          //success reading a face in v/t/n format; proceed accordingly
          t = new Triangle(mesh.vertices.at(i1), mesh.vertices.at(i4), mesh.vertices.at(i7));
          t->t1 = &(mesh.textureCoords.at(i2));
          t->t2 = &(mesh.textureCoords.at(i5));
          t->t3 = &(mesh.textureCoords.at(i8));
          t->n1 = &(mesh.normals.at(i3));
          t->n2 = &(mesh.normals.at(i6));
          t->n3 = &(mesh.normals.at(i9));
          t->computeT();
        } else if (sscanf(line.c_str(), "%d//%d %d//%d %d//%d", &i1, &i3, &i4, &i6, &i7, &i9) == 6) {
          //success reading a face in v//n format; proceed accordingly
          t = new Triangle(mesh.vertices.at(i1), mesh.vertices.at(i4), mesh.vertices.at(i7));
          t->n1 = &(mesh.normals.at(i3));
          t->n2 = &(mesh.normals.at(i6));
          t->n3 = &(mesh.normals.at(i9));
        } else if (sscanf(line.c_str(), "%d/%d %d/%d %d/%d", &i1, &i2, &i4, &i5, &i7, &i8) == 6) {
          //success reading a face in v/t format; proceed accordingly
          t = new Triangle(mesh.vertices.at(i1), mesh.vertices.at(i4), mesh.vertices.at(i7));
          t->t1 = &(mesh.textureCoords.at(i2));
          t->t2 = &(mesh.textureCoords.at(i5));
          t->t3 = &(mesh.textureCoords.at(i8));
          t->computeT();
        } else if (sscanf(line.c_str(), "%d %d %d", &i1, &i4, &i7) == 3) {
          //success reading a face in v format; proceed accordingly
          t = new Triangle(mesh.vertices.at(i1), mesh.vertices.at(i4), mesh.vertices.at(i7));
        } else {
          //error reading face data
          std::cerr << "Error on line " << lineNumber << ": could not read in triangle face parameter" << std::endl;
          std::cerr << "expected specification:" << std::endl;
          std::cerr << "  f [v1{/vt1/vn1}] [v2{/vt2/vn2}] [v3{/vt3/vn3}]" << std::endl;
          throw new std::invalid_argument("input file not formatted correctly");
        }
      } catch(std::out_of_range e) {
        std::cerr << "Error on line " << lineNumber << ": index out of range" << std::endl;
        std::cerr << e.what() << std::endl;
        throw new std::invalid_argument("input file not formatted correctly");
      }
      t->material = materials[materials.size() - 1];
      if(textures.size() > 0) t->texture = textures[textures.size()-1];
      if(bumps.size() > 0) t->bump = bumps[bumps.size()-1];
      objects.push_back(t);
    } else if (w == "texture") {
      if(!(ss >> w)) {
        std::cerr << "Error on line " << lineNumber << ": not enough arguments specified for texture parameter" << std::endl;
        std::cerr << "expected specification:" << std::endl;
        std::cerr << "  texture [ppmFile]" << std::endl;
        throw new std::invalid_argument("input file not formatted correctly");
      }
      textures.push_back(new Image(w));
    } else if (w == "bump") {
      if(!(ss >> w)) {
        std::cerr << "Error on line " << lineNumber << ": not enough arguments specified for bump parameter" << std::endl;
        std::cerr << "expected specification:" << std::endl;
        std::cerr << "  bump [ppmFile]" << std::endl;
        throw new std::invalid_argument("input file not formatted correctly");
      }
      bumps.push_back(new Image(w));
    }
    lineNumber++;
  }
}
