#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

#include "sphere.h"
#include "cylinder.h"
#include "triangle.h"
#include "scene.h"

int main(int argc, char** argv) {

  // checks the number of input arguments
  if(argc < 2) {
    std::cout << "Error: not enough input arguments" << std::endl;
    std::cout << "usage:" << std::endl;
    std::cout << "  [executable_name] [input_file]" << std::endl;
    return 1;
  }

  // gets the filename from the input
  std::string filename = argv[1];
  std::ifstream fin(filename);
  if(!fin.is_open()) { // ensures that the file exists and can be opened
    std::cout << "Error: Could not open specified file" << std::endl;
    return 1;
  }

  std::vector<int> dim({-1, -1}); // the width and height of the output
  ppmFormat format = P3;

  Scene scene;
  scene.camera.hfov = 120; // default fov

  // default material with the color white
  scene.materials.push_back(new Material(Color(1)));

  scene.readFromFile(fin, dim, format);

  fin.close(); // input file not needed after more

  // check that dim was updated from its default value
  if(dim[0] == -1) {
    // if it was not updated, then the input did not contain imsize
    std::cout << "Error: file does not contain imsize parameter" << std::endl;;
    return 1;
  }

  // check vdir and up are not co-linear
  if(scene.camera.direction.cross(scene.camera.up).magnitude() < EPS) {
    std::cout << "Error: view direction and up direction are co-linear" << std::endl;;
    return 1;
  }

  // render the scene and get the image
  Image* img = scene.render(dim[0], dim[1]);

  // create output file with same name but different extension
  std::string outputfile = filename.substr(0, filename.rfind('.')) + ".ppm";
  img->save(outputfile, format); // save image
  delete img;
  return 0;
}
