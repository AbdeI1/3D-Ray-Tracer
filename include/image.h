#ifndef IMAGE_H_
#define IMAGE_H_

#include "color.h"
#include <fstream>

enum ppmFormat {
  P3,
  P6
};

class Image {
  private:
    Color* img = nullptr; // array that holds all the color data
    int w, h; // width and height of image
  public:
    Image() : Image(100, 100) {}
    Image(int ww, int hh);
    Image(std::string filename);
    Image(const Image& o);
    ~Image();
    Image& operator=(const Image& o);
    Color* operator[](int i) const; // access the ith row of img
    int width() { return w; }
    int height() { return h; }
    void save(std::string filename, ppmFormat format = P3); // saves image to filename with the ppm format specified
    Color bilinearInterpolate(double u, double v);
};

#endif
