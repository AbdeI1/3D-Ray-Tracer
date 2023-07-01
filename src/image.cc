#include "image.h"

Image::Image(int ww, int hh) : w(ww), h(hh) {
  img = new Color[w*h];
  for(int i = 0; i < w*h; i++) img[i] = Color(0.0);
}

Image::Image(std::string filename) {
  std::ifstream fin(filename);
  if(!fin.is_open()) { // ensures that the file exists and can be opened
    std::cout << "Error: Could not open specified file" << std::endl;
    throw std::runtime_error("could not open file: " + filename);
  }

  std::string format;
  double maxColor;
  fin >> format >> w >> h >> maxColor;

  img = new Color[w*h];

  if(format == "P3") {
    double x;
    for(int i = 0; i < h; i++)
      for(int j = 0; j < w; j++)
        for(int k = 0; k < 3; k++) {
          fin >> x;
          (*this)[i][j][k] = x/maxColor;
        }
  } else if(format == "P6") {
    int c = fin.get();
    for(int i = 0; i < h; i++)
      for(int j = 0; j < w; j++)
        for(int k = 0; k < 3; k++) {
          c = fin.get();
          if(maxColor >= 256) c = c*256 + fin.get();
          (*this)[i][j].setI(k, c, maxColor);
        }
  } else throw std::runtime_error("unexpected ppm format: " + format);
}

Image::Image(const Image& o) {
  (*this) = o;
}

Image::~Image() {
  if(img) delete[] img;
}

Image& Image::operator=(const Image& o) {
  if(this != &o) {
    if(img) delete[] img;
    w = o.w;
    h = o.h;
    img = new Color[w*h];
    for(int i = 0; i < h; i++) {
      for(int j = 0; j < w; j++) {
        (*this)[i][j] = o[i][j];
      }
    }
  }
  return *this;
}

Color* Image::operator[](int i) const {
  return (img + i*w);
}

void Image::save(std::string filename, ppmFormat format) {
  std::ofstream fout(filename);

  switch(format) {
    case P3: fout << "P3\n"; break;
    case P6: fout << "P6\n"; break;
  }
  
  fout << w << " " << h << "\n";
  fout << 255 << std::endl;

  switch(format) {
    case P3:
      for(int i = 0; i < h; i++)
        for(int j = 0; j < w; j++)
          fout << (*this)[i][j].I(0) << " " << (*this)[i][j].I(1) << " " << (*this)[i][j].I(2) << '\n';
      break;
    case P6:
      for(int i = 0; i < h; i++)
        for(int j = 0; j < w; j++)
          fout << (unsigned char)(*this)[i][j].I(0) << (unsigned char)(*this)[i][j].I(1) << (unsigned char)(*this)[i][j].I(2);
      break;
  }
}

Color Image::bilinearInterpolate(double u, double v) {
  if(u < 0) u = 0; if(u > 1) u = 1;
  if(v < 0) v = 0; if(v > 1) v = 1;
  double x = u * (w - 1), y = v * (h - 1);
  int i = (int)y, j = (int)x;
  double beta = y - i, alpha = x - j;
  Color c(0);
  for(int k = 0; k < 3; k++) {
    c[k] = 
      (1-alpha)*(1-beta)*(*this)[i][j][k] + 
      (1-alpha)*(beta)*(*this)[i+1][j][k] + 
      (alpha)*(1-beta)*(*this)[i][j+1][k] + 
      (alpha)*(beta)*(*this)[i+1][j+1][k];
  }
  return c;
}
