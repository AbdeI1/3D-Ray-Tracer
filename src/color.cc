#include "color.h"

Color::Color(double rr, double gg, double bb, double aa) : r(rr), g(gg), b(bb), a(aa) {
  if(r < 0) r = 0; if(r > 1) r = 1;
  if(g < 0) g = 0; if(g > 1) g = 1;
  if(b < 0) b = 0; if(b > 1) b = 1;
  if(a < 0) a = 0; if(a > 1) a = 1;
}

bool Color::operator==(const Color& c) {
  return r == c.r && g == c.g && b == c.b && a == c.a;
}

double& Color::operator[](int i) {
  if(i < 0 || i >= 4) throw std::out_of_range("i not in range of Color");
  if(r < 0) r = 0; if(r > 1) r = 1;
  if(g < 0) g = 0; if(g > 1) g = 1;
  if(b < 0) b = 0; if(b > 1) b = 1;
  if(a < 0) a = 0; if(a > 1) a = 1;
  switch(i) {
    case 0: return r;
    case 1: return g;
    case 2: return b;
    case 3: return a;
  }
  return r;
}

int Color::I(int i, int m) {
  return round((*this)[i]*m);
}

Color& Color::setI(int i, int c, int m) {
  (*this)[i] = (double)c/m;
  return *this;
}

std::ostream& operator<<(std::ostream& strm, const Color& p) {
  strm << "{" << p.r << ", " << p.g << ", " << p.b << ", " << p.a << "}";
  return strm;
}
