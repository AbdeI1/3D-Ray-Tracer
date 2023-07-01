#ifndef COLOR_H_
#define COLOR_H_

#include <cmath>
#include <iostream>

// simple class for colors
class Color {
  public:
    double r, g, b, a;
    Color() : Color(0.0) {}
    Color(double c) : Color(c, 1.0) {}
    Color(double c, double alpha) : Color(c, c, c, alpha) {}
    Color(double rr, double gg, double bb) : Color(rr, gg, bb, 1.0) {}
    Color(double rr, double gg, double bb, double aa);
    bool operator==(const Color& c);
    double& operator[](int i);

    // these two methods set/get parameters as integers between 0 and m
    int I(int i, int m = 255);
    Color& setI(int i, int c, int m = 255);
    
    friend std::ostream& operator<<(std::ostream&, const Color&);
};

#endif
