
#include "auxFun.h"

// anticlockwise or in a line (1) or otherwise (0) for (p1, p2, p3)
bool acw(double p1x1, double p1x2, double p2x1, double p2x2, double p3x1, double p3x2) {
  return (p1x1*p2x2 + p2x1*p3x2 + p3x1*p1x2 - p1x2*p2x1 - p2x2*p3x1 - p3x2*p1x1) >= 0;
}

// clockwise or in a line (1) or otherwise (0) for (p1, p2, p3)
bool cw(double p1x1, double p1x2, double p2x1, double p2x2, double p3x1, double p3x2) {
  return (p1x1*p2x2 + p2x1*p3x2 + p3x1*p1x2 - p1x2*p2x1 - p2x2*p3x1 - p3x2*p1x1) <= 0;
}

// penalized normal
double pNorm(double s, int l, int n) {
  return fabs(s)/sqrt(double(l)) - sqrt(2+2*log(double(n)/double(l)));
}
double ubPenNorm(double s, int l, int n, double q) {
  return s/double(l) + (q + sqrt(2+2*log(double(n)/double(l))))/sqrt(double(l));
}
double lbPenNorm(double s, int l, int n, double q) {
  return s/double(l) - (q + sqrt(2+2*log(double(n)/double(l))))/sqrt(double(l));
}

// Poisson case
double Pois(double s, int l, int n) {
  double mu = 1.;
  double ms = s/double(l);
  return sqrt(2*double(l)*(ms*log(ms/mu) + mu - ms));
}



