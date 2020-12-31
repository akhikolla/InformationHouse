
#ifndef auxFun_h
#define auxFun_h


#include <Rcpp.h>
//#include <ctime>
#include <iostream>
#include <vector>
//#include <assert.h>
#include <numeric>
#include <cmath>
#include <cfloat>

using namespace Rcpp;

// anticlockwise or in a line (1) or otherwise (0) for (p1, p2, p3)
bool acw(double, double, double, double, double, double);

// clockwise or in a line (1) or otherwise (0) for (p1, p2, p3)
bool cw(double, double, double, double, double, double);

// multiscale statistics
//      penalized normal
double pNorm(double, int, int);
//      Poisson case
double Pois(double, int, int);

// upper / lower bounds w.r.t. multiscale statistics
//      penalized normal
double ubPenNorm(double, int, int, double);
double lbPenNorm(double, int, int, double);

#endif /* auxFun_h */
