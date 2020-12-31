#include <assert.h>
#include <cmath>
#include <boost/bind.hpp>
#include <boost/function.hpp>
#define TOL 0.0003

extern "C"
{
  typedef boost::function<double(double x)> bindtype;
  using boost::function;
  double fminbr(double a, double b, bindtype& fn,const double& tol);
  #define SQRT_EPSILON   1.4901161193847656e-10
}
