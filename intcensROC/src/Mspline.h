#ifndef __MSPLINE_H__
#define __MSPLINE_H__

#include <iostream>

using namespace Eigen;
using namespace std;

static double Mspline_3(int m, double t, VectorXd &knot)
{
  if( t<knot[m] || t>knot[m+3] )	
    return 0;
  if (t<=knot[m+1] && t>knot[m])
    return 3*(t-knot[m])*(t-knot[m])
    /(knot[m+1]-knot[m])
    /(knot[m+2]-knot[m])
    /(knot[m+3]-knot[m]);
    if (t<=knot[m+2] &&t>knot[m+1])
    {
      double res1=3*(t-knot[m])*(knot[m+2]-t)
      /(knot[m+2]-knot[m+1])
      /(knot[m+2]-knot[m])
      /(knot[m+3]-knot[m]);
      double res2=3*(knot[m+3]-t)*(t-knot[m+1])
        /(knot[m+2]-knot[m+1])
        /(knot[m+3]-knot[m+1])
        /(knot[m+3]-knot[m]);
        return res1+res2;
    }
    if(t<=knot[m+3] &&t>knot[m+2])                        
      return 3*(knot[m+3]-t)*(knot[m+3]-t)
      /(knot[m+3]-knot[m+2])
      /(knot[m+3]-knot[m+1])
      /(knot[m+3]-knot[m]);
   return 0;
}

static MatrixXd Mspline(int qn, VectorXd &x, VectorXd &knot)
{
  MatrixXd res2(qn, x.size());
  res2.fill(0); // this initialization is followint Yuan's code, I have not understand this
  for(int i = 1; i < qn; i++)
    for(int j = 0; j < x.size(); j++){
      res2(i,j) = Mspline_3(i-1, x(j), knot);
    }
    return res2;
}

// function to return a M spline array.
static VectorXd MsplineX(int qn, double x, VectorXd &knot)
{
  VectorXd res2(qn);
  res2.fill(0);
  for(int i = 1; i < qn; i++)
      res2(i) = Mspline_3(i-1, x, knot);
    return res2;
}




#endif




