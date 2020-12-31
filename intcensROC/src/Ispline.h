#ifndef __ISPLINE_H__
#define __ISPLINE_H__

#include <iostream>

using namespace Eigen;
using namespace std;

static double Ispline_3(int m, double t, VectorXd &knot)
{
  if (t<=knot[m])
    return 0;
  if (t<=knot[m+1] && t>knot[m])
	return (t-knot[m])*(t-knot[m])*(t-knot[m])
    /((knot[m+1]-knot[m])*
      (knot[m+2]-knot[m])*
      (knot[m+3]-knot[m]));
    if (t<=knot[m+2] && t>knot[m+1])
    {
      double res1=(knot[m+1]-knot[m])*(knot[m+1]-knot[m])
      			  /((knot[m+2]-knot[m])*(knot[m+3]-knot[m]));
      
      double n2=-(t*t*t-knot[m+1]*knot[m+1]*knot[m+1])
                 +1.5*(knot[m]+knot[m+2])
                 *(t*t-knot[m+1]*knot[m+1])
                 -3*knot[m]*knot[m+2]*(t-knot[m+1]);
      double d2=	(knot[m+2]-knot[m+1])
                 *(knot[m+2]-knot[m])
                 *(knot[m+3]-knot[m]);
      double res2=n2/d2;
      double n3=	-(t*t*t-knot[m+1]*knot[m+1]*knot[m+1])
                +1.5*(knot[m+1]+knot[m+3])
                *(t*t-knot[m+1]*knot[m+1])
                -3*knot[m+1]*knot[m+3]*(t-knot[m+1]);
      double d3=	 (knot[m+2]-knot[m+1])
              *(knot[m+3]-knot[m+1])
              *(knot[m+3]-knot[m]);
      double res3=n3/d3;
      return res1+res2+res3;
    }
    if (t<knot[m+3] && t>knot[m+2])	
      return 1- (knot[m+3]-t)*(knot[m+3]-t)*(knot[m+3]-t)
      /((knot[m+3]-knot[m+2])
      *(knot[m+3]-knot[m+1])
      *(knot[m+3]-knot[m]));
      if (t>=knot[m+3])    
        return 1;

    return 0;
}

static MatrixXd Ispline(int pn, VectorXd &x, VectorXd &knot)
{
  MatrixXd res(pn, x.size());
  res.fill(1);
  for(int i = 1; i < pn; i++)
  {
    for(int j = 0; j < x.size(); j++)
      res(i,j) = Ispline_3(i-1, x(j), knot);
  }
  return res;
}

// function to return a I spline array.
static VectorXd IsplineX(int pn, double x, VectorXd &knot)
{
  VectorXd res(pn);
  res.fill(1);
  for(int i = 1; i < pn; i++)
  {
      res(i) = Ispline_3(i-1, x, knot);
  }
  return res;
}




#endif
