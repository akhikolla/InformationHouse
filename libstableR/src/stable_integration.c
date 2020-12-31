/* stable/stable_integration.c
 * 
 * Functions to perform numerical integration used when calculating
 * the PDF and CDF of alpha-stable distributions. Based on GSL
 * numerical quadrature methods.
 *
 * Copyright (C) 2013. Javier Royuela del Val
 *                     Federico Simmross Wattenberg
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; version 3 of the License.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; If not, see <http://www.gnu.org/licenses/>.
 *
 *
 *  Javier Royuela del Val.
 *  E.T.S.I. Telecomunicación
 *  Universidad de Valladolid
 *  Paseo de Belén 15, 47002 Valladolid, Spain.
 *  jroyval@lpi.tel.uva.es    
 */
#include "stable.h"
#include "stable_integration.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "methods.h"

int stable_integration_METHODNAME(unsigned short method, char* name)
{
  switch (method)
    {
      case STABLE_QAG1:
        return sprintf(name,
        "QAG2: Adaptative 15 point Gauss-Kronrod rule");
      case STABLE_QAG2:
        return sprintf(name,
        "QAG2: Adaptative 21 point Gauss-Kronrod rule");
      case STABLE_QAG5:
        return sprintf(name,
        "QAG2: Adaptative 51 point Gauss-Kronrod rule");
      case STABLE_QUADSTEP:
        return sprintf(name,
        "QUADSTEP: Adaptative Bisection");
      case STABLE_QROMBPOL:
        return sprintf(name,
        "QROMBPOL: Romberg with Polinomial Extrapolation");
      case STABLE_QROMBRAT:
        return sprintf(name,
        "ROMBRAT: Romberg with Rational Extrapolation");
      case STABLE_QNG:
	return sprintf(name,
        "GSL_QNG: Non-adaptative Gauss-Kronrod rule 10, 21, 43 and 87 points");
    }

  sprintf(name,"Invalid method");
  return -1;
}

void
stable_integration_QAG1(StableDist *dist,double(function)(double, void*),
                   double a, double b,
                   double epsabs, double epsrel, unsigned short limit,
                   double *result, double *abserr)
{
  gsl_function F;

  F.function = function;
  F.params = (void*)dist;
  gsl_integration_qag(&F, a, b, epsabs, epsrel,
                      limit, 1, dist->gslworkspace, result, abserr);
}

void
stable_integration_QAG2(StableDist *dist,double(function)(double, void*),
                   double a, double b,
                   double epsabs, double epsrel, unsigned short limit,
                   double *result, double *abserr)
{
  gsl_function F;

  F.function = function;
  F.params = (void*)dist;
  gsl_integration_qag(&F, a, b, epsabs, epsrel,
                      limit, 2, dist->gslworkspace, result, abserr);
}
void
stable_integration_QAG5(StableDist *dist,double(function)(double, void*),
                   double a, double b,
                   double epsabs, double epsrel, unsigned short limit,
                   double *result, double *abserr)
{
  gsl_function F;

  F.function = function;
  F.params = (void*)dist;
  gsl_integration_qag(&F, a, b, epsabs, epsrel,
                      limit, 5, dist->gslworkspace, result, abserr);
}

void
stable_integration_QUADSTEP(StableDist *dist,double(function)(double, void*),
                   double a, double b,
                   double epsabs, double epsrel, unsigned short limit,
                   double *result, double *abserr)
{
  double fa,fc,fb;

  fa = function(a,(void*)dist);
  fc = function((a+b)*0.5,(void*)dist);
  fb = function(b,(void*)dist);
  *result = quadstep(function,(void*)dist,a,b,fa,fc,fb,
                       epsabs,epsrel,abserr,NULL,NULL);
}

void
stable_integration_QNG(StableDist *dist,double(function)(double, void*),
                   double a, double b,
                   double epsabs, double epsrel, unsigned short limit,
                   double *result, double *abserr)
{
  gsl_function F;
  size_t fcnt=0.0;

  //double c,d;
  //double res_aux=0.0,err_aux=0.0;
  //int warn=0;

  F.function = function;
  F.params = (void*)dist;
  gsl_integration_qng(&F,a,b,epsabs,epsrel,result,abserr,&fcnt);
/*
  if(*abserr<=epsabs || *abserr<=epsrel*fabs(*result)) return;
  else
    {
      *result=0;
      *abserr=0;
      warn=fcnt;
      while(warn<IT_MAX)
        {
          fcnt=0;
          d = (b-a)*0.25;
          c = a+d;
          d = b-d;
          gsl_integration_qng(&F,c,d,epsabs,epsrel,&res_aux,&err_aux,&fcnt);
          warn+=fcnt;
          *result=res_aux;
          *abserr=err_aux*err_aux;
          stable_integration(dist,function,a,c,fabs(*result*epsrel*0.5),epsrel,limit,&res_aux, &err_aux,STABLE_QAG2);
          *result+=res_aux;
          *abserr+=err_aux*err_aux;
          stable_integration(dist,function,d,b,fabs(*result*epsrel*0.5),epsrel,limit,&res_aux, &err_aux,STABLE_QAG2);
          *result+=res_aux;
          *abserr+=err_aux*err_aux;
          *abserr=sqrt(*abserr);
          if(*abserr<epsabs || *abserr<epsrel*fabs(*result))
           {
            printf(" %d ! ",warn);
            break; 
           }
        }
    }
*/
}

void
stable_integration(StableDist *dist,double(function)(double, void*),
                   double a, double b,
                   double epsabs, double epsrel, unsigned short limit,
                   double *result, double *abserr, unsigned short method)
{

  switch (method)
    {
      case STABLE_QAG2:
        stable_integration_QAG2(dist,function,a,b,epsabs,epsrel,limit,result,abserr);
        break;
      case STABLE_QUADSTEP:
        stable_integration_QUADSTEP(dist,function,a,b,epsabs,epsrel,limit,result,abserr);
        break;
      case STABLE_QROMBPOL:
        *result = qromb(function,(void*)dist,a,b,epsabs,epsrel,4,10,1,NULL,NULL,abserr);
        break;
      case STABLE_QROMBRAT:
        *result = qromb(function,(void*)dist,a,b,epsabs,epsrel,4,10,2,NULL,NULL,abserr);
	    break;
      case STABLE_QNG:
        stable_integration_QNG(dist,function,a,b,epsabs,epsrel,limit,result,abserr);
        break;
      case STABLE_QAG1:
        stable_integration_QAG1(dist,function,a,b,epsabs,epsrel,limit,result,abserr);
        break;
      case STABLE_QAG5:
        stable_integration_QAG5(dist,function,a,b,epsabs,epsrel,limit,result,abserr);
        break;
    }
}
