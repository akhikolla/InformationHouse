/* stable/stable_pdf.c
 * 
 * Functions wrappers of GSL routines for random sample generation of
 * alpha-stable random variable.
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
/* According https://cran.r-project.org/doc/manuals/r-release/R-exts.html#Random-numbers,
 * no random number generator other than R's RNG must be used */
// #include <gsl/gsl_randist.h>
#include <math.h>
#include <R.h>
/*
void
stable_rnd_seed(StableDist * dist, unsigned long int s)
{
  gsl_rng_set(dist->gslrand, s);
}
*/

static inline double
stable_rnd_point(StableDist *dist)
{
// This implementation of the Chambers, Mallows and Stuck method
// is based on https://github.com/markveillette/stbl/blob/master/stblrnd.m
  // return dist->mu_1 +
  //        gsl_ran_levy_skew(dist->gslrand, dist->sigma, dist->alpha, dist->beta);

  double rnd;
  double V;
  double W;
  double alpha = dist->alpha;
  double beta  = dist->beta;
  double aux,B,S;

  if (alpha==2) {
    rnd = sqrt(2) * norm_rand();
  }
  else if (alpha==1.0 && beta==0.0) {
    rnd = tan(M_PI_2 * (2.0 * unif_rand() - 1.0));
  }
  else if (alpha==0.5 && fabs(beta)==1) {
    rnd = beta / pow(norm_rand(),2.0);
  }
  else if (beta==0.0) {
    V = M_PI_2 * (2.0 * unif_rand() - 1.0);
    W = -log(unif_rand());
    rnd = sin(alpha*V) / pow(cos(V),1.0/alpha) * pow(cos(V*(1-alpha))/W,(1-alpha)/alpha);
  }
  else if (alpha != 1.0) {
    V = M_PI_2 * (2.0 * unif_rand() - 1.0);
    W = -log(unif_rand());
    aux = beta * tan(M_PI_2*alpha);
    B   = atan(aux);
    S   = pow(1+aux*aux,0.5/alpha);
    rnd = S * sin(alpha*V+B) / pow(cos(V),1.0/alpha) * pow(cos((1-alpha)*V-B)/W,(1-alpha)/alpha);
  }
  else {
    V = M_PI_2 * (2.0 * unif_rand() - 1.0);
    W = -log(unif_rand());
    aux = M_PI_2 + beta*V;
    rnd = 1.0/M_PI_2 * (aux * tan(V) - beta * log(M_PI_2*W*cos(V)/aux));
  }

  if (alpha != 1) {
    rnd = dist->sigma * rnd + dist->mu_1;
  }
  else {
    aux = dist->sigma;
    rnd = aux * rnd + 2.0/M_PI * beta * aux * log(aux) + dist->mu_1;
  }

  return rnd;
}

void
stable_rnd(StableDist *dist, double *rnd, unsigned int n)
{
  //double *rnd;
  int i;

  //rnd = (double*)malloc(n*sizeof(double));
  if (rnd ==NULL) {
    perror("stable_rnd: NULL output pointer");
    return;
  }

  GetRNGstate();
  for(i=0;i<n;i++)
    {
      rnd[i]=stable_rnd_point(dist);
    }
  GetRNGstate();
  return;
}
