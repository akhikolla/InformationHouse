/* stable/stable_fit.c
 *
 * Functions employed by different methods of estimation implemented
 * in Libstable.
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
#include "mcculloch.h"

#include <gsl/gsl_complex.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_fft_real.h>


void stable_fft(double *data, const unsigned int length, double * y)
{
  //int i;

  memcpy ( (void *)y, (const void *) data, length*sizeof(double));

  gsl_fft_real_radix2_transform (y, 1, length);

  return;
}

double stable_loglikelihood(StableDist *dist, double *data, const unsigned int length)
{
  double *pdf=NULL;
  double l=0.0;
  int i;

  pdf=(double*)malloc(sizeof(double)*length);

  stable_pdf(dist,data,length,pdf,NULL);

  for(i=0;i<length;i++)
    {
      if (pdf[i]>0.0) l+=log(pdf[i]);
    }

  free(pdf);
  return l;
}

double stable_loglike_p(stable_like_params *params)
{
  double *pdf;
  double l=0.0;
  int i;

  pdf=(double*)malloc(sizeof(double)*(params->length));

  stable_pdf(params->dist,params->data,params->length,pdf,NULL);

  for(i=0;i<params->length;i++)
    {
      if (pdf[i]>0.0) {l+=log(pdf[i]);}
    }

  free(pdf);

  return l;
}

double stable_minusloglikelihood(const gsl_vector * theta, void * p)
{
/* Cost function to minimize, with the estimation of sigma and mu given by McCulloch at each iteration*/
  double alpha=1, beta=0, sigma=1.0, mu=0.0;
  double minusloglike=0;
  stable_like_params * params = (stable_like_params *) p;

  alpha = gsl_vector_get(theta,0);
  beta = gsl_vector_get(theta,1);

  /* Update sigma and mu with McCulloch. It needs nu_c nu_z*/
  czab(alpha, beta, params->nu_c, params->nu_z, &sigma, &mu);

  /* Check that the parameters are valid */
  if(stable_setparams(params->dist, alpha, beta, sigma, mu, 0) < 0)
    {
      return GSL_NAN;
    }
  else minusloglike = -stable_loglike_p(params);

  if (isinf(minusloglike) || isnan(minusloglike)) minusloglike=GSL_NAN;

  return minusloglike;
}

int compare (const void * a, const void * b)
{
/* qsort compare function */
  return ((*(double *)b < *(double *)a) - (*(double *)a < *(double *)b));
}

static inline void get_original(const gsl_vector *s,double *a,double *b,double *c,double *m)
{
  *a = M_2_PI*atan(gsl_vector_get(s,0))+1.0;
  *b = M_2_PI*atan(gsl_vector_get(s,1));
  *c = exp(gsl_vector_get(s,2));
  *m = gsl_vector_get(s,3);
}
static inline void set_expanded(gsl_vector *s,const double a,const double b,const double c,const double m)
{
  gsl_vector_set(s,0,tan(M_PI_2*(a-1.0)));
  gsl_vector_set(s,1,tan(M_PI_2*b));
  gsl_vector_set(s,2,log(c));
  gsl_vector_set(s,3,m);
}

double stable_minusloglikelihood_whole(const gsl_vector * theta, void * p)
{
/* Whole cost function to minimize in a 4D parameter space */
  double alpha=1, beta=0, sigma=1.0, mu=0.0;
  double minusloglike=0;
  stable_like_params * params = (stable_like_params *) p;

  get_original(theta,&alpha,&beta,&sigma,&mu);

  /* Check that the parameters are valid */
  if(stable_setparams(params->dist, alpha, beta, sigma, mu, 0) < 0)
    {
      perror("setparams error");
      return GSL_NAN;
    }
  else minusloglike = -stable_loglike_p(params);

  if (isinf(minusloglike) || isnan(minusloglike)) minusloglike=GSL_NAN;

  return minusloglike;
}

void stable_fit_init(StableDist *dist, const double * data, const unsigned int length, double *pnu_c,double *pnu_z)
{
  /* McCulloch estimation */

  double *sorted=NULL;
  double alpha0, beta0, sigma0, mu0;

  /* We need to sort the data to get percentiles */
  sorted = (double*)malloc(length*sizeof(double));
  memcpy ( (void *)sorted, (const void *) data, length*sizeof(double));
  qsort  ( sorted, length, sizeof(double), compare);

  /* Estimate the parameters. */
  stab((const double *) sorted,length,0,&alpha0,&beta0,&sigma0,&mu0);

  /* Set parameters in the distribution */
  if(stable_setparams(dist,alpha0,beta0,sigma0,mu0, 0)<0)
    {
      perror("INITIAL ESTIMATED PARAMETER ARE NOT VALID");
      return;
    }

  /* Get pnu_c and pnu_z needed for mle2d estimation */
  cztab(sorted, length, pnu_c, pnu_z);

  free(sorted);
  return;
}

int stable_fit_iter(StableDist *dist, const double * data, const unsigned int length,const double nu_c,const double nu_z)
{
  const gsl_multimin_fminimizer_type *T;
  gsl_multimin_fminimizer *s;

  gsl_multimin_function likelihood_func;

  gsl_vector *theta, *ss;

  unsigned int iter = 0;
  int status=0;
  double size=0;

  double a=1,b=0.0,c=1,m=0.0;
  stable_like_params par;

  par.dist=dist;
  par.data=(double *)data;
  par.length=length;
  par.nu_c=nu_c;
  par.nu_z=nu_z;

  /* Init: Dist must be initialized with McCulloch estimation */
  theta=gsl_vector_alloc(2);
  gsl_vector_set (theta, 0, dist->alpha);
  gsl_vector_set (theta, 1, dist->beta);

  #ifdef DEBUG
    Rprintf("%lf, %lf\n",gsl_vector_get (theta, 0),gsl_vector_get (theta, 1));
  #endif

  /* Initial steps */
  ss = gsl_vector_alloc (2);
  gsl_vector_set_all (ss, 0.01);

  /* Cost function */
  likelihood_func.n = 2; // Dimension 2 (alpha y beta)
  likelihood_func.f = &stable_minusloglikelihood;
  likelihood_func.params = (void *) (&par);  // Parametros de la funcion

  T = gsl_multimin_fminimizer_nmsimplex2rand;

  s = gsl_multimin_fminimizer_alloc (T, 2); /* Dimension 2*/

  gsl_multimin_fminimizer_set (s, &likelihood_func, theta, ss);

  #ifdef DEBUG
    Rprintf("5\n");
  #endif

  /* Iterar */
  do
    {
      iter++;
      status = gsl_multimin_fminimizer_iterate(s);
      size = gsl_multimin_fminimizer_size (s);
      status = gsl_multimin_test_size (size, 0.02);
    } while (status == GSL_CONTINUE && iter < 200);

//  if (status!=GSL_SUCCESS)
//    {
//      printf("Minimizer warning: %s\n",gsl_strerror(status));
//      fflush(stdout);
//    }

  /* Recover alpha and beta estimations */
  gsl_vector_free(theta);

  a = gsl_vector_get (s->x, 0);
  b = gsl_vector_get (s->x, 1);

  /* Estimate sigma and beta with McCulloch */
  czab(a, b, nu_c, nu_z, &c, &m);

  // Store estimation
  if (stable_setparams(dist,a,b,c,m,0)<0)
   {
    perror("FINAL ESTIMATED PARAMETER ARE NOT VALID\n");
   }

  gsl_vector_free(ss);
  gsl_multimin_fminimizer_free (s);

  return status;
}

int stable_fit(StableDist *dist, const double *data, const unsigned int length)
{
  double nu_c=0.0,nu_z=0.0;
  int status = 0;

  stable_fit_init(dist,data,length,&nu_c,&nu_z);
  status=stable_fit_iter(dist,data,length,nu_c,nu_z);

  return status;
}

int stable_fit_iter_whole(StableDist *dist, const double * data, const unsigned int length)
{
  const gsl_multimin_fminimizer_type *T;
  gsl_multimin_fminimizer *s;

  gsl_multimin_function likelihood_func;

  gsl_vector *theta, *ss;

  unsigned int iter = 0;
  int status=0;
  double size=0;

  double a=1,b=0.0,c=1,m=0.0;
  stable_like_params par;

  par.dist=dist;
  par.data=(double *)data;
  par.length=length;
  par.nu_c=0;
  par.nu_z=0;

  /* Inital params (with McCulloch) */
  theta=gsl_vector_alloc(4);
  set_expanded(theta,dist->alpha,dist->beta,dist->sigma,dist->mu_1);

  #ifdef DEBUG
    Rprintf("%lf, %lf, %lf, %lf\n",gsl_vector_get (theta, 0),gsl_vector_get (theta, 1),gsl_vector_get (theta, 2),gsl_vector_get (theta, 3));
  #endif

  /* Initial steps */
  ss = gsl_vector_alloc (4);
  gsl_vector_set_all (ss, 0.01);

  /* Cost function to minimize */
  likelihood_func.n = 4; // 4 Dimensions (alpha, beta, sigma, mu_0)
  likelihood_func.f = &stable_minusloglikelihood_whole;
  likelihood_func.params = (void *) (&par);  // Cost function arguments

  /* Minimizer creation */
  T = gsl_multimin_fminimizer_nmsimplex2rand;

  s = gsl_multimin_fminimizer_alloc (T, 4); /* 4 dimensions */

  /* Set cost function, initial guess and initial steps */
  gsl_multimin_fminimizer_set (s, &likelihood_func, theta, ss);


  #ifdef DEBUG
    Rprintf("5\n");
  #endif

  /* Start iterations */
  do
    {
      iter++;
      status = gsl_multimin_fminimizer_iterate(s);
      if (status!=GSL_SUCCESS) {
        perror("Minimizer warning:\n");
      }

      size   = gsl_multimin_fminimizer_size (s);
      status = gsl_multimin_test_size (size, 0.002);
      /*
      printf(" %03d\t size = %f a_ = %f  b_ = %f  c_ = %f  m_ = %f f_ = %f \n",iter,size,gsl_vector_get (s->x, 0),gsl_vector_get (s->x, 1),
                                                           gsl_vector_get (s->x, 2),gsl_vector_get (s->x, 3), gsl_multimin_fminimizer_minimum(s));
      */

    } while (status == GSL_CONTINUE && iter < 200);



  if (status!=GSL_SUCCESS)
    {
      perror("Minimizer warning");
    }

  /* Get last estimation */

  gsl_vector_free(theta);
  theta = gsl_multimin_fminimizer_x (s);
  get_original(theta,&a,&b,&c,&m);

  /* Set estimated parameters to the distribution and check if their have valid values*/
  if (stable_setparams(dist,a,b,c,m,0)<0)
   {
    perror("FINAL ESTIMATED PARAMETER ARE NOT VALID\n");
   }

  gsl_vector_free(ss);
  gsl_multimin_fminimizer_free (s);

  return status;
}

int stable_fit_whole(StableDist *dist, const double *data, const unsigned int length)
{
//  double nu_c=0.0,nu_z=0.0;
  int status=0;
//  stable_fit_init(dist,data,length,&nu_c,&nu_z);
//  printf("McCulloch %d sampless: %f %f %f %f\n",length,dist->alpha,dist->beta,dist->sigma,dist->mu_1);

  status = stable_fit_iter_whole(dist,data,length);

  return status;
}



double * load_rand_data(char * filename, int N)
{
  FILE * f_data;
  double * data;
  int i;

  if ((f_data = fopen(filename,"rt")) == NULL)
   {
    perror("Error when opening file with random data");
   }

  data=(double*)malloc(N*sizeof(double));

  for(i=0;i<N;i++)
   {
    if (EOF==fscanf(f_data,"%le\n",data+i))
     {
      perror("Error when reading data");
     }
   }
  return data;
}

int stable_fit_mle(StableDist *dist, const double *data, const unsigned int length) {
  return stable_fit_whole(dist,data,length);
}

int stable_fit_mle2d(StableDist *dist, const double *data, const unsigned int length) {
  return stable_fit(dist,data,length);
}

