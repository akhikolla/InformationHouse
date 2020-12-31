/* stable/stable_cdf.c
 *
 * Code for computing the CDF of an alpha-estable distribution.
 * Expresions presented in [1] are employed.
 *
 * [1] Nolan, J. P. Numerical Calculation of Stable Densities and
 *     Distribution Functions Stochastic Models, 1997, 13, 759-774
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

#include "methods.h"
#include <pthread.h>

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define SUBS_def 2

double stable_cdf_g1(double theta, void *args)
{
  StableDist *dist = (StableDist *)args;
  double g, V, aux;

  aux = (dist->beta_*theta+M_PI_2)/cos(theta);
  V = sin(theta)*aux/dist->beta_ + log(aux) + dist->k1;

  #ifdef DEBUG
  integ_eval++;
  #endif

  g = V + dist->xxipow;
  // Taylor: exp(-x) ~ 1-x en x ~ 0
  // If g < 1.52e-8 -> exp(-g) = (1-g) with double precision.
  // This way, we avoid evaluating an extra exp(x).
  if((g=exp(g)) < 1.522e-8)
    {
      return (1.0-g);
    }

  g = exp(-g);

//  if (g > 1.0) return 1.0;
//  else if(g < 0.0 || isnan(g)) return 0.0;
//  else return g;
  return g;
}

double stable_cdf_g2(double theta, void *args)
{
  StableDist *dist = (StableDist *)args;
  double g, cos_theta,aux,V;

  cos_theta = cos(theta);
  aux = (dist->theta0_+theta)*dist->alpha;
  V = log(cos_theta/sin(aux))*dist->alphainvalpha1 +
       + log(cos(aux-theta)/cos_theta) + dist->k1;

  #ifdef DEBUG
  integ_eval++;
  #endif

  g = V + dist->xxipow;
  // Taylor: exp(-x) ~ 1-x en x ~ 0
  // If g < 1.52e-8 -> exp(-g) = (1-g) with double precision.
  // This way, we avoid evaluating an extra exp(x).
  if((g=exp(g)) < 1.522e-8)
    {
      return (1.0-g);
    }

  g = exp(-g);

//  if (g > 1.0) return 1.0;
//  else if(g < 0.0 || isnan(g)) return 0.0;
//  else return g;
  return g;
}

double stable_cdf_g(double theta, void *args)
{
  StableDist * dist = (StableDist *)args;
  if (dist->ZONE == ALPHA_1) {
    return stable_cdf_g1(theta, args); }
  else if (dist->ZONE == CAUCHY) {
    return -1.0; }
  else {
    return stable_cdf_g2(theta, args); }
}

double stable_cdf_g_aux1(double theta, void *args)
{
  StableDist *dist = (StableDist *)args;
  double g, V,aux;

  aux = (dist->beta_*theta+M_PI_2)/cos(theta);
  V = sin(theta)*aux/dist->beta_ + log(aux) + dist->k1;

  g = V + dist->xxipow;

  #ifdef DEBUG
  integ_eval++;
  #endif

//  if (isnan(g)) { return -HUGE_VAL; }
//  else return g;
  return g;
}

double stable_cdf_g_aux2(double theta, void *args)
{
  StableDist *dist = (StableDist *)args;
  double g, cos_theta,aux,V;

  cos_theta = cos(theta);
  aux = (dist->theta0_+theta)*dist->alpha;
  V = log(cos_theta/sin(aux))*dist->alphainvalpha1 +
       + log(cos(aux-theta)/cos_theta) + dist->k1;

  g = V + dist->xxipow;

  #ifdef DEBUG
  integ_eval++;
  #endif
//  if (g > 1.0) return 1.0;
//  else if(g < 0.0 || isnan(g)) return 0.0;
//  else return g;
  return g;
}

double stable_cdf_g_aux(double theta, void *args)
{
  StableDist *dist = (StableDist *)args;

  if (dist->ZONE==ALPHA_1)
    return stable_cdf_g_aux1(theta,args);
  else
    return stable_cdf_g_aux2(theta,args);
}

void * thread_init_cdf(void *ptr_args)
{
  StableArgsCdf *args = (StableArgsCdf *)ptr_args;
  int counter_ = 0;

  while (counter_ < args->Nx)
    {
      args->cdf[counter_]=(*(args->ptr_funcion))(args->dist,args->x[counter_],
                                                 &(args->err[counter_]));
      counter_++;
    }
  pthread_exit(NULL);

  return NULL;
}

void stable_cdf(StableDist *dist, const double* x, const unsigned int Nx, double *cdf, double *err)
{
  int Nx_thread[THREADS],
      initpoint[THREADS],
      k,flag=0;
  void *status;
  pthread_t threads[THREADS];
  StableArgsCdf args[THREADS];

  if (err==NULL) {flag=1;err=(double*)malloc(Nx*sizeof(double));}

  // Distribute evaluation points between threads
  Nx_thread[0] = Nx/THREADS;
  if (0 < Nx%THREADS) Nx_thread[0]++;

  initpoint[0] = 0;
  for(k=1;k<THREADS;k++)
    {
      Nx_thread[k] = Nx/THREADS;
      if (k < Nx%THREADS) Nx_thread[k]++;
      initpoint[k] = initpoint[k-1] + Nx_thread[k-1];
    }

  // Create threads with a copy of the distribution

  for(k=0; k<THREADS; k++)
    {
      args[k].ptr_funcion = dist->stable_cdf_point;
      args[k].dist = stable_copy(dist);
      args[k].cdf  = cdf+initpoint[k];
      args[k].x    = x+initpoint[k];
      args[k].Nx   = Nx_thread[k];
      args[k].err  = err+initpoint[k];

      if(pthread_create(&threads[k], NULL, thread_init_cdf, (void *)&args[k]))
        {
          perror("Error while creating thread");
          if (flag==1) free(err);
          return;
        }
    }

/* Wait until every thread finishes */
  for(k=0; k<THREADS; k++)
    {
      pthread_join(threads[k], &status);
    }

/* Free dist copies */
  for(k=0; k<THREADS; k++)
    {
      stable_free(args[k].dist);
    }
  if (flag==1) free(err);
}

/******************************************************************************/
/*   CDF evaluatio                                       */
/******************************************************************************/

double
stable_integration_cdf(StableDist *dist, double(*integrand)(double,void*),
                  double(*auxiliar)(double,void*), double *err)
{
    int k,warnz[SUBS_def],method_[SUBS_def];
    double cdf=0,cdf1=0,err1=0;
    double theta[SUBS_def+1],g[SUBS_def+1];

    theta[0]=-dist->theta0_+THETA_TH/5;
    g[0]=stable_cdf_g(theta[0],(void*)dist);

    theta[SUBS_def]=M_PI_2-THETA_TH/5;
    g[SUBS_def]=stable_cdf_g(theta[SUBS_def],(void*)dist);

    method_[0] = METHOD2;
    method_[1] = METHOD2;

    if (dist->alpha>1.0 || (dist->alpha==1 && dist->beta_<0)) //Entonces max a la derecha
      {
        for (k=SUBS_def-1;k>=0;k--)
          {
            if (k>0) {
            theta[k]=zbrent(auxiliar,(void*)dist,theta[0],theta[k+1],
                            -log(g[k+1]*1e-2),1e-9*(theta[k+1]-theta[0]),&warnz[k]);
            }


            g[k]=stable_cdf_g(theta[k],(void*)dist);

            stable_integration(dist,integrand,
                               theta[k],theta[k+1],
                               max(cdf*relTOL,absTOL)/SUBS_def,relTOL,IT_MAX,
                               &cdf1,&err1,method_[SUBS_def-k-1]);
            cdf+=cdf1;
            *err+=err1*err1;
          }
          /*stable_integration(&F,&theta[0],SUBS_def+1,0.0,relTOL,
                               IT_MAX,dist->gslworkspace,&cdf,&err);*/
      }

    else if(dist->alpha<1.0 || (dist->alpha==1 && dist->beta_>0)) //Entonces max a la izqda
      {
        for (k=1;k<=SUBS_def;k++)
          {
            if(k<SUBS_def){
            theta[k]=zbrent(auxiliar,(void*)dist,theta[k-1],theta[SUBS_def],
                            -log(g[k-1]*1e-2),1e-9*(theta[SUBS_def]-theta[k-1]),&warnz[k]);
            }

            g[k]=stable_cdf_g(theta[k],(void*)dist);

            stable_integration(dist,integrand,
                               theta[k-1],theta[k],
                               max(cdf*relTOL,absTOL)/SUBS_def,relTOL,IT_MAX,
                               &cdf1,&err1,method_[k-1]);
            cdf+=cdf1;
            *err+=err1*err1;
          }
      }
  *err=sqrt(*err);
  //freopen("data_integrand.txt","w",file_integ);
  /*fprintf(FINTEG,"%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t\n",
          x,theta[0],theta[SUBS_def/2],theta[SUBS_def],g[0],g[SUBS_def/2],g[SUBS_def],pdf);*/

  return cdf;
}


/******************************************************************************/
/*   CDF in particular cases                                                */
/******************************************************************************/

double
stable_cdf_point_GAUSS(StableDist *dist, const double x, double *err)
{
  double x_=(x-dist->mu_0)/dist->sigma;
  *err = 0.0;

  return 0.5+0.5*gsl_sf_erf(x_*0.5);
}

double
stable_cdf_point_CAUCHY(StableDist *dist, const double x, double *err)
{
  double x_=(x-dist->mu_0)/dist->sigma;
  *err = 0.0;

  return 0.5+M_1_PI*atan(x_);
}

double
stable_cdf_point_LEVY(StableDist *dist, const double x, double *err)
{
  double xxi=(x-dist->mu_0)/dist->sigma-dist->xi;

  if (xxi>0 && dist->beta>0)
    return gsl_sf_erfc(sqrt(0.5/xxi));
  else if (xxi<0 && dist->beta<0)
    return gsl_sf_erfc(sqrt(-0.5/xxi));
  else return 0.0;
}

/******************************************************************************/
/*   CDF general case                                                       */
/******************************************************************************/
double
stable_cdf_point_ALPHA_1(StableDist *dist, const double x, double *err)
{
  double cdf=0;
  double x_;

  double(*integrand)(double,void *) = &stable_cdf_g1;
  double(*auxiliar)(double,void *)   = &stable_cdf_g_aux1;

  x_=(x-dist->mu_0)/dist->sigma;

  *err=0.0;

  if (dist->beta < 0.0)
   {
    x_ = -x_;
    dist->beta_ = -dist->beta;
   }
  else dist->beta_ = dist->beta;

  //dist->xxipow = exp(-PI*x_*0.5/dist->beta_);
  dist->xxipow = (-M_PI*x_*0.5/dist->beta_);

  cdf = stable_integration_cdf(dist,integrand,auxiliar,err);

  if (dist->beta>0)
    cdf = dist->c3*cdf;
  else
    cdf = 1.0 - dist->c3*cdf;

  return cdf;
}

double
stable_cdf_point_STABLE(StableDist *dist, const double x, double *err)
{
  double cdf=0;
  double x_, xxi;

  double(*integrand)(double,void *) = &stable_cdf_g2;
  double(*auxiliar)(double,void *)= &stable_cdf_g_aux2;
  x_=(x-dist->mu_0)/dist->sigma;
  xxi=x_-dist->xi;
  *err=0.0;

        /* If too close to xi*/
    if (fabs(xxi) < XXI_TH)
      {
        cdf = M_1_PI*(M_PI_2-dist->theta0);
        return cdf;
      }
    else if (xxi<0) /*F(x<xi,alpha,beta) = 1-F(-x,alpha,-beta)*/
      {
        dist->theta0_=-dist->theta0;/*theta0(alpha,-beta)=-theta0(alpha,beta)*/
        dist->beta_=-dist->beta;
	if (fabs(dist->theta0_+M_PI_2)<THETA_TH) return 0.0;
      }
    else
      {
        dist->theta0_=dist->theta0;
        dist->beta_=dist->beta;
        if (fabs(dist->theta0_+M_PI_2)<THETA_TH) return 1.0;
      }

    dist->xxipow=dist->alphainvalpha1*log(fabs(xxi));

    // Only get here if alpha==1 or stable evaluation zone.
    cdf = stable_integration_cdf(dist,integrand,auxiliar,err);

  if (xxi>0)
    cdf = dist->c1+dist->c3*cdf;
  else if (dist->alpha>1.0)
    cdf = - dist->c3*cdf;
  else
    cdf = 0.5 - (dist->theta0 + cdf)*M_1_PI;

  return cdf;
}

/******************************************************************************/
/*   CDF general point                                                 */
/******************************************************************************/

double
stable_cdf_point(StableDist *dist, const double x, double *err)
{
  double temp;

  if (err==NULL) err=&temp;

  return (dist->stable_cdf_point)(dist,x,err);
}
