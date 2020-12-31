/* stable/stable_dist.c
 *
 * Main Libstable source file. Definition of the StableDist structures
 * and auxiliary functions to manage alpha-stable distributions.
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
#include <gsl/gsl_errno.h>

#include "stable.h"

#include <pthread.h>
#include <math.h>
#include <unistd.h>

/*----------------------------------------------------------------------------*/
/*                             Public part                                    */
/*----------------------------------------------------------------------------*/

unsigned int stable_get_THREADS() { return THREADS; }

//TODO: Get number of cores in differents platforms following C99 standar
#ifdef __WIN32
void stable_set_THREADS(unsigned int value) {
// Had some problems to pass CRAN tests with this section:
/*
  SYSTEM_INFO sysinfo;
  GetSystemInfo(&sysinfo);
  THREADS = (unsigned int)sysinfo.dwNumberOfProcessors;
*/
  if (value == 0) THREADS = 12;
  else THREADS = value;
}
#else
void stable_set_THREADS(unsigned int value) {
/*
  if (value == 0) THREADS = sysconf(_SC_NPROCESSORS_ONLN);
  else THREADS = value;
  //printf("\nCPUs = %u\n",THREADS);
*/
  if (value == 0) THREADS = 12;
  else THREADS = value;
}
#endif


int stable_get_METHOD1() { return METHOD1; }
void stable_set_METHOD1(int value) { METHOD1 = value; }
int stable_get_METHOD2() { return METHOD2; }
void stable_set_METHOD2(int value) { METHOD2 = value; }
int stable_get_METHOD3() { return METHOD3; }
void stable_set_METHOD3(int value) { METHOD3 = value; }

unsigned int stable_get_IT_MAX() { return IT_MAX; }
void stable_set_IT_MAX(unsigned int value) { IT_MAX = value;}

unsigned int stable_get_INV_MAXITER() { return INV_MAXITER; }
void stable_set_INV_MAXITER(unsigned int value) { INV_MAXITER = value;}

double stable_get_relTOL() { return relTOL; }
void stable_set_relTOL(double value) {
/* Switch between integration methods according to relTOL */
  relTOL = value;
  if (value < 1e-12) {
    stable_set_METHOD1(STABLE_QNG);
    stable_set_METHOD2(STABLE_QAG2);
    stable_set_METHOD3(STABLE_QAG1);
  }
  else if (value < 1e-8) {
    stable_set_METHOD1(STABLE_QNG);
    stable_set_METHOD2(STABLE_QAG2);
    stable_set_METHOD3(STABLE_QAG1);
  }
  else {
    stable_set_METHOD1(STABLE_QNG);
    stable_set_METHOD2(STABLE_QAG1);
    stable_set_METHOD3(STABLE_QAG1);
  }
}

double stable_get_absTOL() { return absTOL; }
void stable_set_absTOL(double value)
{
  absTOL = value;
  //FACTOR = 1e-16;
}

/* Parameter thresholds */

/* When abs(alpha - 1)<ALPHA_TH alpha is set to 1 */
double stable_get_ALPHA_TH() { return ALPHA_TH; }
void stable_set_ALPHA_TH(double value) { ALPHA_TH = value; }

/* When 1-abs(beta)<BETA_TH beta is set to sign(beta)*1.0 */
/* When alpha = 1 and abs(beta)<BETA_TH beta is set to 0.0*/
double stable_get_BETA_TH() { return BETA_TH; }
void stable_set_BETA_TH(double value) { BETA_TH = value; }

/* When abs(x-xxi)<XXI_TH x is set to XXI */
double stable_get_XXI_TH() { return XXI_TH; }
void stable_set_XXI_TH(double value) { XXI_TH = value; }

/* When theta get closer than THETA_TH to integration interval limits theta is set to the limit value */
double stable_get_THETA_TH() { return THETA_TH; }
void stable_set_THETA_TH(double value) { THETA_TH = value; }

/* Debug purposes*/
FILE * stable_get_FINTEG() { return FINTEG; }
FILE * stable_set_FINTEG(char *name)
{
  FINTEG = fopen(name,"wt");
  return FINTEG;
}

/* Log-file configuration */
FILE * stable_get_FLOG() { return FLOG; }
FILE * stable_set_FLOG(char * name)
{
  FLOG = fopen(name,"wt");
  return FLOG;
}

int stable_setparams(StableDist *dist,
                     double alpha, double beta, double sigma, double mu,
                     int parametrization)
{
  int zona;

  if(dist==NULL)
    {
      perror("ERROR");
      return(-1);
    }
  if((zona = stable_checkparams(alpha,beta,sigma,mu,parametrization)) == NOVALID)
    {
    //  printf ("No valid parameters: %lf %lf %lf %lf %d\n",
    //         alpha, beta, sigma, mu, parametrization);
      return zona;
    }

  dist->alpha = alpha;
  dist->beta = beta;
  dist->sigma = sigma;

  switch (zona)
    {
      case STABLE_B1:
        dist->beta = (dist->beta > 0) ? 1.0 : -1.0;
      case STABLE:
        dist->alphainvalpha1 = alpha/(alpha-1.0);
        dist->xi = -beta*tan(0.5*alpha*M_PI);
        dist->theta0 = atan(-dist->xi)/alpha;
        //dist->k1 = pow(1.0+dist->xi*dist->xi,-0.5/(alpha-1.0));
        dist->k1 = -0.5/(alpha-1.0)*log(1.0+dist->xi*dist->xi);
        dist->S = pow(1.0+dist->xi*dist->xi,0.5/alpha);
        dist->Vbeta1 = dist->k1 - dist->alphainvalpha1 * log(dist->alpha)
                                + log(fabs(dist->alpha-1.0));
        dist->stable_pdf_point = &stable_pdf_point_STABLE;
        dist->stable_cdf_point = &stable_cdf_point_STABLE;

        if (alpha < 1.0)
          {
            dist->c1 = 0.5-dist->theta0*M_1_PI;
            dist->c2_part = alpha/((1.0-alpha)*M_PI);
            dist->c3 = M_1_PI;
          }
        else
          {
            dist->c1 = 1.0;
            dist->c2_part = alpha/((alpha-1.0)*M_PI);
            dist->c3 = -M_1_PI;
          }

        if (alpha>1)
         {
          AUX1=log(log(8.5358/(relTOL))/0.9599);/*3.76;*/
          AUX2=log(relTOL);/*-40;*/
         }
        else
         {
          AUX1=log(relTOL);/*-40;*/
          AUX2=log(log(8.5358/(relTOL))/0.9599);/*3.76;*/
         }

        break;

      case ALPHA_1_B1:
        dist->beta = (dist->beta > 0) ? 1.0 : -1.0;
      case ALPHA_1:
        dist->alpha=1;
        dist->c2_part = 0.5/fabs(beta);
        dist->alphainvalpha1 = 0.0;
        dist->xi = 0.0;
        dist->theta0 = M_PI_2;
        dist->k1 = log(2.0*M_1_PI);
        //dist->k1 = 2.0*M_1_PI;
        dist->S = 2.0*M_1_PI;
        dist->c1 = 0.0;
        dist->c3 = M_1_PI;
        dist->Vbeta1=2.0*M_1_PI/M_E;
        dist->stable_pdf_point = &stable_pdf_point_ALPHA_1;
        dist->stable_cdf_point = &stable_cdf_point_ALPHA_1;
        //XXI_TH = 10*EPS;
        if (beta<0)
         {
          AUX1=log(log(8.5358/(relTOL))/0.9599);/*4;*/
          AUX2=log(relTOL);/*-25;*/
         }
        else
         {
          AUX1=log(relTOL);/*-25;*/
          AUX2=log(log(8.5358/(relTOL))/0.9599);/*4;*/
         }
        break;

      case CAUCHY:
        dist->beta=0;
        dist->alpha=1;
        dist->c2_part = 0.0;
        dist->alphainvalpha1 = 0.0;
        dist->xi = 0.0;
        dist->theta0 = M_PI_2;
        dist->k1 = log(2.0*M_1_PI);
        //dist->k1 = 2.0*M_1_PI;
        dist->S = 2.0*M_1_PI;
        dist->c1 = 0.0;
        dist->c3 = M_1_PI;
        dist->Vbeta1=2.0*M_1_PI/M_E;
        dist->stable_pdf_point = &stable_pdf_point_CAUCHY;
        dist->stable_cdf_point = &stable_cdf_point_CAUCHY;
        break;

      case GAUSS:
        dist->alpha=2;
        dist->beta=0.0;
        dist->alphainvalpha1 = 2.0;
        dist->xi = 0.0;
        dist->theta0 = 0.0;
        dist->k1 = log(2.0);
        //dist->k1 = 2.0;
        dist->S = 2.0;
        dist->c1 = 1.0;
        dist->c2_part = 2.0*M_1_PI;
        dist->c3 = -M_1_PI;
        dist->Vbeta1=0.25;
        dist->stable_pdf_point = &stable_pdf_point_GAUSS;
        dist->stable_cdf_point = &stable_cdf_point_GAUSS;
        break;

      case LEVY:
        dist->alpha = 0.5;
        dist->beta = (2.0*(beta>0)-1.0);//  1 ó -1
        dist->alphainvalpha1 = -1.0;
        dist->xi = -dist->beta;
        dist->theta0 = 0.5*M_PI;
        dist->k1 = 0.0;
        dist->S = 1.0;
        dist->c1 = 0.0;
        dist->c2_part = 0.5*M_1_PI;
        dist->c3 = M_1_PI;
        dist->Vbeta1 = dist->k1 - dist->alphainvalpha1 * log(dist->alpha)
                                + log(fabs(dist->alpha-1.0));
        dist->stable_pdf_point = &stable_pdf_point_LEVY;
        dist->stable_cdf_point = &stable_cdf_point_LEVY;
        break;
    }

  if (parametrization == 0)
    {
      dist->mu_0 = mu;
      if (dist->alpha==1)
        dist->mu_1 = mu-dist->beta*2*M_1_PI*dist->sigma*log(dist->sigma);
      else dist->mu_1 = mu+dist->xi*dist->sigma;
    }
  else if (parametrization == 1)
    {
      dist->mu_1 = mu;
      if (dist->alpha==1)
        dist->mu_0 = mu+dist->beta*2*M_1_PI*dist->sigma*log(dist->sigma);
      else dist->mu_0 = mu-dist->xi*dist->sigma;
    }

  dist->theta0_ = dist->theta0;
  dist->beta_ = dist->beta;
  dist->xxipow = 0.0;
  dist->ZONE = zona;

  return zona;
}

int stable_checkparams(double alpha, double beta, double sigma, double mu,
                       int parametrization)
{
  /*Check parameters*/
  if (0.0 >= alpha || alpha > 2.0)
    {
      //printf("Alpha must lie between 0.0 and 2.0.");
      return NOVALID;
    }
  else if (beta < -1.0 || beta > 1.0)
    {
      //printf("Beta must lie between -1.0 and 1.0.");
      return NOVALID;
    }
  else if (sigma <= 0.0)
    {
      //printf("Sigma must be positive.");
      return NOVALID;
    }
  else if (isnan(mu) || isinf(mu))
    {
      //printf("Mu must be real.");
      return NOVALID;
    }
  else if (parametrization != 0 && parametrization != 1)
    {
      //printf("Only parametrizations 0 and 1 are accepted.");
      return NOVALID;
    }

  /*ZONE determination*/
  if ((2.0 - alpha) <= ALPHA_TH)
    return GAUSS;  //GAUSS
  else if (fabs(alpha-0.5) <= ALPHA_TH && fabs((fabs(beta)-1.0)) <= BETA_TH)
    return LEVY; //LEVY
  else if (fabs(alpha-1.0) <= ALPHA_TH && fabs(beta) <= BETA_TH)
    return CAUCHY;
  else if (fabs(alpha-1.0) <= ALPHA_TH) {
//    if (fabs(fabs(beta)-1) <= BETA_TH) return ALPHA_1_B1;
    return ALPHA_1;
  }
  else {
//    if (fabs(fabs(beta)-1) <= BETA_TH) return STABLE_B1;
    return STABLE;
  }
  /*When alpha=1,beta=1 Landau distribution is obtained,
    but it is calculated as in the general case*/
}

StableDist * stable_create(double alpha, double beta, double sigma, double mu,
                           int parametrization)
{
  /*gsl_error_handler_t * old_handler;
  old_handler = */gsl_set_error_handler (&error_handler);

  StableDist * dist = (StableDist *) malloc(sizeof (StableDist));

  if (dist == NULL)
    {
      perror("Error during distribution creation.");
      return NULL;
    }
  if((stable_setparams(dist,alpha,beta,sigma,mu,parametrization)) == NOVALID)
    {
      perror ("Error during distribution creation.");
      return NULL;
    }
  // gsl_rng_env_setup(); //leemos las variables de entorno
  dist->gslworkspace = gsl_integration_workspace_alloc(IT_MAX);
  // dist->gslrand = gsl_rng_alloc (gsl_rng_default);

  //Allow the distribution to use THREADS threads.
  stable_set_THREADS(THREADS);

  return dist;
}

StableDist * stable_copy(StableDist *src_dist)
{
  StableDist *dist;

  dist = stable_create(src_dist->alpha, src_dist->beta,
                       src_dist->sigma, src_dist->mu_0, 0);
  return dist;
}

void stable_free(StableDist *dist)
{
  if (dist == NULL)
    return;
  gsl_integration_workspace_free(dist->gslworkspace);
  // gsl_rng_free(dist->gslrand);
  free(dist);
}
