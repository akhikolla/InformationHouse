/* stable/stable.h
 *
 * Main header file of Libstable. Contains all declarations of the
 * usable functions in the library.
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

#ifndef _stable_H_
#define _stable_H_

#include <stdlib.h>
#include <string.h>
#include <stdio.h>

/*
#ifdef __WIN32
#include <windows.h>
#endif
*/

#include <gsl/gsl_integration.h>
// #include <gsl/gsl_rng.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_vector.h>

#define TINY 1e-50
#define EPS 2.2204460492503131E-16

#define max(a,b) (((a) > (b)) ? (a) : (b))
#define min(a,b) (((a) < (b)) ? (a) : (b))

/******************************************************************************/
/*          Library parameters                                                */
/******************************************************************************/

extern FILE * FLOG;            // Log file (optional)
extern FILE * FINTEG;          // Integrand evaluations output file (debug purposes)

extern unsigned short THREADS;    // threads of execution (0 => total available)
extern unsigned short IT_MAX;     // Maximum # of iterations in quadrature methods
extern unsigned short SUBS;       // # of integration subintervals
extern unsigned short METHOD1;    // Integration method on main subinterval
extern unsigned short METHOD2;    // Integration method on second subinterval
extern unsigned short METHOD3;    // Integration method on third subinterval

extern unsigned short INV_MAXITER; // Maximum # of iterations inversion method

extern double relTOL;     // Relative error tolerance
extern double absTOL;     // Absolut error tolerance
//extern double FACTOR;   //
extern double ALPHA_TH;    // Alpha threshold
extern double BETA_TH;    // Beta threshold
extern double EXP_MAX;    // Exponent maximum value
extern double XXI_TH;     // Zeta threshold
extern double THETA_TH;   // Theta threshold

extern double AUX1; // Auxiliary values
extern double AUX2;

#ifdef DEBUG
extern unsigned int integ_eval; // # of integrand evaluations
#endif

/******************************************************************************/
/******************************************************************************/
// Particular cases
enum
  {
    NOVALID =-1,
    STABLE,
    ALPHA_1,
    GAUSS ,
    CAUCHY,
    LEVY,
    STABLE_B1,
    ALPHA_1_B1
  };

// Function to evaluate
enum
  {
    CDF,
    PDF
  };

// Quadrature methods
enum
  {
    STABLE_QAG2 = 0,
    STABLE_QUADSTEP,
    STABLE_QROMBPOL,
    STABLE_QROMBRAT,
    STABLE_QNG,
    STABLE_QAG1,
    STABLE_QAG5,
    STABLE_VECT
  };


/************************************************************************
 ************************************************************************
 * Scalar methods                                                       *
 ************************************************************************
 ************************************************************************/
/* Scalar methods are those that, for each thread of execution, obtain a
   single evaluation of the desired function, at a single point.
*/

/******************************************************************************/
/*    Stable distribution structure.                                          */
/******************************************************************************/
struct StableDistStruct
  {
    /* Parameters:
    0-parametrization describen in Nolan, 1997 is employed by default
        alpha : stability index
        beta : skewness parameter
        scale: scale parameter
        mu_0 : 0-parametrization location parameter
        mu_1 : correspondig 1-parametrization location parameter    */
    double alpha;
    double beta;
    double sigma;
    double mu_0;
    double mu_1;


    /* Particular cases indicator (Gauss, Cauchy, Levy distribution, alpha==1, etc.) */
    int ZONE;

    /* Pointers to pdf and cdf evaluation functions */
    double(*stable_pdf_point)(struct StableDistStruct *, const double, double *);
    double(*stable_cdf_point)(struct StableDistStruct *, const double, double *);

    /* Precalculated values. */
    double alphainvalpha1;  /* alpha/(alpha-1)*/
    double xi;            /* -beta*tan(alpha*pi/2)*/
    double theta0;        /* 1/alpha*atan(beta*(tan(alpha*pi/2))=atan(-xi)/alpha;*/
    double c1, c2_part, c3;  /* additive and multiplicative constants*/
    double k1;     /* cos(alpha*theta0)^(1/(alpha-1)) = (1+xi^2)^(-0.5/(alpha-1));*/
    double S;     /* (1+xi^2)^(1/(2*alpha));*/
    double Vbeta1; /*pow(1/dist->alpha,dist->alphainvalpha1) *
                     (dist->alpha-1)*pow(-cos(dist->alpha*PI_2),1/(dist->alpha-1))*/

    /* These ones change from point to point of evaluation */
    double theta0_; /* theta0_ = +-theta0 */
    double beta_;
    double xxipow;  /* (x-xi)^(alpha/(alpha-1))*/

    /* gsl integration workspace */
    gsl_integration_workspace * gslworkspace;

    /* gsl random numbers generator */
    // gsl_rng * gslrand;
  };

typedef struct StableDistStruct StableDist;
/******************************************************************************/
/*        Auxiliary functions                                                 */
/******************************************************************************/
unsigned int stable_get_THREADS();
void   stable_set_THREADS(unsigned int threads);

unsigned int stable_get_IT_MAX();
void   stable_set_IT_MAX(unsigned int itmax);

unsigned int stable_get_INV_MAXITER();
void   stable_set_INV_MAXITER(unsigned int invmaxiter);

int stable_get_METHOD1();
void stable_set_METHOD1(int method);

int stable_get_METHOD2();
void stable_set_METHOD2(int method);

int stable_get_METHOD3();
void stable_set_METHOD3(int method);

double stable_get_relTOL();
void   stable_set_relTOL(double reltol);

double stable_get_absTOL();
void   stable_set_absTOL(double abstol);

double stable_get_ALPHA_TH();
void   stable_set_ALPHA_TH(double alphath);

double stable_get_BETA_TH();
void   stable_set_BETA_TH(double betath);

double stable_get_XXI_TH();
void   stable_set_XXI_TH(double xxith);

double stable_get_THETA_TH();
void   stable_set_THETA_TH(double thetath);

FILE * stable_get_FINTEG();
FILE * stable_set_FINTEG(char * filename);

FILE * stable_get_FLOG();
FILE * stable_set_FLOG(char * filename);


StableDist *stable_create(double alpha, double beta, double sigma, double mu,
                          int parametrization);

StableDist *stable_copy(StableDist *src_dist);

void stable_free(StableDist *dist);

int stable_setparams(StableDist *dist,
                     double alpha, double beta, double sigma, double mu,
                     int parametrization);

int stable_checkparams(double alpha, double beta, double sigma, double mu,
                       int parametrization);

void error_handler (const char * reason, const char * file,
               int line, int gsl_errno);

/******************************************************************************/
/*   PDF in particular cases                                                  */
/******************************************************************************/

double stable_pdf_point_GAUSS(StableDist *dist, const double x, double *err);

double stable_pdf_point_CAUCHY(StableDist *dist, const double x, double *err);

double stable_pdf_point_LEVY(StableDist *dist, const double x, double *err);

/******************************************************************************/
/*   PDF otherwise                                                            */
/******************************************************************************/

double stable_pdf_point_STABLE(StableDist *dist, const double x, double *err);

double stable_pdf_point_ALPHA_1(StableDist *dist, const double x, double *err);

double stable_pdf_point(StableDist *dist, const double x, double *err);

void stable_pdf(StableDist *dist, const double* x, const unsigned int Nx,
                double *pdf, double *err);

/******************************************************************************/
/*   PDF integrand functions                                                  */
/******************************************************************************/

double stable_pdf_g(double theta, void *dist);
double stable_pdf_g_aux1(double theta, void *args);
double stable_pdf_g_aux2(double theta, void *args);

/******************************************************************************/
/*   CDF in particular cases                                                  */
/******************************************************************************/

double stable_cdf_point_GAUSS(StableDist *dist, const double x, double *err);

double stable_cdf_point_CAUCHY(StableDist *dist, const double x, double *err);

double stable_cdf_point_LEVY(StableDist *dist, const double x, double *err);

/******************************************************************************/
/*   CDF otherwise                                                            */
/******************************************************************************/

double stable_cdf_point_STABLE(StableDist *dist, const double x, double *err);

double stable_cdf_point_ALPHA_1(StableDist *dist, const double x, double *err);

double stable_cdf_point(StableDist *dist, const double x, double *err);

void stable_cdf(StableDist *dist, const double* x, const unsigned int Nx,
                double *cdf, double *err);

/******************************************************************************/
/*   CDF integrad functions                                                   */
/******************************************************************************/

double stable_cdf_g(double theta, void *dist);

/******************************************************************************/
/*   CDF^{-1} (quantiles)                                                     */
/******************************************************************************/

double stable_q_point(StableDist * dist, const double q, double * err);
void   stable_q(StableDist *dist, const double* q, const unsigned int Nq,
                double * inv, double * err);

/************************************************************************
 ************************************************************************
 * Vectorial methods                                                    *
 ************************************************************************
 ************************************************************************/
/*
  Alternative non-parallelized methods of evaluation have been implemented.
  These methods exploit the fact that some calculations are shared between
  different points of evaluation when evaluating the PDF or CDF, so these
  calculations can be realized just once.

  The performance achieved is high, sometimes comparable with parallel
  methods when little precision is required. However, achievable precision
  with these methods is low and non desired behavior of the PDF and CDF
  evaluation is observed.
*/

/* Stable distribution structure for vectorial methods*/

typedef struct
  {
/* Parameters:
    0-parametrization describen in Nolan, 1997 is employed by default
        alpha : stability index
        beta : skewness parameter
        scale: scale parameter
        mu_0 : 0-parametrization location parameter
        mu_1 : correspondig 1-parametrization location parameter    */
    double alpha;
    double beta;
    double sigma;
    double mu_0;
    double mu_1;

    /* Particular cases indicator (Gauss, Cauchy, Levy distribution, alpha==1, etc.) */
    int ZONE;

    /* Precalculated values */
    double alphainvalpha1;     /* alpha/(alpha-1)*/
    double xi;               /* -beta*tan(alpha*pi/2)*/
    double theta0;           /* 1/alpha*atan(beta*(tan(alpha*pi/2))=atan(-xi)/alpha;*/
    double c1, c2_part, c3;  /* additive and multiplicative constants*/
    double k1;               /* cos(alpha*theta0)^(1/(alpha-1)) = (1+xi^2)^(-0.5/(alpha-1));*/
    double S;                /* (1+xi^2)^(1/(2*alpha));*/
    double Vbeta1;           /*pow(1/dist->alpha,dist->alphainvalpha1) *
                                   (dist->alpha-1)*pow(-cos(dist->alpha*PI_2),1/(dist->alpha-1))*/

     /* These ones change from point to point of evaluation */
    double theta0_; /* theta0_ = +-theta0 */
    double beta_;
    double *xxipow;  /* (x-xi)^(alpha/(alpha-1))*/

    gsl_integration_workspace * gslworkspace;
    // gsl_rng * gslrand;
  }
StableDistV;

typedef struct
  {
    double (*ptr_funcion)(StableDist *dist, const double x, double *err);
    StableDist *dist;
    const double *x;
    int Nx;
    double *pdf;
    double *err;
  }
StableArgsPdf;

typedef struct
  {
    double (*ptr_funcion)(StableDist *dist, const double x, double *err);
    StableDist *dist;
    const double *x;
    int Nx;
    double *cdf;
    double *err;
  }
StableArgsCdf;

/************************************************************************
 ************************************************************************
 * Parameter estimation                                                 *
 ************************************************************************
 ************************************************************************/


/******************************************************************************/
/*        Parameter estimation structure                                      */
/******************************************************************************/
typedef struct
  {
    StableDist *dist;
    double *data;
    unsigned int length;
    double nu_c;
    double nu_z;
  }
stable_like_params;

/* Estimation functions */

void stable_fit_init(StableDist *dist, const double *data,
             const unsigned int length,  double *nu_c,double *nu_z);

int stable_fit_koutrouvelis(StableDist *dist, const double *data, const unsigned int length);

int stable_fit(StableDist *dist, const double *data, const unsigned int length);

int stable_fit_mle(StableDist *dist, const double *data, const unsigned int length);

int stable_fit_mle2d(StableDist *dist, const double *data, const unsigned int length);

int stable_fit_whole(StableDist *dist, const double *data, const unsigned int length);

/* Auxiliary functions */

gsl_complex stable_samplecharfunc_point(const double* x,
             const unsigned int N, double t);

void stable_samplecharfunc(const double* x, const unsigned int Nx,
             const double* t, const unsigned int Nt, gsl_complex * z);

void stable_fft(double *data, const unsigned int length, double * y);

double stable_loglikelihood(StableDist *dist, double *data, const unsigned int length);

//stable_like_params

int stable_fit_iter_whole(StableDist *dist, const double * data, const unsigned int length);

int stable_fit_iter(StableDist *dist, const double * data,
             const unsigned int length, const double nu_c, const double nu_z);

double stable_loglike_p(stable_like_params *params);

double stable_minusloglikelihood(const gsl_vector * theta, void * p);

/************************************************************************
 ************************************************************************
 * Random numbers generation                                            *
 ************************************************************************
 ************************************************************************/
void stable_rnd(StableDist *dist, double* rnd, unsigned int n);

// static inline double stable_rnd_point(StableDist *dist);

// void stable_rnd_seed(StableDist * dist, unsigned long int s);
#endif //stable_H
