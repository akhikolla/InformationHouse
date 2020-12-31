/* stable/stable_koutrouvelis.h
 *
 * Koutrouvelis method for parameter estimation of alpha-stable
 * distributions. Based on [1] and MATLAB code available in [2].
 *
 * [1] Koutrouvelis, I. A. An Iterative Procedure for the Estimation
 *     of the Parameters of Stable Laws Communications in Statistics
 *     - Simulation and Computation, 1981, 10, 17-28
 * [2] Mark S. Veillete. Alpha-Stable Distributions in MATLAB
 *     http://math.bu.edu/people/mveillet/research.html
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
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_fit.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit.h>
#include "stable.h"

static inline double sign(double x) {
  return ((.0 < x) - (x < .0));
}

double var(const double * data, const int N) {
  double acum = .0, acum2 = .0;
  double v;

  int i;
  for(i=0;i<N;i++) {
    acum  += data[i];
    acum2 += data[i]*data[i];
  }

  v = (1.0/(N-1.0)) * (acum2 - (1.0/N)*acum*acum);
  return v;
}

int chooseK(double,int);

int chooseL(double,int);

void setcovYY(const double * t, int K, int N, double alpha, double beta, double gam, double **covYY);

void setcovZZ(const double * t, int K, int N, double alpha, double beta, double gam, double **covZZ);

double ecfRoot(const double * data, const int N);

void stable_samplecharfunc(const double* x, const unsigned int Nx,
             const double* t, const unsigned int Nt, gsl_complex *z);

gsl_complex stable_samplecharfunc_point(const double* x,
            const unsigned int N, double t);


int stable_fit_koutrouvelis(StableDist * dist, const double * data, const unsigned int N) {

  int maxiter = 10;
  double xTol = 0.01;

  double * s = NULL;
  gsl_complex * phi = NULL;
  gsl_complex * phi2 = NULL;

  int i;
  int stat;

  double alpha = dist->alpha;
  double beta  = dist->beta;
  double sigma = dist->sigma;
  double mu1   = dist->mu_1;

  double alphaold = alpha;
  double mu1old   = mu1;

  double diff     = .0;
  double diffbest = DBL_MAX;
  double diffmax  = .0;

  double alphabest = alpha;
  double betabest  = beta;
  double sigmabest = sigma;
  double mu1best   = mu1;

  double * t = NULL;
  double * w = NULL;
  double * y = NULL;
  double * p = NULL;

  double * t2 = NULL;
  double * w2 = NULL;
  double * y2 = NULL;
  double * p2 = NULL;

  gsl_matrix * X = NULL;
  gsl_matrix * covmat = NULL;
  gsl_vector * cvout = NULL;
  gsl_vector * ydat = NULL;
  gsl_vector * weights = NULL;

  double sumsq = .0;
  double c0 = .0, c1 = .0;
  double cov00 = .0, cov01 = .0, cov11 = .0;
  double sigmanew = 0;

  double **covYY = NULL;// **covZZ;
  double **covZZ = NULL;

  int iter = 0;
  int K = 0;
  int L = 0;
  int row, col;
  double t_0=0;
  double step_u=0;
  double estshift;

  gsl_multifit_linear_workspace * linws = NULL;

  if (sigma == 0) {
    sigma = sqrt(var(data,N));
	stable_setparams(dist,alpha,beta,sigma,mu1,1);
  }

  s = (double *) malloc(N*sizeof(double));
  for (i=0;i<N;i++) {
    s[i] = (data[i]-mu1)/sigma;
  }

  covmat = gsl_matrix_alloc(2,2);
  cvout = gsl_vector_alloc(2);

  for (iter = 0; iter<maxiter; iter ++) {
    if (iter <= 1) {
      K = chooseK(alpha,N);

      t = (double *) realloc(t,K*sizeof(double));
      p = (double *) realloc(p,K*sizeof(double));
      w = (double *) realloc(w,K*sizeof(double));
      y = (double *) realloc(y,K*sizeof(double));

      phi = (gsl_complex*) realloc(phi,K*sizeof(gsl_complex));

      for(i=0;i<K;i++) {
        t[i]=((double)i+1.0)*M_PI/25.0;
        w[i]=log(fabs(t[i]));
      }
    }

    if (iter==1) {
      covYY = (double **)malloc(K*sizeof(double *));
      covYY[0] = (double*)malloc(K*K*sizeof(double));
      for (row=0;row<K;row++) {
        covYY[row] = covYY[0] + row*K;
        for (col=0; col<K; col++) {
            covYY[row][col] = 0.0;
        }
      }
    }

    stable_samplecharfunc(s,N,t,K,phi);

    for(i=0;i<K;i++) {
      y[i] = log(-2.0*gsl_complex_logabs(phi[i]));
    }

    if (iter==0) { //Use ordinary least squares regression

      stat = gsl_fit_linear (w, 1, y, 1, K, &c0, &c1, &cov00, &cov01, &cov11, &sumsq);
      alpha = c1;
      sigmanew = pow(exp(c0)/2.0,1.0/alpha);
      sigma = sigma*sigmanew;

    }
    else { //Use weighted least squares regression

      setcovYY(t, K, N, alpha, beta, 1.0, covYY);
      for (i=0;i<K;i++) {
        p[i] = 1.0/covYY[i][i];
      }
      stat = gsl_fit_wlinear(w, 1, p, 1, y, 1, K, &c0, &c1, &cov00, &cov01, &cov11, &sumsq);
      alpha = c1;
      sigmanew = pow(exp(c0)/2.0,1.0/alpha);
      sigma = sigma*sigmanew;

    }


    /*************** rescale data ******************/

    for (i=0;i<N;i++) {
      s[i] = s[i]/sigmanew;
    }
    if (alpha<0) alpha = 0;
    if (alpha>2) alpha = 2;
    if (beta<-1) beta = -1;
    if (beta> 1) beta = 1;
    if (sigma<0) sigma = 0;


    /********** refine beta and mu **************/

    if (iter <= 1) {
      L = chooseL(alpha,N);

      t_0 = ecfRoot(s,N);
      step_u = M_PI/50;
      if (t_0/L<step_u) step_u = t_0/L;

      t2 = (double *) realloc(t2,L*sizeof(double));
      p2 = (double *) realloc(p2,L*sizeof(double));
      w2 = (double *) realloc(w2,L*sizeof(double));
      y2 = (double *) realloc(y2,L*sizeof(double));

      phi2 = (gsl_complex*)realloc(phi2,L*sizeof(gsl_complex));

      for(i=0;i<L;i++) {
        t2[i]=((double)i+1.0)*step_u;
        w2[i]=sign(t2[i])*pow(fabs(t2[i]),alpha);
      }
    }
    if (iter==0){

      linws = gsl_multifit_linear_alloc (L,2);
    }
    else if (iter ==1) {
      //Reserve covariance matrix mem
      covZZ = (double **)malloc(L*sizeof(double *));
      covZZ[0] = (double*)malloc(L*L*sizeof(double));
      for (row=0;row<L;row++) {
        covZZ[row] = covZZ[0] + row*L;
        for (col=0; col<L; col++) {
            covZZ[row][col] = 0.0;
        }
      }
      gsl_multifit_linear_free(linws);
      linws = gsl_multifit_linear_alloc (L,2);
    }


    stable_samplecharfunc(s,N,t2,L,phi2);


    for(i=0;i<L;i++) {
      y2[i] = gsl_complex_arg(phi2[i]);
    }

    X    = gsl_matrix_alloc(L,2);
    ydat = gsl_vector_alloc(L);

    for(row=0;row<L;row++) {
      gsl_matrix_set(X,row,0,t2[row]);
      gsl_matrix_set(X,row,1,w2[row]);
      gsl_vector_set(ydat,row,y2[row]);
    }

    if (iter==0) {

      stat = gsl_multifit_linear(X,ydat,cvout,covmat,&sumsq,linws);

    }
    else {

      setcovZZ(t2, L, N, alpha,  beta,  1.0, covZZ);
      weights = gsl_vector_alloc(L);
      for (i=0;i<L;i++) {
        p2[i] = 1.0/covZZ[i][i];
        gsl_vector_set(weights,i,p2[i]);
      }
      stat = gsl_multifit_wlinear(X,weights,ydat,cvout,covmat,&sumsq,linws);

    }
    if (alpha>1.98 || alpha < 0.05 || fabs(alpha-1.0)<0.05 || isnan(alpha)) {
        beta = 0.0;
    } else {
        beta = gsl_vector_get(cvout,1)/tan(alpha*M_PI*0.5);
    }
    estshift = gsl_vector_get(cvout,0);
    mu1    = mu1 + sigma * estshift;


    /*** remove estimated shift ***/
    for (i=0;i<N;i++) {
      s[i] = s[i] - estshift;
    }

    if (isnan(alpha) || isnan(beta) || isnan(sigma) || isnan(mu1))
    {
      iter++;
      break;
    }


    /** check for convergence or blow up**/
    diff = pow(alpha - alphaold,2) + pow(mu1 - mu1old,2);
    if (iter <= 2 && diff > diffmax) {
      diffmax = diff;
    }

#ifdef DEBUG
     Rprintf("blow up on iter %d\n",iter);
#endif

    if (fabs(diff) < diffbest) {
      alphabest = alpha;
      betabest  = beta;
      sigmabest = sigma;
      mu1best   = mu1;

      diffbest = diff;

      if (diff < xTol) {
        gsl_matrix_free(X);
        gsl_vector_free(ydat);
        if (iter>0) gsl_vector_free(weights);
        iter++;
        break;
      }
    }

    alphaold = alpha;
    mu1old   = mu1;

    gsl_matrix_free(X);
    gsl_vector_free(ydat);
    if (iter>0) {
      gsl_vector_free(weights);
    }
  }


  if (maxiter > 0 && iter >= 0) {
      alpha = alphabest;
      beta = betabest;
      sigma = sigmabest;
      mu1   = mu1best;
  }

  if (alpha<0) alpha = 0;
  if (alpha>2) alpha = 2;
  if (beta<-1) beta = -1;
  if (beta> 1) beta = 1;
  if (sigma<0) sigma = 0;

  stable_setparams(dist,alpha,beta,sigma,mu1,1);

  free(t);
  free(p);
  free(w);
  free(y);
  free(phi);
  free(t2);
  free(p2);
  free(w2);
  free(y2);
  free(phi2);
  gsl_matrix_free(covmat);
  gsl_vector_free(cvout);
  free(s);

  if (iter>1) {
    free(covYY[0]);
    free(covYY);
    free(covZZ[0]);
    free(covZZ);
  }

  gsl_multifit_linear_free(linws);

  (void)stat;  // Just to shut off unused-but-set-variable warning (not elegant but portable)
  return 0;
}  // end stable_fit_koutrouvelis

int chooseK(double alpha, int N) {
  double a[]={1.9,1.5,1.3,1.1,.9,.7,.5,.3};
  int n[] = {200, 800, 1600};

  double Kmat[8][3]={{ 9, 9, 9},
                     {11,11,11},
                     {22,16,14},
                     {24,18,15},
                     {28,22,18},
                     {30,24,20},
                     {86,68,56},
                     {134,124,118}};

  int i,j;
  double xi;
  double xj;
  double Kp;

  if (alpha < 0.3) alpha = 0.3;
  if (alpha > 1.9) alpha = 1.9;
  if (N<200) N=200;
  if (N>1600) N=1600;

  i=1;
  j=1;
  while (i<7 && a[i]>=alpha) {
    ++i;
  }
  while (j<2 && n[j]<=N) {
    ++j;
  }

  xi = 1.0 - (alpha-a[i])/(a[i-1]-a[i]);
  xj = 1.0 - (double)(n[j]-N)/(double)(n[j]-n[j-1]);

  /* bilinear interpolation */
  Kp =       xj * (xi * Kmat[i][ j ] + (1.0-xi) * Kmat[i-1][ j ] ) +
       (1.0-xj) * (xi * Kmat[i][j-1] + (1.0-xi) * Kmat[i-1][j-1] );

  Kp = floor(Kp+.5);

  return (int)Kp;
}

int chooseL(double alpha, int N) {
  double a[]={1.9,1.5,1.1,.9,.7,.5,.3};
  int n[] = {200, 800, 1600};

  double Lmat[7][3]={{ 9,10,11},
                     {12,14,15},
                     {16,18,17},
                     {14,14,14},
                     {24,16,16},
                     {40,38,36},
                     {70,68,66}};

  int i,j;
  double xi, xj;
  double Lp;

  if (alpha < 0.3) alpha = 0.3;
  if (alpha > 1.9) alpha = 1.9;
  if (N<200) N=200;
  if (N>1600) N=1600;

  i=1;
  j=1;
  while (i<6 && a[i]>=alpha) {
    ++i;
  }
  while (j<2 && n[j]<=N) {
    ++j;
  }

  xi = 1.0 - (alpha-a[i])/(a[i-1]-a[i]);
  xj = 1.0 - (double)(n[j]-N)/(double)(n[j]-n[j-1]);

  /* bilinear interpolation */
  Lp =       xj * (xi * Lmat[i][ j ] + (1.0-xi) * Lmat[i-1][ j ] ) +
       (1.0-xj) * (xi * Lmat[i][j-1] + (1.0-xi) * Lmat[i-1][j-1] );

  Lp = floor(Lp + .5);

  return (int)Lp;
}

void setcovYY(const double * t, int K, int N, double alpha, double beta, double gam, double **covYY) {
  double w = tan(alpha*M_PI_2);
  double calpha = pow(gam,alpha);
  int    j = 0, k = 0;

  double * talpha = (double *)malloc(K*sizeof(double));

  for (j=0;j<K;j++) {
    talpha[j]=pow(fabs(t[j]),alpha);
  }

  double tj,tja,tk,tka;
  double tjmtka,tjptka,stj,stk;
  double A,B,D,E;//F,G,H;
  for (j=0;j<K;j++) {
    tj = t[j];
    tja = talpha[j];
    stj = sign(tj);
    for (k=0;k<K;k++) {
      tk  = t[k];
      tka = talpha[k];
      stk = sign(tk);
      tjmtka = pow(fabs(tj-tk),alpha);
      tjptka = pow(fabs(tj+tk),alpha);

      A =  calpha * (tja + tka - tjmtka);
      B =  calpha * beta * (-tja * stj * w + tka * stk * w + tjmtka * sign(tj - tk) * w );
      D =  calpha * (tja + tka - tjptka);
      E =  calpha * beta * ( tja * stj * w + tja * stk * w - tjptka * sign(tj + tk) * w );

      covYY[j][k] = (exp(A)*cos(B)+exp(D)*cos(E)-2.0) / ( 2.0 * N * pow(gam,2.0*alpha)*pow(fabs(tj*tk),alpha) );

    }
  }

  free(talpha);

  return;
}
void setcovZZ(const double * t, int K, int N, double alpha, double beta, double gam, double **covZZ) {
  double w = tan(alpha*M_PI_2);
  double calpha = pow(gam,alpha);
  int    j = 0, k = 0;

  double * talpha = (double *)malloc(K*sizeof(double));

  for (j=0;j<K;j++) {
    talpha[j]=pow(fabs(t[j]),alpha);
  }

  double tj,tja,tk,tka;
  double tjmtka,tjptka,stj,stk;
  double B,E,F,G,H;
  for (j=0;j<K;j++) {
    tj = t[j];
    tja = talpha[j];
    stj = sign(tj);
    for (k=0;k<K;k++) {
      tk  = t[k];
      tka = talpha[k];
      stk = sign(tk);
      tjmtka = pow(fabs(tj-tk),alpha);
      tjptka = pow(fabs(tj+tk),alpha);


      B =  calpha * beta * (-tja * stj * w + tka * stk * w + tjmtka * sign(tj - tk) * w );

      E =  calpha * beta * ( tja * stj * w + tja * stk * w - tjptka * sign(tj + tk) * w );
      F =  calpha * (tja + tka);
      G = -calpha * tjmtka;
      H = -calpha * tjptka;

      covZZ[j][k] =  exp(F)*( exp(G)*cos(B)-exp(H)*cos(E) ) /(2.0*N);
    }
  }

  free(talpha);

  return;
}

double ecfRoot(const double * data, int N) {

  double m=0;
  int i;
  for(i=0;i<N;i++) {
    m += fabs(data[i]);
  }
  m = m/N;

  double t = 0;
  gsl_complex val = stable_samplecharfunc_point(data, N, t);
  int iter=0;
  while (iter<10000 && fabs(GSL_REAL(val))>1e-3) {
    t += GSL_REAL(val)/m;
    val = stable_samplecharfunc_point(data, N, t);
  }

  return t;
}

gsl_complex stable_samplecharfunc_point(const double* x,
            const unsigned int N, double t)
{
  unsigned int i;
  double zr=0.0,zi=0.0;
  gsl_complex z;

  for(i=0;i<N;i++)
   {
      zr+=cos(t*x[i]);
      zi+=sin(t*x[i]);
    }

  GSL_SET_COMPLEX(&z,zr/N,zi/N);

  return z;
}

void stable_samplecharfunc(const double* x, const unsigned int Nx,
             const double* t, const unsigned int Nt, gsl_complex *z)
{
  unsigned int ix,it;
  double zr=0.0,zi=0.0,w;

  for(it=0;it<Nt;it++)
   {
     w = t[it];
     zr = .0;
     zi = .0;
     for(ix=0;ix<Nx;ix++)
       {
         zr+=cos(w*x[ix]);
         zi+=sin(w*x[ix]);
       }

     GSL_SET_COMPLEX(&z[it],zr/Nx,zi/Nx);
   }
  return;
}


