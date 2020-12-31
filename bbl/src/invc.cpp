#include <cmath>
#include <Rcpp.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_rng.h>
#include "bfgs.h"

const int Lbig=10;

using namespace std;

void invC(const vector<vector<short> > &ai, const vector<int> &frq,
          const vector<short> &L, double &E, double &lnz, vector<vector<double> > &h,
          vector<vector<vector<double> > > &J, double eps, double priorCount){

  int nsnp=L.size();
  int ndim=0;
  vector<vector<double> > f1(nsnp);
  bool naive= (eps==0);
  vector<vector<vector<double> > > f2;
  for(int i=0; i<nsnp; i++){
    vector<vector<double> > tmp;
    f12(i, ai, frq, f1[i], tmp, L, naive, priorCount);
    if(!naive) f2.push_back(tmp);
    ndim += L[i];
  }

  gsl_matrix *A;
  gsl_matrix *Ai;
  gsl_permutation *perm;

  if(eps>0){
    A=gsl_matrix_alloc(ndim,ndim);   // serial version using GSL
    Ai=gsl_matrix_alloc(ndim,ndim);
    perm=gsl_permutation_alloc(ndim);
  }

  double tr=0;
  for(int i=0;i<nsnp;i++){ 
    int Li=L[i];
    for(int l=0;l<Li;l++)
      tr+=f1[i][l]*(1-f1[i][l]);
  }
  tr/=ndim;

  if(eps>0){
    int idx=0;
    for(int i=0;i<nsnp;i++){
      int Li= L[i];
      for(int l0=0;l0<Li;l0++){
        int jdx = 0;
        for(int j=0;j<nsnp;j++){ 
          int Lj = L[j];
          for(int l1=0;l1<Lj;l1++){
            double x=eps*(f2[i][j][Lj*l0+l1]-f1[i][l0]*f1[j][l1]);
            if(eps==0.0) continue;
            if(i==j && l0==l1)
              x += (1-eps)*tr;
            gsl_matrix_set(A, idx, jdx++, x);
          }
        }
        idx++;
      }
    }
    int s;
    gsl_linalg_LU_decomp(A,perm,&s);
    gsl_linalg_LU_invert(A,perm,Ai);
  }

  h.resize(nsnp);
  J.resize(nsnp);
  lnz=0;
  int idx=0;
  for(int i=0;i<nsnp;i++){
    int Li= L[i];
    h[i].resize(Li);
    J[i].resize(nsnp);
    for(int j=0;j<nsnp;j++){
      int Lj = L[j];
      J[i][j].resize(Li*Lj);
    }
    double f=0.0;
    double s0=0.0;
    for(int l0=0; l0<Li; l0++) s0 += f1[i][l0];
    lnz += -log(1-s0);
    for(int l0=0;l0<Li;l0++){
      f=log(f1[i][l0]/(1.0-s0));
      if(eps>0.0){
        int jdx=0;
        for(int j=0;j<nsnp;j++){
          int Lj = L[j];
          for(int l1=0;l1<Lj;l1++){
            if(i!=j){
              double x=gsl_matrix_get(Ai, idx, jdx);
              J[i][j][Lj*l0+l1]=-x;
              f+=x*f1[j][l1];
              lnz += 0.5*x*f1[i][l0]*f1[j][l1];
            }
            jdx++;
          }
        }
      }
      h[i][l0]=f;
      idx++;
    }
  }
  E=0;
  int nind=ai.size();
  for(int k=0; k<nind; k++){
    for(int i=0; i<nsnp; i++){
      int a0=ai[k][i];
      if(a0==0) continue;
      E += h[i][a0-1];
      if(eps>0){
        for(int j=i+1;j<nsnp;j++){
          int a1=ai[k][j];
          if(a1==0) continue;
          E += J[i][j][L[j]*(a0-1)+a1-1];
        }
      }
    }
  }
  E = E/nind - lnz;

  if(eps>0){
    gsl_matrix_free(A);
    gsl_matrix_free(Ai);
    gsl_permutation_free(perm);
  }
}
