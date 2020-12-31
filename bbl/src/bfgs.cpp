#include <cmath>
#include <string>
#include <vector>
#include <iostream>
#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_multimin.h>
#include <Rcpp.h>
#include "bfgs.h"

using namespace std;

double pan3(vector<double> &peff, int nsnp, int i0, 
          const vector<short> &L, const vector<short> &ci, vector<double> h1, 
          const vector<vector<double> > &J1, bool naive, bool lzhalf){

  peff.resize(L[i0]);
  vector<double> peff2 = peff;
  for(int a=0;a<L[i0];a++){
    double e= h1[a];
    double e2=e;
    if(!naive){
      for(int j=0;j< nsnp;j++){
        if(j==i0) continue;
        int b=ci[j];
        if(b==0) continue;
        e+= J1[j][L[j]*a+b-1];
        if(lzhalf) e2+= J1[j][L[j]*a+b-1]/2.0;
      }
    }
    peff[a]=e;
    if(lzhalf) peff2[a]=e2;
  }
  double max=0;
  double max2=0;
  for(int a=0;a<L[i0];a++){
    if(peff[a]>max) max=peff[a];
    if(lzhalf && (peff2[a] > max2)) max2=peff2[a];
  }
  double z = exp(-max);
  double z2= exp(-max2);
  for(int a=0;a<L[i0];a++){
    peff[a] = exp(peff[a]-max);
    z+=peff[a];
    if(lzhalf){
      peff2[a] = exp(peff2[a]-max2);
      z2+=peff2[a];
    }
  }
  for(int a=0;a<L[i0];a++){
    peff[a]/=z;
    if(lzhalf) peff2[a]/=z2;
  }
  
  double lz=0;
  if(lzhalf) lz += log(z2)+max2;
  else lz += log(z)+max;
  
  return lz;
}

double pan2(int nsnp, int i0, const vector<short> &L, 
            const vector<short> &ci, 
            const vector<double> &h1, const vector<vector<double> > &J1, 
            double &lz, bool naive, bool lzhalf){

  vector<double> peff(L[i0]);
  lz = pan3(peff, nsnp, i0, L, ci, h1, J1, naive, lzhalf);

  int a0=ci[i0];
  if(a0>0)
    return peff[a0-1];
  double p=1;
  for(int l=0; l<L[i0]; l++)
    p -= peff[l];
   
  return p;
}

double lnl_psl(const gsl_vector *v,void *params){  // evaluates log likelihood

  double ln;
  Param *par=(Param *)params;
  int i0=par->i0;
  vector<short> L=par->L;
  double lambda=par->lambda;
  double lambdah=par->lambdah;
  int nsnp=(par->ai)[0].size();

  vector<double> h1(L[i0]);
  vector<vector<double> > J1(nsnp);
  
  if(!par->naive)
    for(int i=0; i<nsnp; i++) J1[i].resize(L[i0]*L[i]);

  int m=0;
  for(int l0=0;l0<L[i0];l0++){
    h1[l0]=gsl_vector_get(v,m++);
    for(int i=0; i<nsnp; i++){ 
      if(i==i0 || par->naive) continue;
      if(!(par->qj)[i]) continue;
      for(int l1=0;l1<L[i];l1++)
        J1[i][L[i]*l0+l1]=gsl_vector_get(v,m++);
    }
  }

  int nind=int((par->ai).size());
  ln=0;
  par->lzp = 0;
  double wsum = 0;
  for(int n=0;n<nind;n++){
    double lz=0;
    double p=pan2(nsnp,i0,L,(par->ai)[n],h1,J1,lz, par->naive, par->lzhalf);
    ln += -log(p)*(par->frq)[n];
    par->lzp += lz;
    wsum += (par->frq)[n];
  }
  ln /= wsum;
  par->lzp /= wsum;

  for(int l=0;l<L[i0];l++)
    ln+= lambdah*h1[l]*h1[l]/2;
  
  if(par->naive) return ln;
  
  for(int i=0; i<nsnp; i++){
    if(i==i0) continue;
    if(!(par->qj)[i]) continue;
    for(int l=0;l<L[i0]*L[i];l++)
      ln+=lambda*J1[i][l]*J1[i][l]/2;
  }
  return ln;
}


void dlnl_psl(const gsl_vector *v,void *params,gsl_vector *df){   // first derivatives

  Param *par=(Param *)params;
  vector<short> L=par->L;
  int nsnp=(par->ai)[0].size();
  double lambda=par->lambda;
  double lambdah=par->lambdah;
  int i0=par->i0;

  vector<double> s1(L[i0]);
  vector<vector<double> > s2(nsnp);

  vector<double> h1(L[i0]);
  vector<vector<double> > J1(nsnp);
  
  if(!par->naive)
    for(int i=0; i<nsnp; i++) J1[i].resize(L[i0]*L[i]);

  int m=0;
  for(int l0=0;l0<L[i0];l0++){
    h1[l0]=gsl_vector_get(v,m++);
    if(par->naive) continue;
    for(int i=0; i<nsnp; i++){ 
      if(i==i0) continue;
      if(!(par->qj)[i]) continue;
      for(int l1=0; l1<L[i]; l1++)
        J1[i][L[i]*l0+l1]=gsl_vector_get(v,m++);
    }
  }

  int nind=int((par->ai).size());

  for(int l=0;l<L[i0];l++) s1[l]=0;
  if(!par->naive){
    for(int i=0; i<nsnp; i++){
      s2[i].resize(L[i0]*L[i]);
      for(int l=0;l<L[i0]*L[i];l++) s2[i][l]=0;
    }
  }

  double wsum = 0.0;
  for(int k=0;k<nind;k++)
    wsum += (par->frq)[k];
  
  for(int k=0;k<nind;k++){
    vector<double> peff(L[i0]);
    double lz=0;
    pan3(peff, nsnp, i0, L, (par->ai)[k], h1, J1, par->naive, par->lzhalf);
    for(int l0=0;l0<L[i0];l0++){
      double f=peff[l0]*(par->frq)[k]/wsum;
      s1[l0]+= f;
      if(par->naive) continue;
      for(int j=0;j<nsnp;j++){
        if(j==i0) continue;
        if(!(par->qj)[j]) continue;
        short a=(par->ai)[k][j];
        if(a==0) continue;
        s2[j][L[j]*l0+a-1] += f;
      }
    }
  }

  for(int l0=0;l0<L[i0];l0++){
    s1[l0] += lambdah*h1[l0] - (par->f1)[l0];
    if(par->naive) continue;
    for(int j=0;j<nsnp;j++){
      if(j==i0) continue;
      if(!(par->qj)[j]) continue;
      for(int l1=0;l1<L[j];l1++)
        s2[j][L[j]*l0+l1]+=-(par->f2)[j][L[j]*l0+l1]
             +lambda*J1[j][L[j]*l0+l1];
    }
  }

  m=0;
  for(int l0=0;l0<L[i0];l0++){
    gsl_vector_set(df,m++,s1[l0]);
    if(par->naive) continue;
    for(int i=0; i<nsnp; i++){ 
      if(i==i0) continue;
      if(!(par->qj)[i]) continue;
      for(int l1=0;l1<L[i];l1++)
        gsl_vector_set(df,m++,s2[i][L[i]*l0+l1]);
    }
  }
}

void ln_dln_psl(const gsl_vector *x,void *params,double *f,gsl_vector *df){

  *f=lnl_psl(x,params);
  dlnl_psl(x,params,df);

}

double lpr_psl(int i0, const vector<vector<short> > &ai, 
               const vector<int> &frq, const vector<bool> &qj,
               const vector<short> &L, double lambda,
               double lambdah, vector<double> &h, 
               vector<vector<double> > &J, int nprint, 
               unsigned int Imax, double Tol, int verbose, double &lz, 
               bool naive, bool &failed, bool lzhalf){

  size_t iter=0;
  int status;

  int nind=int(ai.size());
  int nsnp=int(ai[0].size());
  vector<double> f1(nsnp);
  vector<vector<double> > f2(nsnp);

  f12(i0, ai, frq, f1, f2, L, naive, 0);

  const gsl_multimin_fdfminimizer_type *T;
  gsl_multimin_fdfminimizer *s;

  gsl_vector *x;
  gsl_multimin_function_fdf my_func;

  int ndim = L[i0];
  if(!naive){
    for(int i=0; i<nsnp; i++)
      if(i!=i0 && qj[i]) ndim += L[i0]*L[i];
  }
  
  my_func.n=ndim;         
  my_func.f=lnl_psl;
  my_func.df=dlnl_psl;
  my_func.fdf=ln_dln_psl;

  x=gsl_vector_alloc(ndim);
  T=gsl_multimin_fdfminimizer_vector_bfgs2;  // BFGS2 optimizer
  s=gsl_multimin_fdfminimizer_alloc(T,ndim);

  Param par={i0, ai, frq, qj, L, lambda, lambdah, f1, f2, lz, naive, lzhalf};

  my_func.params=&par;
  gsl_vector_set_zero(x);  // initial guess

  gsl_multimin_fdfminimizer_set(s,&my_func,x,0.1,0.1);
  iter=0;
  do{
      iter++;
      status=gsl_multimin_fdfminimizer_iterate(s);
      if(iter%nprint==0 && verbose>1)
        Rcpp::Rcout << "  iteration # " << iter << ": " << s->f << endl;
      if(status){
        Rcpp::Rcerr << " GSL status code " << status << endl;
        failed=true;
        break;
      }
      status=gsl_multimin_test_gradient(s->gradient,Tol);
  }while(status==GSL_CONTINUE && iter< Imax); 
  if(iter==Imax)
    Rcpp::Rcerr << "BFGS2 iteration failed to converge after " 
         << Imax << " iterations\n";
    if(verbose > 0) Rcpp::Rcout << "  Predictor " << i0+1 
         << ": " << iter
         << " iterations, likelihood = " << s->f << endl;

  h.resize(L[i0]);
  J.resize(nsnp);
  double min=0;
  for(int i=0; i<nsnp; i++) J[i].resize(L[i0]*L[i]);
  int m=0;
  for(int l0=0;l0<L[i0];l0++){
    h[l0]=gsl_vector_get(s->x,m++);
    if(naive) continue;
    for(int i=0; i<nsnp; i++) for(int l1=0;l1<L[i];l1++){
      if(i==i0 || !qj[i]) J[i][L[i]*l0+l1]=0;
      else 
        J[i][L[i]*l0+l1]=gsl_vector_get(s->x,m++);
    }
  }
  min=-nind*(s->f);

  gsl_multimin_fdfminimizer_free(s);
  gsl_vector_free(x);

  return min;
}


void f12(int i0, const vector<vector<short> > &si, 
         const vector<int> &frq, vector<double> &f1, vector<vector<double> > &f2, 
         const vector<short> &L, bool naive, double pcount){

  int n = si.size();
  int m = si[0].size();
  int Li0 = L[i0];  // upper bound for x
  f1.resize(Li0);
  f2.resize(m);

  for(int l=0; l<Li0; l++)
    f1[l] = pcount/(1+Li0);
  
  if(!naive){
    for(int i=0; i<m; i++){
      int Li = L[i];
      f2[i].resize(Li0*Li);
      for(int l0=0; l0<Li0; l0++) for(int l1=0; l1<Li; l1++){
        int id = Li*l0+l1;
        f2[i][id]= (i==i0 ? pcount/(Li0+1) : pcount/(Li0+1)/(Li+1));
      }
    }
  }
  
  double wsum=0.0;
  for(int k=0; k<n; k++){
    wsum += frq[k];
    short a=si[k][i0];
    if(a==0) continue;
    f1[a-1] += frq[k];
    if(naive) continue;
    for(int j=0; j<m; j++){
      if(j==i0) continue;
      short b=si[k][j];
      if(b==0) continue;
      f2[j][L[j]*(a-1)+b-1] += frq[k];
    }
  }
  
  double nc = wsum + pcount;
  for(int l0=0; l0<f1.size(); l0++){ 
    f1[l0]/=nc;
    if(naive) continue;
    for(int j=0; j<m; j++){
      int Lj = L[j];
      if(i0==j){
        for(int l1=0; l1<Lj; l1++){
          if(l0==l1)
            f2[j][Lj*l0+l1]=f1[l0];
          else
            f2[j][Lj*l0+l1]=0;
        }
      }
      else{
        for(int l1=0; l1<Lj; l1++)
          f2[j][Lj*l0+l1]/=nc;
      }
    }
  }
}
