
#include "MAVEfast.h"
#include "MAVEfast.h"
#include "MAVEfast_emxutil.h"
#include <vector>
#include <map>
#include <stdio.h>
#include <algorithm>
#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace std;


// [[Rcpp::export]]
List MAVEfastCpp(NumericMatrix x,NumericMatrix y,CharacterVector method, double max_dim, double screen) {

  emxArray_real_T *emxArray_x,*emxArray_y,*emxArray_BB,*emxArray_ky,*emxArray_BBvs,*emxArray_idx,*emxArray_C;
  emxArray_char_T *emxArray_method;
  std::string std_method = Rcpp::as<std::string>(method);
  emxInit_real_T(&emxArray_x,2);
  emxInit_real_T(&emxArray_y,2);
  emxInit_real_T(&emxArray_BB,3);
  emxInit_real_T(&emxArray_ky,2);
  emxInit_char_T(&emxArray_method,2);
  emxInit_real_T(&emxArray_BBvs,3);
  emxInit_real_T(&emxArray_idx,2);
  emxInit_real_T(&emxArray_C,2);
  emxArray_x->size[0]=x.nrow();
  emxArray_x->size[1]=x.ncol();
  emxArray_y->size[0]=y.nrow();
  emxArray_y->size[1]=y.ncol();
  emxArray_method->size[0]=1;
  emxArray_method->size[1] = std_method.size();

  emxEnsureCapacity((emxArray__common *)emxArray_x, 0, (int)sizeof(double));
  emxEnsureCapacity((emxArray__common *)emxArray_y, 0, (int)sizeof(double));
  emxEnsureCapacity((emxArray__common *)emxArray_method, 0, (int)sizeof(char));

  for(int i=0;i<x.size();++i){
    emxArray_x->data[i] = x[i];
  }
  //for(int i=0;i<5;++i) printf("%.3lf ",x[i]); printf("\n");
  for(int i=0;i<y.size();++i){
    emxArray_y->data[i] = y[i];
  }
  for(int i=0;i<std_method.size();++i){
    emxArray_method->data[i] = std_method[i];
  }

  MAVEfast(emxArray_x, emxArray_y, emxArray_method, max_dim, screen, emxArray_BB, emxArray_ky, emxArray_BBvs, emxArray_idx,emxArray_C);

  NumericVector nx(emxArray_x->size[0]*emxArray_x->size[1]);
  NumericVector ky(emxArray_ky->size[0]*emxArray_ky->size[1]);
  NumericVector BB(emxArray_BB->size[0]*emxArray_BB->size[1]*emxArray_BB->size[2]);
  NumericVector BBvs(emxArray_BBvs->size[0]*emxArray_BBvs->size[1]*emxArray_BBvs->size[2]);
  NumericVector idx(emxArray_idx->size[0]*emxArray_idx->size[1]);
  NumericVector C(emxArray_C->size[0]*emxArray_C->size[1]);

  for(int i=0;i<emxArray_x->size[0]*emxArray_x->size[1];++i){
    nx[i]=emxArray_x->data[i];
  }

  for(int i=0;i<emxArray_BB->size[0]*emxArray_BB->size[1]*emxArray_BB->size[2];++i){
    BB[i]=emxArray_BB->data[i];
  }

    for(int i=0;i<emxArray_BBvs->size[0]*emxArray_BBvs->size[1]*emxArray_BBvs->size[2];++i){
        BBvs[i]=emxArray_BBvs->data[i];
    }

  for(int i=0;i<emxArray_idx->size[0]*emxArray_idx->size[1];++i){
    idx[i] = emxArray_idx->data[i];
  }

  for(int i=0;i<emxArray_ky->size[0];++i){
    for(int j=0;j<emxArray_ky->size[1];++j){
      ky[j+i*emxArray_ky->size[1]] = emxArray_ky->data[j+i*emxArray_ky->size[1]];
    }
  }

  for(int i=0;i<emxArray_C->size[0]*emxArray_C->size[1];++i){
    C[i] = emxArray_C->data[i];
  }

  IntegerVector dimC(2);
  dimC[0] = emxArray_C->size[0];
  dimC[1] = emxArray_C->size[1];
  C.attr("dim") = dimC;

  IntegerVector dimBB(3);
  IntegerVector dimBBvs(3);
  for(int i=0;i<3;++i){
    dimBB[i] = emxArray_BB->size[i];
    dimBBvs[i] = emxArray_BBvs->size[i];
  }
  BB.attr("dim") = dimBB;
  BBvs.attr("dim") = dimBBvs;

  IntegerVector dimidx(2);
  dimidx[0] = 1;
  dimidx[1] = emxArray_idx->size[0]*emxArray_idx->size[1];
  idx.attr("dim") = dimidx;

  IntegerVector dimky(2);
  dimky[0] = emxArray_ky->size[0];
  dimky[1] = emxArray_ky->size[1];
  ky.attr("dim") = dimky;

  IntegerVector dimnx(2);
  dimnx[0] = emxArray_x->size[0];
  dimnx[1] = emxArray_x->size[1];
  nx.attr("dim")=dimnx;
  Rcpp::List result = Rcpp::List::create(Rcpp::Named("dir")=BB,
                                         Rcpp::Named("ky")=ky,
                                         Rcpp::Named("idx") = idx);
                                         //Rcpp::Named("x")=nx,
                                         //Rcpp::Named("dirvs")=BBvs);
  return result;

}
