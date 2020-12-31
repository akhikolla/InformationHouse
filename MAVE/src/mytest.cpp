#include "rt_nonfinite.h"
#include "MAVEfast.h"
#include "MAVEfast_emxutil.h"
#include "power.h"
#include "eye.h"
#include "mldivide.h"
#include "repmat.h"
#include "exp.h"
#include "sort1.h"
#include "diag.h"
#include "eig.h"
#include "kron.h"
#include "sum.h"
#include "strcmp.h"
#include "mean.h"
#include "std.h"
#include "norm.h"
#include "floor.h"
#include "abs.h"
#include "unique.h"
#include "sortIdx.h"
#include "quantile.h"
#include "rdivide.h"
#include "sqrt.h"
#include "upper.h"
#include "bsxfun.h"

#include "MAVEfast.h"
#include "eig.h"
#include "MAVEfast_emxutil.h"
#include <vector>
#include <map>
#include <stdio.h>
#include <algorithm>
#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace std;



List mytestCpp(NumericMatrix C) {

  emxArray_real_T *emxArray_C;
  emxArray_creal_T *emxArray_V, *emxArray_D;
  emxInit_real_T(&emxArray_C,2);
  emxInit_creal_T(&emxArray_V,2);
  emxInit_creal_T(&emxArray_D,2);

  emxArray_C->size[0] = C.nrow(); emxArray_C->size[1] = C.ncol();
  emxEnsureCapacity((emxArray__common *)emxArray_C, 0, sizeof(double));
  for(int i=0;i<C.size();++i){
    emxArray_C->data[i] = C[i];
  }

  eig(emxArray_C,emxArray_V,emxArray_D);

  NumericMatrix D(emxArray_D->size[0],emxArray_D->size[1]);
  NumericMatrix V(emxArray_V->size[0],emxArray_V->size[1]);

  for(int i=0;i<C.size();++i){
    C[i] = emxArray_C->data[i];
    D[i] = abs(emxArray_D->data[i].re)+abs(emxArray_D->data[i].im);
    V[i] = abs(emxArray_V->data[i].re)+abs(emxArray_V->data[i].im);
  }

  Rcpp::List result = Rcpp::List::create(Rcpp::Named("C")=C,
                                         Rcpp::Named("D")=D,
                                         Rcpp::Named("V")=V);
  return result;

}
