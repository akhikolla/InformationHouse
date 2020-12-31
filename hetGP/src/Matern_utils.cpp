#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix matern5_2_1args(NumericMatrix X1){
  int nr = X1.nrow();
  int nc = X1.ncol();
  
  NumericMatrix s(nr, nr);
  NumericMatrix r(nr, nr);
  s.fill(1.);
  double tmp;
  
  // First compute polynomial term and distance
  const double* ptrX1 = (const double*) &X1(1,0);
  const double* ptrX2 = (const double*) &X1(0,0);
  double* ptrs = &s(0,1);
  double* ptrr = &r(0,1);
  
  for(int i = 1; i < nr; i++, ptrX1++){
    for(int j = 0; j < i; j++, ptrs++, ptrr++){
      for(int k = 0; k < nc; k++){
        tmp = sqrt(5.) * std::abs(*ptrX1 - *ptrX2);
        *ptrs *= (1 + tmp + tmp * tmp /3.);
        *ptrr -= tmp;
        ptrX1 += nr;
        ptrX2 += nr;
      }
      ptrX1 -= nr*nc;
      ptrX2 -= nr*nc - 1;
    }
    
    ptrX2 -= i;
    ptrs += (nr - i);
    ptrr += (nr - i);
  }
  
  // Now multiply by exponential part
  double* ptrs2 = &s(1,0); //symmetric
  ptrs = &s(0,1);
  ptrr = &r(0,1);
  for(int i = 1; i < nr; i++){
    for(int j = 0; j < i; j++, ptrs++, ptrr++){
      *ptrs *= exp(*ptrr);
      *ptrs2 = *ptrs;
      ptrs2 += nr;
    }
    ptrs += (nr - i);
    ptrr += (nr - i);
    ptrs2 += 1 - i*nr;
  }
  
  return s;
}


// [[Rcpp::export]]
NumericMatrix d_matern5_2_1args_theta_k_iso(NumericMatrix X1, double theta){
  int nr = X1.nrow();
  int nc = X1.ncol();
  NumericMatrix s(nr, nr);
  double tmp;
  
  const double* ptrX1 = (const double*) &X1(1,0);
  const double* ptrX2 = (const double*) &X1(0,0);
  double* ptrs = &s(0,1);
  double* ptrs2 = &s(1,0); //symmetric
  
  for(int i = 1; i < nr; i++, ptrX1++){
    for(int j = 0; j < i; j++, ptrs++){
      for(int k = 0; k < nc; k++){
        tmp = std::abs(*ptrX1 - *ptrX2) / theta;
        *ptrs -= ((10./3. - 5.) * tmp - 5 * sqrt(5.)/3. * tmp * tmp) / (1 + sqrt(5.) * tmp + 5./3. * tmp * tmp) * tmp/theta;
        ptrX1 += nr;
        ptrX2 += nr;
      }
      *ptrs2 = *ptrs;
      ptrs2 += nr;
      ptrX1 -= nr*nc;
      ptrX2 -= nr*nc - 1;
    }
    
    ptrX2 -= i;
    ptrs += (nr - i);
    ptrs2 += 1 - i*nr;
  }
  return s;
}

// [[Rcpp::export]]
NumericMatrix d_matern5_2_1args_theta_k(NumericMatrix X1, double theta){
  // X1 has just one column here
  int nr = X1.nrow();
  NumericMatrix s(nr, nr);
  double tmp;
  
  const double* ptrX1 = (const double*) &X1(1,0);
  const double* ptrX2 = (const double*) &X1(0,0);
  double* ptrs = &s(0,1);
  double* ptrs2 = &s(1,0); //symmetric
  
  for(int i = 1; i < nr; i++, ptrX1++){
    for(int j = 0; j < i; j++, ptrs++){
      tmp = std::abs(*ptrX1 - *ptrX2) / theta;
      *ptrs -= ((10./3. - 5.) * tmp - 5 * sqrt(5.)/3. * tmp * tmp) / (1 + sqrt(5.) * tmp + 5./3. * tmp * tmp) * tmp/theta;
      
      *ptrs2 = *ptrs;
      ptrs2 += nr;
      ptrX2 ++;
    }
    
    ptrX2 -= i;
    ptrs += (nr - i);
    ptrs2 += 1 - i*nr;
  }
  return s;
}

// [[Rcpp::export]]
NumericMatrix d_matern5_2_1args_kthetag(NumericMatrix X1, double kt){
  int nr = X1.nrow();
  int nc = X1.ncol();
  NumericMatrix s(nr, nr);
  double tmp;
  
  const double* ptrX1 = (const double*) &X1(1,0);
  const double* ptrX2 = (const double*) &X1(0,0);
  double* ptrs = &s(0,1);
  double* ptrs2 = &s(1,0); //symmetric
  
  for(int i = 1; i < nr; i++, ptrX1++){
    for(int j = 0; j < i; j++, ptrs++){
      for(int k = 0; k < nc; k++){
        tmp = std::abs(*ptrX1 - *ptrX2) / kt;
        *ptrs -= ((10./3. - 5.) * tmp - 5 * sqrt(5.)/3. * tmp * tmp) / (1 + sqrt(5.) * tmp + 5./3. * tmp * tmp) * tmp/kt;
        ptrX1 += nr;
        ptrX2 += nr;
      }
      *ptrs2 = *ptrs;
      ptrs2 += nr;
      ptrX1 -= nr*nc;
      ptrX2 -= nr*nc - 1;
    }
    
    ptrX2 -= i;
    ptrs += (nr - i);
    ptrs2 += 1 - i*nr;
  }
  return s;
}

// [[Rcpp::export]]
NumericMatrix matern5_2_2args(NumericMatrix X1, NumericMatrix X2){
  int nr1 = X1.nrow();
  int nr2 = X2.nrow();
  int dim = X1.ncol();
  
  NumericMatrix s(nr1, nr2);
  s.fill(1.);
  NumericMatrix r(nr1, nr2);
  double tmp;
  
  double* ptrs = &s(0,0);
  double* ptrr = &r(0,0);
  const double* ptrX2 = (const double*) &X2(0,0);
  const double* ptrX1 = (const double*) &X1(0,0);
  
  // Polynomial part
  for(int i = 0; i < nr2; i++){
    for(int j = 0; j < nr1; j++, ptrs++, ptrr++){
      for(int k = 0; k < dim; k++){
        tmp = sqrt(5.) * std::abs(*ptrX1 - *ptrX2);
        *ptrs *= (1 + tmp + tmp * tmp / 3.);
        *ptrr -= tmp;
        ptrX1 += nr1;
        ptrX2 += nr2;
      }
      ptrX2 -= nr2*dim;
      ptrX1 -= nr1*dim - 1;
      
    }
    ptrX2++;
    ptrX1 -= nr1;
  }
  
  ptrs = &s(0,0);
  ptrr = &r(0,0);
  
  // Exponential part
  for(int i = 0; i < nr2; i++){
    for(int j = 0; j < nr1; j++, ptrs++, ptrr++){
      *ptrs *= exp(*ptrr);
    }
  }
  
  return s;
}


// [[Rcpp::export]]
NumericMatrix d_matern5_2_2args_theta_k_iso(NumericMatrix X1, NumericMatrix X2, double theta){
  int nr1 = X1.nrow();
  int nr2 = X2.nrow();
  int dim = X1.ncol();
  NumericMatrix s(nr1, nr2);
  double tmp;
  
  double* ptrs = &s(0,0);
  const double* ptrX2 = (const double*) &X2(0,0);
  const double* ptrX1 = (const double*) &X1(0,0);
  
  for(int i = 0; i < nr2; i++){
    for(int j = 0; j < nr1; j++, ptrs++){
      for(int k = 0; k < dim; k++){
        tmp = std::abs(*ptrX1 - *ptrX2) / theta;
        *ptrs -= ((10./3. - 5.) * tmp - 5 * sqrt(5.)/3. * tmp * tmp) / (1 + sqrt(5.) * tmp + 5./3. * tmp * tmp) * tmp / theta;
        ptrX1 += nr1;
        ptrX2 += nr2;
      }
      ptrX2 -= nr2*dim;
      ptrX1 -= nr1*dim - 1;
    }
    ptrX2++;
    ptrX1 -= nr1;
  }
  return s;
}


// [[Rcpp::export]]
NumericMatrix d_matern5_2_2args_kthetag(NumericMatrix X1, NumericMatrix X2, double kt){
  int nr1 = X1.nrow();
  int nr2 = X2.nrow();
  int dim = X1.ncol();
  NumericMatrix s(nr1, nr2);
  double tmp;
  
  double* ptrs = &s(0,0);
  const double* ptrX2 = (const double*) &X2(0,0);
  const double* ptrX1 = (const double*) &X1(0,0);
  for(int i = 0; i < nr2; i++){
    for(int j = 0; j < nr1; j++, ptrs++){
      for(int k = 0; k < dim; k++){
        tmp = std::abs(*ptrX1 - *ptrX2) / kt;
        *ptrs -= ((10./3. - 5.) * tmp - 5 * sqrt(5.)/3. * tmp * tmp) / (1 + sqrt(5.) * tmp + 5./3. * tmp * tmp) * tmp/kt;
        ptrX1 += nr1;
        ptrX2 += nr2;
      }
      ptrX2 -= nr2*dim;
      ptrX1 -= nr1*dim - 1;
      
    }
    ptrX2++;
    ptrX1 -= nr1;
  }
  return s;
}

// [[Rcpp::export]]
NumericMatrix partial_d_dist_abs_dX_i1_i2(NumericMatrix X1, int i1, int i2){
  int nr = X1.nrow();
  NumericMatrix s(nr, nr);
  double tmp;
  
  for(int i = 0; i < nr; i++){
    if(i == (i1 - 1))
      continue;
    tmp = (X1(i1 - 1, i2 - 1) - X1(i, i2 - 1)) ;
    if(tmp > 0){
      s(i1 - 1, i) = s(i, i1 - 1) = ((10./3. - 5.) * tmp - 5 * sqrt(5.)/3. * tmp * tmp) / (1. + sqrt(5.) * tmp + 5./3. * tmp * tmp);
    }else{
      if(tmp == 0){
        s(i1 - 1, i) = s(i, i1 - 1) = 0;
      }else{
        tmp = std::abs(tmp);
        s(i1 - 1, i) = s(i, i1 - 1) = -((10./3. - 5.) * tmp - 5. * sqrt(5.)/3. * tmp * tmp) / ((1. + sqrt(5.) * tmp + 5./3. * tmp * tmp));
      }
    }
  }
  return s;
}

// [[Rcpp::export]]
NumericMatrix partial_d_dist_abs_dX1_i1_i2_X2(NumericMatrix X1, NumericMatrix X2, int i1, int i2){
  int nr = X2.nrow();
  NumericMatrix s(X1.nrow(), nr);
  double tmp;
  
  for(int i = 0; i < nr; i++){
    tmp = X1(i1-1, i2-1) - X2(i, i2-1);
    if(tmp > 0){
      s(i1 - 1, i) = ((10./3. - 5.) * tmp - 5 * sqrt(5.)/3. * tmp * tmp) / (1. + sqrt(5.) * tmp + 5./3. * tmp * tmp);
    }else{
      if(tmp == 0){
        s(i1 - 1, i) = 0;
      }else{
        tmp = std::abs(tmp);
        s(i1 - 1, i) = -((10./3. - 5.) * tmp - 5. * sqrt(5.)/3. * tmp * tmp) / ((1. + sqrt(5.) * tmp + 5./3. * tmp * tmp));
      }
    }
  }
  return s;
}



// [[Rcpp::export]]
NumericMatrix matern3_2_1args(NumericMatrix X1){
  int nr = X1.nrow();
  int nc = X1.ncol();
  
  NumericMatrix s(nr, nr);
  NumericMatrix r(nr, nr);
  s.fill(1.);
  double tmp;
  
  // First compute polynomial term and distance
  const double* ptrX1 = (const double*) &X1(1,0);
  const double* ptrX2 = (const double*) &X1(0,0);
  double* ptrs = &s(0,1);
  double* ptrr = &r(0,1);
  
  for(int i = 1; i < nr; i++, ptrX1++){
    for(int j = 0; j < i; j++, ptrs++, ptrr++){
      for(int k = 0; k < nc; k++){
        tmp = sqrt(3.) * std::abs(*ptrX1 - *ptrX2);
        *ptrs *= (1 + tmp);
        *ptrr -= tmp;
        ptrX1 += nr;
        ptrX2 += nr;
      }
      ptrX1 -= nr*nc;
      ptrX2 -= nr*nc - 1;
    }
    
    ptrX2 -= i;
    ptrs += (nr - i);
    ptrr += (nr - i);
  }
  
  // Now multiply by exponential part
  double* ptrs2 = &s(1,0); //symmetric
  ptrs = &s(0,1);
  ptrr = &r(0,1);
  for(int i = 1; i < nr; i++){
    for(int j = 0; j < i; j++, ptrs++, ptrr++){
      *ptrs *= exp(*ptrr);
      *ptrs2 = *ptrs;
      ptrs2 += nr;
    }
    ptrs += (nr - i);
    ptrr += (nr - i);
    ptrs2 += 1 - i*nr;
  }
  
  return s;
}


// [[Rcpp::export]]
NumericMatrix d_matern3_2_1args_theta_k_iso(NumericMatrix X1, double theta){
  int nr = X1.nrow();
  int nc = X1.ncol();
  NumericMatrix s(nr, nr);
  double tmp;
  
  const double* ptrX1 = (const double*) &X1(1,0);
  const double* ptrX2 = (const double*) &X1(0,0);
  double* ptrs = &s(0,1);
  double* ptrs2 = &s(1,0); //symmetric
  
  for(int i = 1; i < nr; i++, ptrX1++){
    for(int j = 0; j < i; j++, ptrs++){
      for(int k = 0; k < nc; k++){
        tmp = std::abs(*ptrX1 - *ptrX2) / theta;
        *ptrs -= 3*tmp / (1 + sqrt(3.) * tmp) * tmp / theta;
        ptrX1 += nr;
        ptrX2 += nr;
      }
      *ptrs2 = *ptrs;
      ptrs2 += nr;
      ptrX1 -= nr*nc;
      ptrX2 -= nr*nc - 1;
    }
    
    ptrX2 -= i;
    ptrs += (nr - i);
    ptrs2 += 1 - i*nr;
  }
  return s;
}

// [[Rcpp::export]]
NumericMatrix d_matern3_2_1args_theta_k(NumericMatrix X1, double theta){
  // X1 has one column
  int nr = X1.nrow();
  NumericMatrix s(nr, nr);
  double tmp;
  
  const double* ptrX1 = (const double*) &X1(1,0);
  const double* ptrX2 = (const double*) &X1(0,0);
  double* ptrs = &s(0,1);
  double* ptrs2 = &s(1,0); //symmetric
  
  for(int i = 1; i < nr; i++, ptrX1++){
    for(int j = 0; j < i; j++, ptrs++){
      tmp = std::abs(*ptrX1 - *ptrX2) / theta;
      *ptrs -= 3*tmp / (1 + sqrt(3.) * tmp) * tmp / theta;
      *ptrs2 = *ptrs;
      ptrs2 += nr;
      ptrX2++;
    }
    
    ptrX2 -= i;
    ptrs += (nr - i);
    ptrs2 += 1 - i*nr;
  }
  return s;
}

// [[Rcpp::export]]
NumericMatrix d_matern3_2_1args_kthetag(NumericMatrix X1, double kt){
  int nr = X1.nrow();
  int nc = X1.ncol();
  NumericMatrix s(nr, nr);
  double tmp;
  
  const double* ptrX1 = (const double*) &X1(1,0);
  const double* ptrX2 = (const double*) &X1(0,0);
  double* ptrs = &s(0,1);
  double* ptrs2 = &s(1,0); //symmetric
  
  for(int i = 1; i < nr; i++, ptrX1++){
    for(int j = 0; j < i; j++, ptrs++){
      for(int k = 0; k < nc; k++){
        tmp = std::abs(*ptrX1 - *ptrX2) / kt;
        *ptrs -= 3*tmp / (1 + sqrt(3.) * tmp) * tmp/kt;
        ptrX1 += nr;
        ptrX2 += nr;
      }
      *ptrs2 = *ptrs;
      ptrs2 += nr;
      ptrX1 -= nr*nc;
      ptrX2 -= nr*nc - 1;
    }
    
    ptrX2 -= i;
    ptrs += (nr - i);
    ptrs2 += 1 - i*nr;
  }
  return s;
}

// [[Rcpp::export]]
NumericMatrix matern3_2_2args(NumericMatrix X1, NumericMatrix X2){
  int nr1 = X1.nrow();
  int nr2 = X2.nrow();
  int dim = X1.ncol();
  
  NumericMatrix s(nr1, nr2);
  s.fill(1.);
  NumericMatrix r(nr1, nr2);
  double tmp;
  
  double* ptrs = &s(0,0);
  double* ptrr = &r(0,0);
  const double* ptrX2 = (const double*) &X2(0,0);
  const double* ptrX1 = (const double*) &X1(0,0);
  
  // Polynomial part
  for(int i = 0; i < nr2; i++){
    for(int j = 0; j < nr1; j++, ptrs++, ptrr++){
      for(int k = 0; k < dim; k++){
        tmp = sqrt(3.) * std::abs(*ptrX1 - *ptrX2);
        *ptrs *= (1 + tmp);
        *ptrr -= tmp;
        ptrX1 += nr1;
        ptrX2 += nr2;
      }
      ptrX2 -= nr2*dim;
      ptrX1 -= nr1*dim - 1;
      
    }
    ptrX2++;
    ptrX1 -= nr1;
  }
  
  ptrs = &s(0,0);
  ptrr = &r(0,0);
  
  // Exponential part
  for(int i = 0; i < nr2; i++){
    for(int j = 0; j < nr1; j++, ptrs++, ptrr++){
      *ptrs *= exp(*ptrr);
    }
  }
  
  return s;
}


// [[Rcpp::export]]
NumericMatrix d_matern3_2_2args_theta_k_iso(NumericMatrix X1, NumericMatrix X2, double theta){
  int nr1 = X1.nrow();
  int nr2 = X2.nrow();
  int dim = X1.ncol();
  NumericMatrix s(nr1, nr2);
  double tmp;
  
  double* ptrs = &s(0,0);
  const double* ptrX2 = (const double*) &X2(0,0);
  const double* ptrX1 = (const double*) &X1(0,0);
  
  for(int i = 0; i < nr2; i++){
    for(int j = 0; j < nr1; j++, ptrs++){
      for(int k = 0; k < dim; k++){
        tmp = std::abs(*ptrX1 - *ptrX2) / theta;
        *ptrs -= 3*tmp / (1 + sqrt(3.) * tmp) * tmp / theta;
        ptrX1 += nr1;
        ptrX2 += nr2;
      }
      ptrX2 -= nr2*dim;
      ptrX1 -= nr1*dim - 1;
    }
    ptrX2++;
    ptrX1 -= nr1;
  }
  return s;
}


// [[Rcpp::export]]
NumericMatrix d_matern3_2_2args_kthetag(NumericMatrix X1, NumericMatrix X2, double kt){
  int nr1 = X1.nrow();
  int nr2 = X2.nrow();
  int dim = X1.ncol();
  NumericMatrix s(nr1, nr2);
  double tmp;
  
  double* ptrs = &s(0,0);
  const double* ptrX2 = (const double*) &X2(0,0);
  const double* ptrX1 = (const double*) &X1(0,0);
  for(int i = 0; i < nr2; i++){
    for(int j = 0; j < nr1; j++, ptrs++){
      for(int k = 0; k < dim; k++){
        tmp = std::abs(*ptrX1 - *ptrX2) / kt;
        *ptrs -= 3*tmp / (1 + sqrt(3.) * tmp) * tmp/kt;
        ptrX1 += nr1;
        ptrX2 += nr2;
      }
      ptrX2 -= nr2*dim;
      ptrX1 -= nr1*dim - 1;
      
    }
    ptrX2++;
    ptrX1 -= nr1;
  }
  return s;
}
