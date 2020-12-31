#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

arma::mat shiftedlegendre(int);
void outer_cpp(arma::mat& A, const arma::vec& x, const arma::vec& y);
void zero_lower_diagtri(arma::mat& A);

// [[Rcpp::export]]
arma::mat Lmoments_calc(const arma::mat& data, unsigned int rmax = 4) {
  const int n = data.n_rows;
  const int p = data.n_cols;
  const int _rmax = (int) std::min((double) rmax, (double) n);
  
  arma::mat x(n, p, arma::fill::zeros);
  arma::mat L(p, _rmax, arma::fill::zeros);
  for (int i = 0; i < p; i++) {
    x.col(i) = arma::sort(data.col(i));
  }
  
  if (_rmax == 1) {
    for(int i = 0; i < p; i++) {
      L(i, 0) = arma::mean(x.col(i));
    }
    return L;
  }
  
  arma::mat bcoef(n, _rmax, arma::fill::zeros);
  arma::cube bcoefm(n, p, _rmax, arma::fill::zeros);
  arma::mat b(p, _rmax, arma::fill::zeros);
  
  bcoef.col(0) = arma::linspace(0.0, 1.0, n);
  for (int i = 0; i < p; i++) {
    bcoefm.slice(0).col(i) = bcoef.col(0); 
    b(i, 0) = arma::mean(x.col(i)); 
    b(i, 1) = arma::mean(bcoefm.slice(0).col(i) % x.col(i)); 
  }
  L.col(0) = b.col(0);
  
  if (_rmax > 2) {
    for(int r = 2; r <= _rmax - 1; r++) {
      bcoef.col(r - 1) = bcoef.col(r - 2) % arma::linspace(-(r - 1.0) / (n - r), 1.0, n);
      for(int j = 0; j < p; j++) {
        bcoefm.slice(r - 1).col(j) = bcoef.col(r - 1);
        b(j, r) = arma::mean(bcoefm.slice(r - 1).col(j) % x.col(j));
      }
    }
  }
  
  for(int r = 1; r <= (_rmax - 1); r++) {
    L.col(r).zeros();
    for(int k = 0; k <= r; k++) {
      L.col(r) = L.col(r) + std::pow(-1.0, r - k) * tgamma(r + k + 1.0)/std::pow(tgamma(k + 1.0), 2.0)/tgamma(r - k + 1.0) * b.col(k);
    }
  }
  
  return L;
}

// [[Rcpp::export]]
Rcpp::List Lmomcov_calc(const arma::mat& data, unsigned int rmax = 4) {
  
  const int n = data.n_rows;
  const int p = data.n_cols;
  const int _rmax = (int) std::min((double) rmax, std::floor(n / 2.0));
  
  if (_rmax <= 1) {
    return Rcpp::List::create(NA_REAL);
  }
  
  arma::mat C = shiftedlegendre(_rmax);
  arma::mat x(n, p, arma::fill::zeros);
  for (int i = 0; i < p; i++) {
    x.col(i) = arma::sort(data.col(i));
  }
  
  arma::mat bcoef(n, _rmax, arma::fill::zeros);
  arma::mat b1coef(n, _rmax, arma::fill::zeros);
  arma::cube b2coef(n, _rmax, _rmax, arma::fill::zeros);
  arma::mat b20coef(n, _rmax, arma::fill::zeros);
  arma::cube bcoefm(n, p, _rmax, arma::fill::zeros);
  arma::mat b(p, _rmax, arma::fill::zeros);
  
  bcoef.col(0) = arma::linspace(0.0, 1.0, n);
  for (int i = 0; i < p; i++) {
    bcoefm.slice(0).col(i) = bcoef.col(0); 
    b(i, 0) = arma::mean(x.col(i)); 
    b(i, 1) = arma::mean(bcoefm.slice(0).col(i) % x.col(i)); 
  }
  
  b20coef.col(0) = arma::linspace(-1.0 / (n - 2), 1.0, n);
  b1coef.col(0) = arma::linspace(0.0, (n - 1.0) / (n - 2), n);
  
  if (_rmax > 2) {
    for(int r = 2; r <= _rmax - 1; r++) {
      bcoef.col(r - 1) = bcoef.col(r - 2) % arma::linspace(-(r - 1.0) / (n - r), 1.0, n);
      for(int j = 0; j < p; j++) {
       bcoefm.slice(r - 1).col(j) = bcoef.col(r - 1);
       b(j, r) = arma::mean(bcoefm.slice(r - 1).col(j) % x.col(j));
      }
    }
  }
  
  if (_rmax > 1) {
    for(int k = 1; k <= (_rmax - 1); k++) {
      if (k > 1) {
        b20coef.col(k - 1) = b20coef.col(k - 2) % arma::linspace(-k / (n - 1.0 - k), 1.0, n);
        b1coef.col(k - 1) = b1coef.col(k - 2) % arma::linspace((-k + 1.0) / (n - 1.0 - k), (n - k) / (n - 1.0 - k), n);
      }
      b2coef.slice(k - 1).col(0) = arma::linspace((-k - 1.0) / (n - k - 2.0), 1.0, n);
      if (_rmax > 2) {
        for(int l = 2; l <= (_rmax - 1); l++) {
          b2coef.slice(k - 1).col(l - 1) = b2coef.slice(k - 1).col(l - 2) % arma::linspace((-k - l) / (n - k - l - 1.0), 1.0, n); 
        }
      }
    }
  }
  
  Rcpp::List covmatrixlist = Rcpp::List(p);
  arma::mat theta(_rmax, _rmax, arma::fill::zeros);
  arma::mat xx(n, n, arma::fill::zeros);
  arma::mat term1(n, n, arma::fill::zeros);
  arma::mat term2(n, n, arma::fill::zeros);
  arma::vec one(n);
  one.ones();
  double jointbb = 0.0;
  
  for(int i = 0; i < p; i++) {
    
    outer_cpp(xx, x.col(i), x.col(i));
    zero_lower_diagtri(xx);
    
    for(int k = 0; k <= _rmax - 1; k++) {
      for(int l = 0; l <= _rmax - 1; l++) {
        if (k > 0 && l > 0) {
          outer_cpp(term1, b1coef.col(k - 1), b2coef.slice(k - 1).col(l - 1));
          outer_cpp(term2, b1coef.col(l - 1), b2coef.slice(l - 1).col(k - 1));
        }
        if (k == 0 && l > 0) {
          outer_cpp(term1, one, b20coef.col(l - 1));
          outer_cpp(term2, b1coef.col(l - 1), one);
        }
        if (k > 0 && l == 0) {
          outer_cpp(term1, b1coef.col(k - 1), one);
          outer_cpp(term2, one, b20coef.col(k - 1));
        }
        if (k == 0 && l == 0) {
          outer_cpp(term1, one, one);
          outer_cpp(term2, one, one);
        }
        zero_lower_diagtri(term1);
        zero_lower_diagtri(term2);
        jointbb = accu((term1 + term2) % xx) / (n * (n - 1.0));
        theta(k, l) = b(i, k) * b(i, l) - jointbb;
      }
      
    }
    covmatrixlist(i) = C * theta * C.t();
  }
  return covmatrixlist;
  // return Rcpp::List::create(Rcpp::Named("res") = covmatrixlist,
  //                           Rcpp::Named("bcoef") = bcoef,
  //                           Rcpp::Named("b1coef") = b1coef,
  //                           Rcpp::Named("b2coef") = b2coef,
  //                           Rcpp::Named("b20coef") = b20coef,
  //                           Rcpp::Named("bcoefm") = bcoefm,
  //                           Rcpp::Named("b") = b,
  //                           Rcpp::Named("theta") = theta,
  //                           Rcpp::Named("xx") = xx,
  //                           Rcpp::Named("term1") = term1,
  //                           Rcpp::Named("term2") = term2,
  //                           Rcpp::Named("C") = C);
}

void outer_cpp(arma::mat& A, const arma::vec& x, const arma::vec& y) {
  for(unsigned int j = 0; j < A.n_cols; j++) {
    for (unsigned int i = 0; i < A.n_rows; i++) {
      A(i, j) = y(j) * x(i);
    }
  }
}

void zero_lower_diagtri(arma::mat& A) {
  for(unsigned int j = 0; j < A.n_cols; j++) {
    for (unsigned int i = 0; i < A.n_rows; i++) {
      if (j <= i) {
        A(i, j) = 0.0;
      }
    }
  }
}

// [[Rcpp::export]]
arma::mat shiftedlegendre(int rmax) {
  if (rmax <= 0) {
    throw std::range_error("'rmax' must be > 0");
  }
  arma::mat C(rmax, rmax, arma::fill::zeros);
  C(0, 0) = 1.0;
  
  if (rmax > 1) {
    C(0, 1) = -1.0;
    C(1, 1) = 2.0;
    if (rmax > 2) {
      int kn;
      for(int k = 2; k < rmax; k++) {
        kn = (k + 1) - 2;
        arma::vec sameorder = (-(2.0 * kn + 1.0) * C(arma::span(0, k - 1), k - 1) - kn * C(arma::span(0, k - 1), k - 2)) / (kn + 1.0);
        arma::vec uporder = 2.0 * (2.0 * kn + 1.0) / (kn + 1.0) * C(arma::span(0, k - 1), k - 1);
        C(arma::span(0, k - 1), k) = sameorder;
        C(arma::span(1, k), k) = uporder + C(arma::span(1, k), k);
      }
    }
  }
  return C.t();
}

// [[Rcpp::export]]
arma::mat t1lmoments_calc(const arma::mat& data, unsigned int rmax = 4) {
  
  const int rmax_out = (int) std::min((double) rmax, (double) 4);
  
  const int n = data.n_rows;
  const int p = data.n_cols;
  
  arma::mat x(n, p, arma::fill::zeros);
  arma::mat L(p, rmax);
  
  for (int j = 0; j < p; j++) {
    x.col(j) = arma::sort(data.col(j));
  }
  arma::vec i = arma::linspace(3, n, n - 2);
  
  arma::vec s11(n - 1, arma::fill::ones);
  arma::vec s12(n - 1, arma::fill::ones);
  arma::vec s13(n - 1, arma::fill::ones);
  arma::vec s14(n - 1, arma::fill::ones);
  arma::vec s21(n - 1, arma::fill::ones);
  arma::vec s22(n - 1, arma::fill::ones);
  arma::vec s23(n - 1, arma::fill::ones);
  arma::vec s31(n - 1, arma::fill::ones);
  arma::vec s32(n - 1, arma::fill::ones);
  arma::vec s41(n - 1, arma::fill::ones);
  
  //std::cout << "First checkpoint ok" << std::endl;
  
  s11(arma::span(1, n - 2)) = (i - 1) / (i - 2) % (n - i) / (n - i + 1);
  s12(arma::span(1, n - 2)) = (i - 1) / (i - 2) % (n - i - 1) / (n - i + 1);
  s13(arma::span(1, n - 2)) = (i - 1) / (i - 2) % (n - i - 2) / (n - i + 1);
  s14(arma::span(1, n - 2)) = (i - 1) / (i - 2) % (n - i - 3) / (n - i + 1);
  
  s21(arma::span(1, n - 2)) = (i - 1) / (i - 3) % (n - i) / (n - i + 1);
  s22(arma::span(1, n - 2)) = (i - 1) / (i - 3) % (n - i - 1) / (n - i + 1);
  s23(arma::span(1, n - 2)) = (i - 1) / (i - 3) % (n - i - 2) / (n - i + 1);
  
  s31(arma::span(1, n - 2)) = (i - 1) / (i - 4) % (n - i) / (n - i + 1);
  s32(arma::span(1, n - 2)) = (i - 1) / (i - 4) % (n - i - 1) / (n - i + 1);
  
  s41(arma::span(1, n - 2)) = (i - 1) / (i - 5) % (n - i) / (n - i + 1);
  
  //std::cout << "Second checkpoint ok" << std::endl;
    
  s21(1) = 1.0;
  s22(1) = 1.0;
  s23(1) = 1.0;
  s31(1) = 1.0;
  s31(2) = 1.0;
  s32(1) = 1.0;
  s32(2) = 1.0;
  s41(1) = 1.0;
  s41(2) = 1.0;
  s41(3) = 1.0;
  
  arma::vec c11 = Rf_choose(n - 2, 1) / Rf_choose(n, 3) * arma::cumprod(s11);
  arma::vec c12 = Rf_choose(n - 2, 2) / Rf_choose(n, 4) * arma::cumprod(s12);
  arma::vec c13 = Rf_choose(n - 2, 3) / Rf_choose(n, 5) * arma::cumprod(s13);
  arma::vec c14 = Rf_choose(n - 2, 4) / Rf_choose(n, 6) * arma::cumprod(s14);
  arma::vec c21 = Rf_choose(n - 3, 1) / Rf_choose(n, 4) * arma::cumprod(s21);
  arma::vec c22 = Rf_choose(n - 3, 2) / Rf_choose(n, 5) * arma::cumprod(s22);
  arma::vec c23 = Rf_choose(n - 3, 3) / Rf_choose(n, 6) * cumprod(s23);
  arma::vec c31 = Rf_choose(n - 4, 1) / Rf_choose(n, 5) * cumprod(s31);
  arma::vec c32 = Rf_choose(n - 4, 2) / Rf_choose(n, 6) * cumprod(s32);
  arma::vec c41 = Rf_choose(n - 5, 1) / Rf_choose(n, 6) * cumprod(s41);
  
  c21(0) = 0.0;
  c22(0) = 0.0;
  c23(0) = 0.0;
  c31(0) = 0.0;
  c31(1) = 0.0;
  c32(0) = 0.0;
  c32(1) = 0.0;
  c41(0) = 0.0;
  c41(1) = 0.0;
  c41(2) = 0.0;
  
  for (int j = 0; j < p; j++) {
    L(j, 0) = arma::accu(c11(arma::span(0, n - 3)) % x(arma::span(1, n - 2), j));
    L(j, 1) = arma::accu((c21(arma::span(0, n - 3)) - c12(arma::span(0, n - 3))) % x(arma::span(1, n - 2), j)) / 2.0;
    L(j, 2) = arma::accu((c31(arma::span(0, n - 3)) - 2 * c22(arma::span(0, n - 3)) + c13(arma::span(0, n - 3))) % 
      x(arma::span(1, n - 2), j)) / 3.0;
    L(j, 3) = arma::accu( (c41(arma::span(0, n - 3)) - 3 * c32(arma::span(0, n - 3)) + 
      3 * c23(arma::span(0, n - 3)) - c14(arma::span(0, n - 3))) % x(arma::span(1, n - 2), j)) / 4.0;
  }
  // 
  // Rcpp::List s = Rcpp::List::create(Rcpp::Named("res") = L.cols(0, rmax_out - 1),
  //                                   Rcpp::Named("s11") = s11,
  //                                   Rcpp::Named("s12") = s12,
  //                                   Rcpp::Named("s13") = s13,
  //                                   Rcpp::Named("s14") = s14,
  //                                   Rcpp::Named("s21") = s21,
  //                                   Rcpp::Named("s22") = s22,
  //                                   Rcpp::Named("s23") = s23,
  //                                   Rcpp::Named("s31") = s31,
  //                                   Rcpp::Named("s32") = s32,
  //                                   Rcpp::Named("s41") = s41);
  // Rcpp::List c = Rcpp::List::create(Rcpp::Named("c11") = c11,
  //                                   Rcpp::Named("c12") = c12,
  //                                   Rcpp::Named("c13") = c13,
  //                                   Rcpp::Named("c14") = c14,
  //                                   Rcpp::Named("c21") = c21,
  //                                   Rcpp::Named("c22") = c22,
  //                                   Rcpp::Named("c23") = c23,
  //                                   Rcpp::Named("c31") = c31,
  //                                   Rcpp::Named("c32") = c32,
  //                                   Rcpp::Named("c41") = c41);
  //return Rcpp::List::create(Rcpp::Named("s") = s,
  //                          Rcpp::Named("c") = c);
  return L.cols(0, rmax_out - 1);
}













