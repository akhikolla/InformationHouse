/******************  created by Erlis Ruli 14/01/2014  ************************
****************** Robust Approximate Bayesian Infernece  *********************
******************************************************************************/
#include "extrafuns.h"
#include <RcppArmadillo.h>
using namespace Rcpp;


// [[Rcpp::export]]
arma::vec
funPsi_rls(double mu,
          double sig,
          arma::vec y,
          int n,
          const double c1,
          const double c2,
          const double kc2){
  arma::vec ans(2, arma::fill::zeros);
  double ystd = 0.0;
  for(int i=0; i<n; i++){
    ystd = (y(i) - mu)/sig;
    ans(0) += psi_huber(ystd, c1);
    ans(1) += psi_huber(ystd, c2)*psi_huber(ystd, c2);
  }
  ans(1) = ans(1) - n*kc2;
  return ans;
}

// [[Rcpp::export]]
arma::mat Owen_Psi(double mu,
        double sig,
        arma::vec y,
        int n,
        const double c1,
        const double c2,
        const double kc2){

  arma::mat ans(2, n, arma::fill::zeros);
  double ystd = 0.0;
  for(int i=0; i<n; i++){
    ystd = (y(i) - mu)/sig;
    ans(0,i) = psi_huber(ystd, c1);
    ans(1,i) = psi_huber(ystd, c2)*psi_huber(ystd, c2);
  }
  ans.row(1) = ans.row(1) - kc2;
  return ans;
}

// [[Rcpp::export]]
double
Owen_Lmult(arma::vec lambda,
          double mu,
          double sig,
          arma::vec y,
          int n,
          const double c1,
          const double c2,
          const double kc2){

  arma::mat gr = Owen_Psi(mu, sig, y, n, c1, c2, kc2);
  double eps = 1.0/n, ans = 0.0;
  arma::vec z = 1 + trans(trans(lambda)*gr);

  for(int i=0; i<n; i++){
    if (z(i) < eps){
      ans += log(eps) - 1.5 + 2*z(i)/eps - 0.5*(z(i)/eps)*(z(i)/eps);
    } else {
      ans += log(z(i));
    }
  }
  return -ans;
}


// [[Rcpp::export]]
double
dPrior(double mu, double lsig, double muSD, double sigScale, bool lg = false){

  double dens = Rf_dnorm4(mu, 0.0, muSD, true) +
                dhalfCauchy(exp(lsig), sigScale, true) +
                lsig;
  if (lg == false)
      dens = exp(dens);

  return dens;
}

//[[Rcpp::export]]
arma::vec simData_rls(double mu, double sig, int n){
  arma::vec ans(n);
  RNGScope scope;
  for(int i=0; i<n; i++){
    ans(i) = Rf_rnorm(mu, sig);
  }
  return ans;
}

//[[Rcpp::export]]
arma::mat ABCrls(int nabc,
            double eps,
            arma::vec y,
            const double c1,
            const double c2,
            const double kc2,
            arma::vec hatPar,
            arma::mat hatVcov,
            arma::vec hatGrad,
            arma::mat hatJinv,
            double muSD,
            double sigScale,
            arma::vec thin)
{
  int nacc = 1,
  nPar = 2,
  n = y.n_elem;

  arma::mat ans(nPar, nabc, arma::fill::zeros);
  arma::vec ysim(n, arma::fill::zeros),
  propVal(nPar, arma::fill::zeros),
  simGrad(nPar, arma::fill::zeros);

  double df = 5.0, sdist = 0.0,
  rUnif = 0.0, logAccRat = 0.0;

  RNGScope scope;

  /*start algorithm from hat.theta*/
  ans.col(0) = hatPar;

    /* for i to nabc do*/
  for(int i = 1; i < nabc; i++){
    // draw from the jumping density
    propVal = rmvt(ans.col(i-1), hatVcov, nPar, df);

    // acceptance ratio of prior and jumping density
    logAccRat = dPrior(propVal(0), propVal(1), muSD, sigScale, true) -
                //dmvt(ans.col(i-1), propVal, hatVcov, nPar, df, true) -
                dPrior(ans(0, i-1), ans(1, i-1), muSD, sigScale, true);
                //dmvt(propVal, ans.col(i-1), hatVcov, nPar, df, true));

    rUnif = as_scalar(arma::randu(1));

    if (log(rUnif) > logAccRat) {

      ans.col(i) = ans.col(i-1);

    } else {
      // compute the ABC disntace
      ysim = simData_rls(propVal(0), exp(propVal(1)), n);
      simGrad = funPsi_rls(hatPar(0), exp(hatPar(1)), ysim, n, c1, c2, kc2);
      sdist = arma::as_scalar(trans(simGrad-hatGrad)*hatJinv*(simGrad-hatGrad));

      if (sdist > eps){
              ans.col(i) = ans.col(i-1);
              } else {
    // else accept
      nacc += 1;
      ans.col(i) = propVal;
      //veps(i) = sdist/n_ind;
      }
    }
    // print some information on the simulation process
    if(i%5000 == 0 ) Rprintf("\rABC %d of %d, accpted %d, acc.ratio %f\r", i,
    nabc, nacc, static_cast<double>(nacc)/static_cast<double>(i));

    // for stopping R computations
    R_CheckUserInterrupt();
  }
  Rprintf("\n\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
  Rprintf("ABC-MCMC accepteance ratio was: %3.5f",
    static_cast<double>(nacc) / static_cast<double>(nabc));
  Rprintf("\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");

//  return List::create(Named("sim") = ans,
//                      Named("naccept") = nacc);
  return thinMat(ans, thin);
}

//[[Rcpp::export]]
double
distPsi_rls(arma::vec mu_lsig, arma::vec y, const double c1, const double c2, const double kc2) {
  arma::vec psi = funPsi_rls(mu_lsig(0), exp(mu_lsig(1)), y, y.n_elem, c1, c2, kc2);
  return as_scalar(trans(psi)*psi);
}

//[[Rcpp::export]]
double post_ls(double mu,
              double lsig,
              arma::vec y,
              double muSD,
              double sigScale,
              bool lg = false)
{
  int n = y.n_elem;
  double ans = 0.0;
  for(int i=0; i<n; i++){
    ans += Rf_dnorm4(y(i), mu, exp(lsig), true);
  }
  ans += dPrior(mu, lsig, muSD, sigScale, true);
  if (lg == false){
      ans = exp(ans);
      }
  return ans;
}


