/******************  created by Erlis Ruli 14/10/2015   ************************
 ****************** Robust Approximate Bayesian Infernece  *********************
 ******************************************************************************/

#include <RcppArmadillo.h>

// [[Rcpp::export]]
double
  psi_huber(double x, const double c)
  {
    // Huber's psi = rho'()
    return (x <= -c) ? -c : ((x < c) ? x : c);
  }

//[[Rcpp::export]]
arma::vec vpsi_huber(arma::vec x, const double c, int xLen)
{
  // Huber's psi = rho'()
  arma::vec ans(xLen);
  for(int i=0; i<xLen; i++){
    ans(i) = psi_huber(x(i), c);
  }
  return ans;
}

//[[Rcpp::export]]
double
  psip_huber(double x, const double c)
  {
    // psi' = rho'' : Second derivative of Huber's loss function
    return (fabs(x) >= c) ? 0. : 1.;
  }

// [[Rcpp::export]]
double
  dhalfCauchy(double x, double scale, bool lg = false)
  {
    double dens = log(2 * scale) - log(arma::datum::pi*(x*x + scale*scale));

    if (lg == false)
        dens = exp(dens);

    return dens;
  }


//[[Rcpp::export]]
double dinvgamma (double x, double shape, double scale, bool lg)
{

  double lden = shape * log(scale) - Rf_lgammafn(shape) -
    (shape + 1) * log(x) - (scale/x);

  if(lg == false) {
    lden = exp(lden);
  }

  return lden;
}


// [[Rcpp::export]]
arma::mat thinMat(arma::mat X, arma::vec index){
  int ncol = index.n_elem, nrow = X.n_rows;
  arma::mat ans(nrow, ncol);

  for(int i=0; i<ncol; i++){
    ans.col(i) = X.col(index(i));
  }

  return ans;
}

//[[Rcpp::export]]
arma::vec rmvnorm(arma::vec mu,
                  arma::mat S,
                  int p)
{
  /*
  Multivariate normal random variates Y = mu + AX
  X is iid standard normal
  mu mean vector and A is such that AA^T = S
  */

  arma::mat A;
  chol(A, S); // note: this give A^TA = S, with A upper-triangular
  Rcpp::RNGScope scope;
  bool status = A.is_empty();

  if (status==true) {
    Rprintf("\nCholesky decomposition in rmvnorm failed!");
    return arma::vec(arma::datum::nan);

  } else {

    arma::vec X(p);
    for(int i=0; i<p; i++){

      X(i) = Rf_rnorm(0.0, 1.0);

    }
    // arma::vec X = arma::randn(p);
    return mu + A.t()*X;
    }
}

//[[Rcpp::export]]
arma::vec rmvnorm2(arma::vec mu,
                   arma::mat lowcholS,
                   int p)
{
  /*
   Multivariate normal random variates Y = mu + lowcholSX
   X is iid standard normal
   mu mean vector and lowcholS is such that lowcholS*lowcholS^T = S
   */

  arma::vec X(p);
  Rcpp::RNGScope scope;

  for(int i=0; i<p; i++) {
    X(i) = Rf_rnorm(0.0, 1.0);
  }

  return mu + lowcholS*X;
}

//[[Rcpp::export]]
arma::vec rmvt(arma::vec mu,
               arma::mat lowcholS,
               int p,
               double df)
{
  double chi = Rf_rchisq(df);
  arma::vec mm(p, arma::fill::zeros);

  arma::vec z = rmvnorm2(mm, lowcholS, p);

  // bool status = z.is_empty();

  // if(status == true){
    // return arma::vec(arma::datum::nan);
    // } else {
      return mu + z*sqrt(df/chi);
      // }
}

//[[Rcpp::export]]
double dmvt(arma::vec x,
            arma::vec mu,
            arma::mat S,
            double ldetS,
            int p,
            double df,
            bool lg)
{
  double lans = 0.0;
  lans += Rf_lgammafn(0.5*(df+p)) - Rf_lgammafn(0.5*df) - 0.5*p*(log(df) +
    log(arma::datum::pi));
  // arma::log_det(ldet, sign, S);
  lans += -(0.5*ldetS + 0.5*(df+p)*log(1.0 + arma::as_scalar(trans(x-mu)
                                                                               *inv( S )*(x-mu))/df));
  if(lg == true){
    return lans;

  } else {

    return exp(lans);
  }

}


