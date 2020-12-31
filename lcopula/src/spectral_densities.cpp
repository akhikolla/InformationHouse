#include <Rcpp.h>
using namespace Rcpp;


//' Spectral density of scaled Dirichlet and negative Dirichlet extreme value distributions
//'
//' @param param vector containing the parameters \eqn{\alpha} and \eqn{\rho}, in this order.
//' @param dat matrix of Frechet or Pareto observations, of dimension \eqn{n} by \eqn{d}.
//' @param d dimension of model
//' @param transform logical indicating whether parameters are on the log scale. Default to \code{TRUE}
//'
//' @details The function is provided as a wrapper and takes parameters on the log scale for \eqn{\alpha} (\eqn{\rho}).
//' @return the log-likelihood for the \code{n} sample
// [[Rcpp::export(dirspecdens)]]
NumericVector dirspecdens(NumericVector param, NumericMatrix dat, int d, bool transform = true){
 //Containers for copies
  NumericVector rho(1);
  NumericVector alphavec(d);

  if((d+1) == param.size()){
    if(transform){
      rho[0] = exp(param(d));
    } else{
      rho[0] = param(d);
    }
  } else if(d == param.size()){
    rho[0] = 1.0;
  }
  if(transform){
    alphavec = exp(param[(seq_len(d)-1)]);
  } else{
    alphavec = param[(seq_len(d)-1)];
  }
  //Check parameter inputs are not aberrant
  if(is_true(any(alphavec>100)) || is_true(any(alphavec<0)) || rho[0]<0 || rho[0]>50){
    return(NumericVector::create(-1e10));
  }
  //Container for result
  NumericVector dens = NumericVector::create(0.0);
  NumericVector proddens(1);
  NumericVector alpha = NumericVector::create(sum(alphavec));
  //Test for dimension match
  if(d != dat.ncol()){
    Rcpp::stop("Invalid argument; input matrix should be such that the row of the d-columns sums to one.");
  }
  NumericVector ckons =  lgamma(alphavec+rho[0])-lgamma(alphavec);
  //Loop over data entries
  for(int j=0; j<dat.nrow(); j++){
    proddens[0] = 0.0;
    for(int k=0; k<d; k++){
      dens[0] = dens[0] + alphavec[k]/rho[0]*ckons[k]+(alphavec[k]/rho[0]-1)*log(dat(j,k));
      proddens[0] += exp((ckons[k]+log(dat(j,k)))/rho[0]);
    }

    dens[0] = dens[0] -(rho[0]+alpha[0])*log(proddens[0]);
  }
  dens[0] = dens[0]+dat.nrow()*(lgamma(alpha[0]+rho[0])-(d-1)*log(rho[0])-sum(lgamma(alphavec))-log(static_cast<double>(d)));

  return dens;
}

//' Spectral density of scaled negative Dirichlet extreme value distribution
//'
//' @param param vector containing the parameters \eqn{\alpha} and \eqn{\rho}, in this order.
//' @param dat matrix of Frechet or Pareto observations, of dimension \eqn{n} by \eqn{d}.
//' @param d dimension of model
//' @param transform logical indicating whether parameters are on the log scale. Default to \code{TRUE}
//'
//' @details The function is provided as a wrapper and takes parameters on the log scale for \eqn{\alpha} and \eqn{\rho}.
//' @return the log-likelihood for the \code{n} sample
// [[Rcpp::export(negdirspecdens)]]
NumericVector negdirspecdens(NumericVector param, NumericMatrix dat, int d, bool transform = true){
  //Containers for copies
  NumericVector rho(1);
  NumericVector alphavec(d);
  //Checks for valid inputs (rho must be less than min (alpha)
  if(min(param)!=param[d]){
    return(NumericVector::create(-1e10));
  }
  if((d+1) == param.size()){
    if(transform){
      rho[0] = exp(param[d]);///(1.0+exp(param(d)));
    } else{
      rho[0] = param[d];
    }
  } else if(d == param.size()){
    rho[0] = 1.0;
  }
  if(transform){
    alphavec = exp(param[(seq_len(d)-1)]);
  } else{
    alphavec = param[(seq_len(d)-1)];
  }
  //Check that transformed input makes sense
  if(is_true(any(alphavec>100)) || is_true(any(alphavec<=rho[0])) || rho[0]<0 || rho[0]>50){
    return(NumericVector::create(-1e10));
  }
  //Container for result
  NumericVector dens = NumericVector::create(0.0);
  NumericVector proddens(1);
  NumericVector alpha = NumericVector::create(sum(alphavec));
  //Test for dimension match
  if(d != dat.ncol()){
    Rcpp::stop("Invalid argument; dimensions of alpha vector and dat do not match");
  }
  NumericVector akons =  lgamma(alphavec-rho[0])-lgamma(alphavec);
  //Loop over data entries
  for(int j=0; j<dat.nrow(); j++){
    proddens[0] = 0.0;
    for(int k=0; k<d; k++){
      dens[0] = dens[0] - alphavec[k]/rho[0]*akons[k]+(-alphavec[k]/rho[0]-1)*log(dat(j,k));
      proddens[0] += exp(-(akons[k]+log(dat(j,k)))/rho[0]);
    }

    dens[0] = dens[0] -(alpha[0]-rho[0])*log(proddens[0]);
  }
  dens[0] = dens[0]+dat.nrow()*(lgamma(alpha[0]-rho[0])-(d-1)*log(rho[0])-sum(lgamma(alphavec))-log(static_cast<double>(d)));

  return dens;
}



//' Spectral density of Coles and Tawn extreme value distribution
//'
//' @param param vector containing the parameters \eqn{\alpha}.
//' @param dat matrix of pseudo-dat, of dimension \eqn{d}
//' @param transform logical indicating whether parameters are on the log scale. Default to \code{TRUE}
//' @return the value of the log likelihood for a sample of size \code{n}
// [[Rcpp::export(.ctspecdens)]]
NumericVector ctspecdens(NumericVector param, NumericMatrix dat, bool transform = true){
  int d = param.size();
  NumericVector alphavec = param;
  if(transform){
    alphavec = exp(param);
  }
  //Container for result
  NumericVector dens = NumericVector::create(0.0);
  NumericVector proddens(1);
  //Test for dimension match
  if(d !=dat.ncol()){
    Rcpp::stop("Invalid argument; dimensions of alpha vector and dat do not match");
  }
  //Test for aberrant values
  if(is_true(any(alphavec>50)) || is_true(any(alphavec<0))){
    return(NumericVector::create(-1e10));
  }
  double alpha = sum(alphavec);
  NumericVector c = alphavec/alpha;
  NumericVector Constant = NumericVector::create(lgamma(alpha)-sum(lgamma(alphavec))-log(static_cast<double>(d)));
  for(int j=0; j< dat.nrow(); j++){
    proddens[0] = 0.0;
    for(int k=0; k<d; k++){
      proddens[0] += alphavec[k]*log(c[k])+(alphavec[k]-1)*log(dat(j,k));
    }
    dens[0] = dens[0] +Constant[0]-
      (1.0+alpha)*log(sum(exp(log(c)+log(dat(j,_)))))+proddens[0];
  }
  return dens;
}
