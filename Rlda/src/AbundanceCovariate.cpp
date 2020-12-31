#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "progress.hpp"
#include <iostream>
#include <ctime>
#include <fstream>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]
using namespace arma;
using namespace Rcpp;



/***************************************************************************************************************************/
/*********************************                      UTILS          *****************************************************/
/***************************************************************************************************************************/

arma::mat mvrnorm(int n, arma::vec mu, arma::mat Sigma){
  int p = Sigma.n_cols;
  arma::mat X = reshape(arma::vec(rnorm(p * n)), p, n);
  arma::vec eigval;
  arma::mat eigvec;
  arma::eig_sym(eigval, eigvec, Sigma);
  X = eigvec * arma::diagmat(sqrt(eigval)) * X;
  X.each_col() += mu;
  return(X.t());
}

IntegerVector getznew(arma::mat ynew, int n, int ncomm){
  IntegerVector znew(n);

  for (int i = 0; i < n; i++){
    double max1=-1 * std::numeric_limits<double>::max();
    for (int j = 0; j < ncomm; j++){
      if (ynew(i,j)> max1) {
        max1=ynew(i,j);
        znew[i]=j+1;
      }
    }
  }
  return znew;
}

IntegerVector generatezz(arma::mat logphi, int n, IntegerVector znew, IntegerVector zold,
                        IntegerVector w, NumericVector unifs) {
  double pnew=0;
  double pold=0;
  double k;
  for (int i = 0; i < n; i++) {
    pnew=logphi(znew[i]-1,w[i]-1);
    pold=logphi(zold[i]-1,w[i]-1);

    k=exp(pnew-pold);
    if (unifs[i]<k){
      zold[i]=znew[i];
    }
  }
  return zold;
}

arma::rowvec rdirichletArma(Rcpp::NumericVector parms) {
  arma::rowvec res(parms.size());
  double sample_sum = 0;
  for(int j=0;j<parms.size();j++){
    res(j) = R::rgamma(parms[j], 1);
    sample_sum += res(j);
  }
  for(int j = 0; j<parms.size();j++) {
    res(j) = res(j) / sample_sum ;
  }
  return (res);
}

arma::mat rdirichletArma(int n, Rcpp::NumericVector parms) {
  arma::mat resMat(n,parms.size());
  for(int i=0;i<n;i++){
    arma::rowvec res(parms.size());
    double sample_sum = 0;
    for(int j=0;j<parms.size();j++){
      res(j) = R::rgamma(parms[j], 1);
      sample_sum += res(j);
    }
    for(int j = 0; j<parms.size();j++) {
      res(j) = res(j) / sample_sum ;
    }
    resMat.row(i) =res;
  }
  return (resMat);
}

/***************************************************************************************************************************/
/*********************************            GIBBS SAMPLING FUNCTIONS           *******************************************/
/***************************************************************************************************************************/

arma::mat generateBetas(List resYZW, arma::mat xMat, arma::vec varBetas){
  //Get the yMat
  arma::mat yMat = resYZW(0);
  //Number of communities
  int nCommunity = yMat.n_cols;
  //Number of covariates
  int nCovariates = xMat.n_cols;
  //Initialize the beta matrix
  arma::mat betaMat = arma::mat(nCommunity,nCovariates);
  //Calculate t(X)%*%X
  arma::mat tXXmat = xMat.t()*xMat;

  //Initialize the generation
  for(int c=0;c<(nCommunity-1);c++){
    //Get the response variable
    arma::vec yVec = yMat.col(c);
    //Get the precision
    arma::mat precMat = tXXmat+arma::diagmat(1.0/varBetas);
    //Get the inverse
    arma::mat varMat = arma::inv(precMat);
    //Get the mean
    arma::vec vecMean = varMat*xMat.t()*yVec;
    arma::mat temp = mvrnorm(1,vecMean,varMat);
    //Generate beta
    arma::rowvec tempVec = temp.row(0);
    betaMat.row(c)=tempVec;
  }

  //For the last group, the intercept has to be zero to ensure parameter identifiability:
  int cLast = nCommunity-1;
  //Get the response variable
  arma::vec yVec = yMat.col(cLast);
  //Remove the first column
  arma::mat xMatTemp = xMat.cols(1,xMat.n_cols-1);
  //Remove the first varBeta
  arma::vec varBetaTemp = varBetas.subvec(1,xMat.n_cols-1);
  //Get the precision
  arma::mat precMat = xMatTemp.t()*xMatTemp+arma::diagmat(1.0/varBetaTemp);
  //Get the inverse
  arma::mat varMat = arma::inv(precMat);
  //Get the mean
  arma::vec vecMean = varMat*xMatTemp.t()*yVec;
  //Generate beta
  arma::mat temp = mvrnorm(1,vecMean,varMat);
  //Create the zero matriz
  arma::mat zero(1,1);
  zero.fill(0.0);
  //Join the matrix
  arma::mat temp2 = arma::join_rows(zero,temp);
  arma::rowvec tempVec = temp2.row(0);
  betaMat.row(cLast)= tempVec;
  return(betaMat);
}

List generateZ(List resYZW, arma::mat xMat, arma::mat betasMat, arma::mat phiMat){
  //Get the yMat
  arma::mat yMat = resYZW(0);
  //Get the zOld
  IntegerVector zOld = resYZW(1);
  //Get the wVector
  IntegerVector wVector = resYZW(2);
  //Get the total number of elements
  int nElements = yMat.n_rows;
  //Number of columns
  int nCol = yMat.n_cols;
  //Generate Uniform
  NumericVector vecUnif = Rcpp::runif(nElements,0,1);
  //Estimate the mean matrix
  arma::mat matMean = xMat*betasMat.t();
  //Generate the Y's
  arma::mat newY(nElements,nCol);
  for(int r=0;r<nElements;r++){
    for(int c=0;c<nCol;c++){
      newY(r,c) = R::rnorm(matMean(r,c),1.0);
    }
  }

  //Get the Z matrix
  IntegerVector zNew = getznew(newY, nElements, nCol);

  //Log-Phi
  arma::mat logPhi = arma::log(phiMat);

  //Acept ou Reject the Y's and Z's
  IntegerVector zVec = generatezz(logPhi, nElements, zNew, zOld, wVector, vecUnif);
  //Compare the results
  for(int i=0;i<zVec.length();i++){
    if(zVec(i)!=zNew(i)){
      newY.row(i) = yMat.row(i);
    }
  }

  //'Store the results
  List resTemp = Rcpp::List::create(Rcpp::Named("Y") = newY,
                                    Rcpp::Named("Z")  = zVec ,
                                    Rcpp::Named("W")  = wVector);

  return resTemp;
}


arma::mat generatePhi(List resYZW, arma::mat xMat, double gamma, int nSpecies){
  //Get the yMat
  arma::mat yMat = resYZW(0);
  //Get the zOld
  IntegerVector zVector = resYZW(1);
  //Get the wVector
  IntegerVector wVector = resYZW(2);
  //Number of communities
  int nCommunity = yMat.n_cols;
  //Initialize the Phi matrix
  arma::mat phiMat(nCommunity,nSpecies);
  for(int c=0;c< nCommunity;c++){
    //Create the NumericVector count
    NumericVector vecCount(nSpecies);
    for(int i=0;i<zVector.length();i++){
      if(zVector(i)==(c+1)){
        int iComm = wVector(i)-1;
        vecCount(iComm)= vecCount(iComm)+1;
      }
    }
    //Generate the Dirichlet
    NumericVector parms = vecCount+gamma;
    phiMat.row(c) = rdirichletArma(parms);
  }
  return(phiMat);
}


/***************************************************************************************************************************/
/*********************************            GIBBS SAMPLING PROCEDURE                  ************************************/
/***************************************************************************************************************************/



List lda_multinomial_cov(DataFrame yData, DataFrame xData, IntegerVector specVector, arma::vec varBetas, int n_community, int n_specie, double gamma, int n_gibbs, bool ll_prior=true, bool display_progress=true) {

  //'Convert to matrix
  NumericMatrix yMatNM = internal::convert_using_rfunction(yData, "as.matrix");
  arma::mat yMat = as<arma::mat>(wrap(yMatNM));

  NumericMatrix xMatNM = internal::convert_using_rfunction(xData, "as.matrix");
  arma::mat xMat = as<arma::mat>(wrap(xMatNM));

  //Initialize the Community Vector
  IntegerVector comVector(yMat.n_rows);
  comVector.fill(1);

  //Create the list object
  List listYZW = Rcpp::List::create(Rcpp::Named("Y")  = yMat,
                                    Rcpp::Named("Z")  = comVector ,
                                    Rcpp::Named("W")  = specVector);

  //'Intialize Phi
  NumericVector hyperPhi(n_specie);
  hyperPhi.fill(1);
  arma::mat Phi = rdirichletArma(n_community, hyperPhi);

  //'Initialize Betas
  arma::mat Betas = arma::mat(n_community,xMat.n_cols);
  Betas.fill(1.0);

  //'Initialize the logLikelihood vector
  NumericVector logLikelihoodVec(n_gibbs);

  //'Initialize the results
  List listResults(n_gibbs);

  //'Intialize the progressbar
  Progress p(n_gibbs, display_progress);
  for (int g = 0; g < n_gibbs; ++g) {
    //'Verify if everything is ok
    if (Progress::check_abort() )
      Rcpp::stop("Operation cancelled by interrupt.");

    //Generate listYZW
    listYZW = generateZ(listYZW, xMat, Betas, Phi);

    //Generate Betas
    Betas = generateBetas(listYZW, xMat, varBetas);

    //Generate Phi
    generatePhi(listYZW, xMat, gamma, n_specie);

    List resObjects = Rcpp::List::create(Rcpp::Named("Betas") = Betas,
                                         Rcpp::Named("Phi")  = Phi );

    //Store the objects
    listResults(g) = resObjects;

    //'Increment the progress bar
    p.increment();

  }

  return listResults;
}
