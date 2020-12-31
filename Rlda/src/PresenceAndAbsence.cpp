#include <Rcpp.h>
#include "progress.hpp"
#include <iostream>
#include <ctime>
#include <fstream>
// [[Rcpp::depends(RcppProgress)]]
using namespace Rcpp;


/***************************************************************************************************************************/
/*********************************                      UTILS          *****************************************************/
/***************************************************************************************************************************/



double matchPresence(double value,NumericVector vec){
  for(int i=0;i<vec.length();i++){
    if(value==vec(i)){
      return(i);
    }
  }
  return(-1);
}

NumericVector rdirichletPresence(Rcpp::NumericVector parms) {
  NumericVector res(parms.size());
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

NumericMatrix rdirichletPresence(int n, Rcpp::NumericVector parms) {
  NumericMatrix res(n, parms.size());
  for(int i=0;i<n;i++){
    double sample_sum = 0;
    for(int j=0;j<parms.size();j++){
      res(i, j) = R::rgamma(parms[j], 1);
      sample_sum += res(i, j);
    }
    for(int j = 0; j<parms.size();j++) {
      res(i, j) = res(i, j) / sample_sum ;
    }
  }
  return (res);
}


int whichLessDVPresence(double value, NumericVector prob) {
  int res = -1;

  //'Create the categorical table
  double probcum = 0;

  for (int i = 0; i < prob.length(); i++) {
    probcum = probcum + prob(i);
    if (value < probcum) {
      res = i;
      break;
    }
  }
  return res;
}

NumericVector rmultinomialDVPresence(int size, NumericVector prob) {
  //'Initialize the NumericMatrix result
  NumericVector res(prob.length());
  //'Initialize the values
  for(int j=0;j<prob.length();j++){
    res(j)=0;
  }

  //'Generate the sample
  int count=0;
  //'Generate until get the sample size
  while(count<size){
    //'Draw a uniform
    double random = R::runif(0,1);
    //'Which category was draw ?
    int iPos=whichLessDVPresence(random,prob);
    //'Increment the matrix
    res(iPos)=res(iPos)+1;
    //'Increment the counter
    count=count+1;
  }

  return res;
}


NumericVector invertedCumsumPresence(NumericVector n){
  NumericVector table(n.length());
  table(n.length()-1)=n(n.length()-1);
  for(int i=(n.length()-2);i>-1;i--){
    table(i)=table(i+1)+n(i);
  }
  return(table);
}

NumericVector meltPresence(NumericMatrix mat){
  //'Initialize the NumericVector
  NumericVector vec(mat.nrow()*mat.ncol());
  //'Initialize the position
  int pos=0;
  for(int col=0;col<mat.ncol();col++){
    for(int row=0;row<mat.nrow();row++){
      //'meltPresence the matrix
      vec(pos)=mat(row,col);
      //'Increment the position
      pos=pos+1;
    }
  }
  return(vec);
}

void updatePhiAndThetaPresence(NumericMatrix &PhiGibbs,NumericMatrix Phi,NumericMatrix &ThetaGibbs, NumericMatrix Theta,int gibbs){
    //'meltPresence the Phi and Theta matrix
    PhiGibbs(gibbs,_)=meltPresence(Phi);
    ThetaGibbs(gibbs,_) =meltPresence(Theta);
}


double tnormPresenceAbsence(double lo, double hi,double mu, double sig){
  double q1 = R::pnorm5(lo,mu,sig,1,0);
  double q2 = R::pnorm5(hi,mu,sig,1,0);
  double z = R::runif(q1,q2);
  z = R::qnorm5(z,mu,sig,1,0);
  return(z);
}

double fixMHPresenceAbsence(double lo, double hi,double old1,double new1,double jump){
  double jold=R::pnorm5(hi,old1,jump,1,0)-R::pnorm5(lo,old1,jump,1,0);
  double jnew=R::pnorm5(hi,new1,jump,1,0)-R::pnorm5(lo,new1,jump,1,0);
  return(std::log(jold)-std::log(jnew));
}

/***************************************************************************************************************************/
/*********************************            GIBBS SAMPLING FUNCTIONS           *******************************************/
/***************************************************************************************************************************/

double gammaMHPresenceAbsence(NumericMatrix vMat, double gamma, double jump, int &acept){
  double newGamma = tnormPresenceAbsence(0.0,1.0,gamma,jump);
  double pold = 0.0;
  double pnew = 0.0;
  for(int r=0;r<vMat.nrow();r++){
    for(int c=0;c<(vMat.ncol()-1);c++){
      pold=pold+R::dbeta(vMat(r,c),1.0,gamma,1);
      pnew=pnew+R::dbeta(vMat(r,c),1.0,newGamma,1);
    }
  }
  double pcorrection = fixMHPresenceAbsence(0,1,gamma,newGamma,jump);
  //Acceptance
  double a = std::exp(pnew+pcorrection-pold);
  double z = R::unif_rand();
  if(z<a){
    acept = 1;
    return(newGamma);
  }
  else{
    acept = 0;
    return(gamma);
  }
  return(0);
}

List generateZPresence(NumericMatrix binaryMat, NumericMatrix Phi, NumericMatrix vMat, NumericMatrix &ThetaMat) {
  //'Number of locations
  int nLocations=binaryMat.nrow();
  //'Number of Species
  int nSpecies=binaryMat.ncol();
  //'Number of Communities
  int n_community=Phi.ncol();
  //'Matrix N the number of locations where specie s comes from community c
  NumericMatrix nMat(nSpecies,n_community);
  //'Matrix R the number of these individuals thar are observed
  NumericMatrix rMat(nSpecies,n_community);
  //'Matrix M the the number of species in plot l that come from community c
  NumericMatrix mMat(binaryMat.nrow(),n_community);
  //'For each location sample from a Multinomial
  for(int l=0;l<nLocations;l++){
    //'Create Theta (Size Cx1) based on vMat V_cl \prod_(k=1)^(c-1)(1-V_kl )
    NumericVector Theta(n_community);
    //'Update the Theta \prod_(k=1)^(c-1)(1-V_kl )
    double prod=1;
    for(int c=0;c<n_community;c++){
      double vNumber = vMat(l,c);
      if (c == 0) prod=1;
      if (c >  0) prod=prod*(1.0-vMat(l,c-1));
      Theta(c)=vNumber*prod;
    }
    //'Update the Theta Matrix by reference
    ThetaMat(l,_)=Theta;
    //'For each Specie
    for(int s=0;s<nSpecies;s++){
      //'Calculate the probability vector
      double sumVec=0;
      NumericVector probability = pow(Phi(s,_),binaryMat(l,s))*pow(1.0-Phi(s,_),1.0-binaryMat(l,s))*Theta;
      //'Find the sum
      sumVec=sum(probability);
      //'Normalize the probability
      probability=probability/sumVec;
      //'Presence in locatiol l and specie s (always will be a draw)
      int iSize=1;
      //'Store the results
      NumericVector tmp = rmultinomialDVPresence(iSize, probability);
      //'Find the community
      int iCommunity = (int) matchPresence(1.0,tmp);
      //'N matrix is the number of individuals in plot L that come from community C
      nMat(s,iCommunity)=nMat(s,iCommunity)+1.0;
      //'R matrix is the number of individuals in plot L that come from community C that are observed
      if(binaryMat(l,s)==1.0) rMat(s,iCommunity)=rMat(s,iCommunity)+1.0;
      //'M matrix is the number of species in plot L that come from community C
      mMat(l,iCommunity)=mMat(l,iCommunity)+1.0;
    }
  }

  //'Store the zMat
  List resZ = Rcpp::List::create(Rcpp::Named("N") = nMat,
                                 Rcpp::Named("R") = rMat,
                                 Rcpp::Named("M") = mMat);
  return resZ;
}

NumericMatrix generatePhiPresence(List zList,double alpha0, double alpha1) {
  //'Getting the N matrix
  NumericMatrix nMat = zList[0];
  //'Getting the R matrix
  NumericMatrix rMat = zList[1];
  //'Total number of species
  int nSpecies = nMat.nrow();
  //'Total number of communities
  int n_community = nMat.ncol();
  //'Initialize the Theta matrix
  NumericMatrix PhiMat(nSpecies,n_community);
  //'For each Specie
  for(int s=0;s<nSpecies;s++){
    //'For each community:
    for(int c=0;c<n_community;c++){
      double aBeta = rMat(s,c)+alpha0;
      double bBeta = nMat(s,c)-rMat(s,c)+alpha1;
      PhiMat(s,c)=R::rbeta(aBeta,bBeta);
    }
  }
  return PhiMat;
}


NumericMatrix generateVPresence(List zList,int nLocations, double gamma) {
  //'Getting the M matrix
  NumericMatrix mMat = zList[2];
  //'Total number of communities
  int n_community = mMat.ncol();
  //'Initialize the Theta matrix
  NumericMatrix vMat(nLocations,n_community);
  //'Foreach Specie
  for(int l=0;l<nLocations;l++){
    //'For each community:
    NumericVector nGreater = invertedCumsumPresence(mMat(l,_));
    for(int c=0;c<n_community;c++){
      //'nLC is the number of species in plot l that come from community c
      double nLC = mMat(l,c);
      if(c<n_community-1){
        //'Generate stick-breaking probabilities
        vMat(l,c)=R::rbeta(1.0+nLC,gamma+nGreater(c+1));
      }
      else{
        //'Ensure that the last community has 1
        vMat(l,c)=1.0;
      }
    }
  }
  return vMat;
}

double ll_priorFunctionPresence(NumericMatrix matDATA,int nLocations, int nSpecies,int n_community, NumericMatrix vMat,NumericMatrix Phi, NumericMatrix Theta, double alpha0,double alpha1, double gamma, bool ll_prior=true) {
  //'Initialize the logLikelihoodVec
  double logLikelihood=0;
  //'Calculate the Loglikelihood and Prior
  if(ll_prior){
    //'Initialize the V_{cl} and Phi_{sc} prior
    double priorV=0.0;
    double priorPhi=0.0;
    //'For each location
    for(int l=0;l<nLocations;l++){
      //'Compute the prior for V_{cl}
      for(int c=0;c<n_community;c++){
        if(vMat(l,c)<1)priorV=priorV+R::dbeta(vMat(l,c),1,gamma,1);
      }

      //'For each Specie
      for(int s=0;s<nSpecies;s++){
        //'Compute just one time the prior for Phi_{sc}
        if(l==0){
          //'Compute the prior for Phi_{sc}
          for(int c=0;c<n_community;c++){
            if(Phi(s,c)<1)priorPhi=priorPhi+R::dbeta(Phi(s,c),alpha0,alpha1,1);
          }
        }

        //'Initialize the product between Phi and Theta
        double PhiTheta=0.0;
        //'For each community
        for(int c=0;c<n_community;c++){
          PhiTheta=PhiTheta+Phi(s,c)*Theta(l,c);  //'DV: mudei aqui
        }
        //'Getting the x_{ls} observation
        double xLS=matDATA(l,s);
        if(xLS==1.0){
          logLikelihood=logLikelihood+log(PhiTheta);
        }
        else{
          if(PhiTheta<1){//'Otherwise a large PhiTheta>1 represents low likelihood, then we don't compute this value.
            logLikelihood=logLikelihood+log(1.0-PhiTheta);
          }
        }
      }
    }

    //'Update the logLikelihood
    logLikelihood=logLikelihood+priorPhi+priorV;
  }
  else{
    //'For each location
    for(int l=0;l<nLocations;l++){
      //'For each Specie
      for(int s=0;s<nSpecies;s++){
        //'Initialize the product between Phi and Theta
        double PhiTheta=0.0;
        //'For each community
        for(int c=0;c<n_community;c++){
          PhiTheta=PhiTheta+Phi(s,c)*Theta(l,c); //'DV: mudei aqui
        }
        //'Getting the x_{ls} observation
        double xLS=matDATA(l,s);
        if(xLS==1.0){
          logLikelihood=logLikelihood+log(PhiTheta);
        }
        else{
          if(PhiTheta<1){//'Otherwise a large PhiTheta>1 represents low likelihood, then we don't compute this value.
            logLikelihood=logLikelihood+log(1.0-PhiTheta);
          }
        }
      }
    }
  }

  return(logLikelihood);
}
/***************************************************************************************************************************/
/*********************************            GIBBS SAMPLING PROCEDURE                  ************************************/
/***************************************************************************************************************************/

//' @name GibbsSamplingPresence
//' @title Gibbs Sampling for LDA Presence and Absence
//' @description Compute the Gibbs Sampling for LDA Presence and Absence
//' @param DATA - DataFrame with Presence and Absecence (Zeros and Ones)
//' @param int n_community - Number of communities
//' @param alpha0 - Hyperparameter Beta(alpha0,alpha1)
//' @param alpha1 - Hyperparameter Beta(alpha0,alpha1)
//' @param gamma - Hyperparameter  Beta(1,gamma)
//' @param n_gibbs - Total number of Gibbs Samples
//' @param ll_prior - Likelihood compute with Priors ?
//' @param bool display_progress=true - Should I Show the progressBar ?
//' @return List - With Phi(n_gibbs,n_community*nSpecies), Theta(n_gibbs,nLocations*n_community) and logLikelihood
// [[Rcpp::export]]
List lda_bernoulli(DataFrame data, int n_community, double alpha0, double alpha1, double gamma, int n_gibbs, bool ll_prior=true, bool display_progress=true) {

  //'Convert to matrix
  NumericMatrix matDATA = internal::convert_using_rfunction(data, "as.matrix");

  //'Total number of locations
  int nLocations = matDATA.nrow();

  //'Total number of species
  int nSpecies = matDATA.ncol();

  //'Intialize Phi
  NumericVector hyperPhi(n_community);
  hyperPhi.fill(1);
  NumericMatrix Phi=rdirichletPresence(nSpecies,hyperPhi);

  //'Intialize vMat
  NumericVector hyperV(n_community);
  hyperV.fill(1);
  NumericMatrix vMat=rdirichletPresence(nLocations,hyperV);

  //'Initialize the PhiGibbs
  NumericMatrix PhiGibbs(n_gibbs,n_community*nSpecies);

  //'Initialize the ThetaGibbs
  NumericMatrix ThetaGibbs(n_gibbs,nLocations*n_community);

  //'Initialize the logLikelihood vector
  NumericVector logLikelihoodVec(n_gibbs);

  //Check if gamma existis
  bool bgamma = false;
  if(std::isnan(gamma)){
    bgamma=true;
    gamma = 0.01;
  }
  //Define the initial jump
  double jump = 0.5;
  int acept = 0;

  //'Intialize the progressbar
  Progress p(n_gibbs, display_progress);
  for (int g = 0; g < n_gibbs; ++g) {
    //'Verify if everything is ok
    if (Progress::check_abort() )
      Rcpp::stop("Operation cancelled by interrupt.");
    //'Initialize the Theta matrix
    NumericMatrix ThetaMat(nLocations, n_community);

    //'Generate zList
    List zList  = generateZPresence(matDATA, Phi, vMat, ThetaMat);

    //Generate gamma MH
    if(bgamma){
      if ((g%50==0) & (g<500)){
        double z = acept/50;
        if ((z>0.4) & (jump<100))   jump=jump*2;
        if ((z<0.1) & (jump>0.001)) jump=jump*0.5;
        gamma = gammaMHPresenceAbsence(vMat, gamma, jump,acept);
      }
    }

    //'Generate Phi
    Phi = generatePhiPresence(zList,alpha0,alpha1);

    //'Generate vMat
    vMat = generateVPresence(zList,nLocations, gamma);

    //'Create the final Phi (n_gibbs,n_community*nSpecies) and final Theta (ThetaGibbs)
    updatePhiAndThetaPresence(PhiGibbs, Rcpp::transpose(Phi), ThetaGibbs, ThetaMat, g);

    //'Initialize the logLikelihood
    double logLikelihood=ll_priorFunctionPresence(matDATA, nLocations,
                                                       nSpecies,n_community,
                                                       vMat, Phi, ThetaMat,
                                                       alpha0, alpha1, gamma,
                                                       ll_prior);
    //'Store the logLikelihood
    logLikelihoodVec(g)=logLikelihood;

    //'Increment the progress bar
    p.increment();

  }

  //'Store the results
  List resTemp = Rcpp::List::create(Rcpp::Named("Theta") = ThetaGibbs,
                                    Rcpp::Named("Phi")  = PhiGibbs,
                                    Rcpp::Named("logLikelihood")  =logLikelihoodVec);

  return resTemp;
}



//' @name GibbsSamplingPresence
//' @title Gibbs Sampling for LDA Presence and Absence
//' @description Compute the Gibbs Sampling for LDA Presence and Absence
//' @param DATA - DataFrame with Presence and Absecence (Zeros and Ones)
//' @param int n_community - Number of communities
//' @param alpha0 - Hyperparameter Beta(alpha0,alpha1)
//' @param alpha1 - Hyperparameter Beta(alpha0,alpha1)
//' @param gamma - Hyperparameter  Beta(1,gamma)
//' @param n_gibbs - Total number of Gibbs Samples
//' @param n_burn - Number of elements to burn-in
//' @param ll_prior - Likelihood compute with Priors ?
//' @param bool display_progress=true - Should I Show the progressBar ?
//' @return List - With Phi(n_gibbs,n_community*nSpecies), Theta(n_gibbs,nLocations*n_community) and logLikelihood
List lda_bernoulli_burn(DataFrame data, int n_community, double alpha0, double alpha1, double gamma, int n_gibbs,int n_burn, bool ll_prior=true, bool display_progress=true) {

  //'Convert to matrix
  NumericMatrix matDATA = internal::convert_using_rfunction(data, "as.matrix");

  //'Total number of locations
  int nLocations = matDATA.nrow();

  //'Total number of species
  int nSpecies = matDATA.ncol();

  //'Intialize Phi
  NumericVector hyperPhi(n_community);
  hyperPhi.fill(1);
  NumericMatrix Phi=rdirichletPresence(nSpecies,hyperPhi);

  //'Intialize vMat
  NumericVector hyperV(n_community);
  hyperV.fill(1);
  NumericMatrix vMat=rdirichletPresence(nLocations,hyperV);

  //'Initialize the PhiGibbs
  NumericMatrix PhiGibbs(n_gibbs-n_burn,n_community*nSpecies);

  //'Initialize the ThetaGibbs
  NumericMatrix ThetaGibbs(n_gibbs-n_burn,nLocations*n_community);

  //'Initialize the logLikelihood vector
  NumericVector logLikelihoodVec(n_gibbs-n_burn);
  int cont=0;

  //Check if gamma existis
  bool bgamma = false;
  if(std::isnan(gamma)){
    bgamma=true;
    gamma = 0.01;
  }
  //Define the initial jump
  double jump = 0.5;
  int acept = 0;

  //'Intialize the progressbar
  Progress p(n_gibbs, display_progress);
  for (int g = 0; g < n_gibbs; ++g) {
    //'Verify if everything is ok
    if (Progress::check_abort() )
      Rcpp::stop("Operation cancelled by interrupt.");
    //'Initialize the Theta matrix
    NumericMatrix ThetaMat(nLocations, n_community);

    //'Generate zList
    List zList  = generateZPresence(matDATA, Phi, vMat, ThetaMat);

    //Generate gamma MH
    if(bgamma){
      if ((g%50==0) & (g<500)){
        double z = acept/50;
        if ((z>0.4) & (jump<100))   jump=jump*2;
        if ((z<0.1) &( jump>0.001)) jump=jump*0.5;
        gamma = gammaMHPresenceAbsence(vMat, gamma, jump,acept);
      }
    }

    //'Generate Phi
    Phi = generatePhiPresence(zList,alpha0,alpha1);

    //'Generate vMat
    vMat = generateVPresence(zList,nLocations, gamma);
    if(g>n_burn){

      //'Create the final Phi (n_gibbs,n_community*nSpecies) and final Theta (ThetaGibbs)
      updatePhiAndThetaPresence(PhiGibbs, Rcpp::transpose(Phi), ThetaGibbs, ThetaMat, cont);

      //'Initialize the logLikelihood
      double logLikelihood=ll_priorFunctionPresence(matDATA, nLocations,
                                                    nSpecies,n_community,
                                                    vMat, Phi, ThetaMat,
                                                    alpha0, alpha1, gamma,
                                                    ll_prior);
      //'Store the logLikelihood
      logLikelihoodVec(cont)=logLikelihood;
      cont=cont+1;
    }

    //'Increment the progress bar
    p.increment();

  }

  //'Store the results
  List resTemp = Rcpp::List::create(Rcpp::Named("Theta") = ThetaGibbs,
                                    Rcpp::Named("Theta")  = PhiGibbs ,
                                    Rcpp::Named("logLikelihood")  =logLikelihoodVec);

  return resTemp;
}
