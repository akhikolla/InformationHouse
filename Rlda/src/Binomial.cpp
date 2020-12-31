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

NumericVector matchBinomialResult(double value,NumericVector vec){
  int cont=0;
  for(int i=0;i<vec.length();i++){
    if(vec(i)!=value){
      cont=cont+1;
    }
  }
  NumericVector res(cont);
  cont=0;
  for(int i=0;i<vec.length();i++){
    if(vec(i)!=value){
      res(cont) = vec(i);
      cont=cont+1;
    }
  }
  return(res);
}
NumericVector matchBinomialIndex(double value,NumericVector vec){
  int cont=0;
  for(int i=0;i<vec.length();i++){
    if(vec(i)!=value){
      cont=cont+1;
    }
  }
  NumericVector res(cont);
  cont=0;
  for(int i=0;i<vec.length();i++){
    if(vec(i)!=value){
      res(cont) = i;
      cont=cont+1;
    }
  }
  return(res);
}

NumericVector rdirichletBinomial(Rcpp::NumericVector parms) {
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

NumericMatrix rdirichletBinomial(int n, Rcpp::NumericVector parms) {
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


int whichLessDVBinomial(double value, NumericVector prob) {
  int res = -1;

  //'Create the categorical table
  double probcum = 0;

  for (int i = 0; i < prob.length(); i++) {
    probcum = probcum + prob(i);
    if (value <= probcum) {
      res = i;
      break;
    }
  }
  return res;
}

NumericVector rmultinomialDVBinomial(int size, NumericVector prob) {
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
    int iPos=whichLessDVBinomial(random,prob);
    //'Increment the matrix
    res(iPos)=res(iPos)+1;
    //'Increment the counter
    count=count+1;
  }

  return res;
}

NumericVector invertedCumsumBinomial(NumericVector n){
  NumericVector table(n.length());
  table(n.length()-1)=n(n.length()-1);
  for(int i=(n.length()-2);i>-1;i--){
    table(i)=table(i+1)+n(i);
  }
  return(table);
}

NumericVector meltBinomial(NumericMatrix mat){
  //'Initialize the NumericVector
  NumericVector vec(mat.nrow()*mat.ncol());
  //'Initialize the position
  int pos=0;
  for(int col=0;col<mat.ncol();col++){
    for(int row=0;row<mat.nrow();row++){
      //'meltBinomial the matrix
      vec(pos)=mat(row,col);
      //'Increment the position
      pos=pos+1;
    }
  }
  return(vec);
}

void updateThetaAndPhiBinomial(NumericMatrix &ThetaGibbs,NumericMatrix Theta,NumericMatrix &PhiGibbs, NumericMatrix Phi,int gibbs){
  //'meltBinomial the Theta and Phi matrix
  ThetaGibbs(gibbs,_)=meltBinomial(Theta);
  PhiGibbs(gibbs,_) =meltBinomial(Phi);
}

NumericMatrix aggregateValuesBinomial(GenericVector zMat,int nLocations,int n_community){
  //'Initialize the matrix
  NumericMatrix sum(nLocations,n_community);
  //'Fill the sum matrix
  sum.fill(0.0);
  for(int i=0;i<zMat.size();i++){
    //Get the List
    List zMatTemp =  zMat[i];

    //'Get the community
    NumericVector iCommunity = zMatTemp[3];
    NumericVector iResults = zMatTemp[4];

    //'Get the location
    int iLocation = zMatTemp[1];

    //'Agregate
    for(int l=0;l<iCommunity.size();l++){
      sum(iLocation,iCommunity(l))=sum(iLocation,iCommunity(l))+iResults(l);
    }
  }
  //'Return the results
  return(sum);
}

List aggregateValuesBinomialByReflectanceBinomial(GenericVector zMat,int nBands,int n_community){
  //'Initialize the results
  List res(2);
  //'Initialize the matrix with Refletance
  NumericMatrix sumRef(nBands,n_community);
  //'Fill the sum matrix
  sumRef.fill(0.0);
  //'Initialize the matrix without Refletance
  NumericMatrix sumNonRef(nBands,n_community);
  //'Fill the sum matrix
  sumNonRef.fill(0.0);
  for(int i=0;i<zMat.size();i++){
    //Get the List
    List zMatTemp = zMat[i];
    //'Get the community
    NumericVector iCommunity = zMatTemp[3];
    NumericVector iResults = zMatTemp[4];
    //'Get the band
    int iBand = zMatTemp[0];
    //'Get the reflectance
    int iRefletance = zMatTemp[2];
    if(iRefletance==1){
      for(int l=0;l<iCommunity.size();l++){
        //'Agregate
        sumRef(iBand,iCommunity(l))=sumRef(iBand,iCommunity(l))+iResults(l);
      }
    }
    else{
      for(int l=0;l<iCommunity.size();l++){
        //'Agregate
        sumNonRef(iBand,iCommunity(l))=sumNonRef(iBand,iCommunity(l))+iResults(l);
      }
    }
  }
  //'Return the results
  res[0]=sumRef;
  res[1]=sumNonRef;
  return(res);

}

double sumLargestBinomial(NumericMatrix nSum, int iCommunity, int iLocation){
  double sum=0.0;
  //'For each community
  for(int c=(iCommunity+1);c<nSum.ncol();c++){
      sum=sum+nSum(iLocation,c);
  }
  return(sum);
}

NumericMatrix mmultBinomial(const NumericMatrix& m1, const NumericMatrix& m2){
  if (m1.ncol() != m2.nrow()) stop ("Incompatible matrix dimensions");
  NumericMatrix out(m1.nrow(),m2.ncol());
  NumericVector rm1, cm2;
  for (int i = 0; i < (int) m1.nrow(); ++i) {
    rm1 = m1(i,_);
    for (int j = 0; j < (int) m2.ncol(); ++j) {
      cm2 = m2(_,j);
      out(i,j) = std::inner_product(rm1.begin(), rm1.end(), cm2.begin(), 0.);
    }
  }
  return out;
}

double tnormBinomial(double lo, double hi,double mu, double sig){
  double q1 = R::pnorm5(lo,mu,sig,1,0);
  double q2 = R::pnorm5(hi,mu,sig,1,0);
  double z = R::runif(q1,q2);
  z = R::qnorm5(z,mu,sig,1,0);
  return(z);
}

double fixMHBinomial(double lo, double hi,double old1,double new1,double jump){
  double jold=R::pnorm5(hi,old1,jump,1,0)-R::pnorm5(lo,old1,jump,1,0);
  double jnew=R::pnorm5(hi,new1,jump,1,0)-R::pnorm5(lo,new1,jump,1,0);
  return(std::log(jold)-std::log(jnew));
}

/***************************************************************************************************************************/
/*********************************            GIBBS SAMPLING FUNCTIONS           *******************************************/
/***************************************************************************************************************************/

double gammaMHBinomial(NumericMatrix vMat, double gamma, double jump, int &acept){
  double newGamma = tnormBinomial(0.0,1.0,gamma,jump);
  double pold = 0.0;
  double pnew = 0.0;
  for(int r=0;r<vMat.nrow();r++){
    for(int c=0;c<(vMat.ncol()-1);c++){
      pold=pold+R::dbeta(vMat(r,c),1.0,gamma,1);
      pnew=pnew+R::dbeta(vMat(r,c),1.0,newGamma,1);
    }
  }
  double pcorrection = fixMHBinomial(0,1,gamma,newGamma,jump);
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

GenericVector generateZBinomial(NumericMatrix binomMat,NumericMatrix populMat, NumericMatrix Theta, NumericMatrix Phi) {
  try{
    //'Number of locations
    int nLocations=binomMat.nrow();
    //'Number of Bands
    int nBands=binomMat.ncol();
    //'Create the Z matrix
    GenericVector zMat;
    //'For each band
    for(int b=0;b<nBands;b++){
      //'For each location
      for(int l=0;l<nLocations;l++){

        //Create the Temporary List
        List zMatTemp = Rcpp::List::create(Rcpp::Named("Band"),
                                           Rcpp::Named("Location"),
                                           Rcpp::Named("Reflectance"),
                                           Rcpp::Named("Community"),
                                           Rcpp::Named("Latent"));

        //'Two cases: Reflectance and Not Reflectance
        //'Calculate the probability of belonging to each community:
        NumericVector prob=Phi(b,_)*Theta(l,_);
        prob=prob/Rcpp::sum(prob);

        if(binomMat(l,b)!=0){
          //'Size
          int iSize=binomMat(l,b);
          //'Store the results
          NumericVector tmp = rmultinomialDVBinomial(iSize, prob);

          //'Find the community
          NumericVector iCommunity = matchBinomialIndex(0.0,tmp);
          NumericVector iResults =   matchBinomialResult(0.0,tmp);

          //'Store the results
          zMatTemp["Band"]=b;             //'Store the band
          zMatTemp["Location"]=l;         //'Store the location
          zMatTemp["Reflectance"]=1;      //'Store the reflectance
          zMatTemp["Community"]=iCommunity;    //'Store the community
          zMatTemp["Latent"]=iResults; //'Store the size presence
          zMat.push_back(zMatTemp);
        }
        //'Calculate the probability of belonging to each community:
        prob=(1.0-Phi(b,_))*Theta(l,_);
        prob=prob/Rcpp::sum(prob);

        //Create the Temporary List
        List zMatTemp2 = Rcpp::List::create(Rcpp::Named("Band"),
                                           Rcpp::Named("Location"),
                                           Rcpp::Named("Reflectance"),
                                           Rcpp::Named("Community"),
                                           Rcpp::Named("Latent"));

        if(populMat(l,b)-binomMat(l,b)!=0){
          //'Size
          int iSize=populMat(l,b)-binomMat(l,b);
          //'Store the results
          NumericVector tmp = rmultinomialDVBinomial(iSize, prob);

          //'Find the community
          NumericVector iCommunity = matchBinomialIndex(0.0,tmp);
          NumericVector iResults =   matchBinomialResult(0.0,tmp);

          //'Store the results
          zMatTemp2["Band"]=b;             //'Store the band
          zMatTemp2["Location"]=l;             //'Store the location
          zMatTemp2["Reflectance"]=0;             //'Store the not reflectance
          zMatTemp2["Community"]=iCommunity;    //'Store the community
          zMatTemp2["Latent"]=iResults; //'Store the size absence
          zMat.push_back(zMatTemp2);
        }
      }
    }
    return zMat;
  }
  catch(std::exception &ex) {
    forward_exception_to_r(ex);
  } catch(...) {
    ::Rf_error("c++ exception (unknown reason)");
  }
  return 0;
}

NumericMatrix generateThetaBinomial(GenericVector zMat, NumericMatrix& vMat, int nLocations,int n_community, double gamma) {
  //'Intialize the Theta matrix
  NumericMatrix Theta(nLocations,n_community);
  //'Calculate the number of individuals in each community and locations
  NumericMatrix sumMat = aggregateValuesBinomial(zMat, nLocations, n_community);

  //'For each location
  for(int l=0;l<nLocations;l++){
    //'For each community
    for(int c=0;c<n_community;c++){
      if(c<(n_community-1)){
        //'How many elements belong to community c in location l
        double nLC=sumMat(l,c);
        //'How many elements are larger than community c
        double nLCgreter=sumLargestBinomial(sumMat, c, l);
        vMat(l,c)=R::rbeta(nLC+1.0,nLCgreter+gamma);
      }
      else{
        //'All locations for the last community are one
        vMat(l,c)=1;
      }
    }
  }
  //'Create the Theta matrix
  //'Foreach location
  for(int l=0;l<nLocations;l++){
    NumericVector thetaVec(n_community);
    //'Update the Theta prod_(k=1)^(c-1)(1-V_kl )
    double prod=1;
    //'For each community
    for(int c=0;c<n_community;c++){
      double vNumber = vMat(l,c);
      if (c == 0) prod=1;
      if (c >  0) prod=prod*(1.0-vMat(l,c-1));
      thetaVec(c)=vNumber*prod;
    }
    //'Store each row
    Theta(l,_)=thetaVec;
  }
  return(Theta);
}


NumericMatrix generatePhiBinomial(GenericVector zMat,int nBands,int n_community,double alpha0, double alpha1) {
  //'Initialize the Phi matrix
  NumericMatrix Phi(nBands,n_community);
  //'Get the sum matrices
  List sumList = aggregateValuesBinomialByReflectanceBinomial(zMat,nBands,n_community);
  NumericMatrix sumRef= sumList[0];
  NumericMatrix sumNonRef=sumList[1];
  //'Generate the Phi
  for(int b=0;b<nBands;b++){
    //'For each community
    for(int c=0;c<n_community;c++){
      Phi(b,c)=R::rbeta(alpha0+sumRef(b,c),alpha1+sumNonRef(b,c));
    }
  }
  return(Phi);
}


double ll_priorFunctionBinomial(NumericMatrix matDATA,NumericMatrix matPOP, NumericMatrix vMat, NumericMatrix Theta, NumericMatrix Phi, double alpha0,double alpha1, double gamma, bool ll_prior=true) {
  //'Initialize the logLikelihoodVec
  double logLikelihood=0;
  //'Total number of locations
  int nLocations=matDATA.nrow();
  //'Number of bands
  int nBands=matDATA.ncol();
  //Number of communities
  int n_community = Theta.ncol();
  //'Matrix of probabilities
  NumericMatrix tPhi = Rcpp::transpose(Phi);
  NumericMatrix probs=mmultBinomial(Theta,tPhi);
  //'Calculate the Loglikelihood and Prior

  double priorV=0.0;
  double priorPhi=0.0;

  if(ll_prior){
    //'Initialize the V_{cl} and Theta_{sc} prior
    for(int c=0;c<n_community;c++){
      for(int l=0;l<nLocations;l++){
        if(vMat(l,c)>0 && vMat(l,c)<1)  priorV=priorV+R::dbeta(vMat(l,c),1,gamma,true);
      }

      for (int b=0;b<nBands;b++){
        if(Phi(b,c)>0 && Phi(b,c)<1) priorPhi=priorPhi+R::dbeta(Phi(b,c),alpha0,alpha1,true);
      }
    }
  }

  for(int l=0;l<nLocations;l++){
    for (int b=0;b<nBands;b++){
      logLikelihood=logLikelihood+R::dbinom(matDATA(l,b),matPOP(l,b),probs(l,b),true);
    }
  }

  return(logLikelihood+priorV+priorPhi);
}


/***************************************************************************************************************************/
/*********************************            GIBBS SAMPLING PROCEDURE                  ************************************/
/***************************************************************************************************************************/

//' @name GibbsSamplingBinomial
//' @title Compute the Gibbs Sampling for LDA Binomial
//' @description Compute the Gibbs Sampling for LDA Binomial
//' @param DATA - DataFrame with Presence and Absecence (Binomial)
//' @param POP - DataFrame with Population Size (Binomial)
//' @param int n_community - Number of communities
//' @param alpha0 - Hyperparameter Beta(alpha0,alpha1)
//' @param alpha1 - Hyperparameter Beta(alpha0,alpha1)
//' @param gamma - Hyperparameter  Beta(1,gamma)
//' @param n_gibbs - Total number of Gibbs Samples
//' @param ll_prior - Likelihood compute with Priors ?
//' @param bool display_progress=true - Should I Show the progressBar ?
//' @return List - With Theta(n_gibbs,n_community*nSpecies), Phi(n_gibbs,nLocations*n_community) and logLikelihood
// [[Rcpp::export]]
List lda_binomial(DataFrame data,DataFrame pop, int n_community, double alpha0, double alpha1, double gamma, int n_gibbs, bool ll_prior=true, bool display_progress=true) {

  //'Convert to matrix
  NumericMatrix matDATA = internal::convert_using_rfunction(data, "as.matrix");

  //'Convert to matrix
  NumericMatrix matPOP = internal::convert_using_rfunction(pop, "as.matrix");

  //'Total number of locations
  int nLocations = matDATA.nrow();

  //'Total number of bands
  int nBands = matDATA.ncol();

  //'Intialize Theta
  NumericVector hyperTheta(n_community);
  hyperTheta.fill(1);
  NumericMatrix Theta=rdirichletBinomial(nLocations,hyperTheta);
  NumericMatrix vMat(nLocations,n_community);

  //'Intialize Phi
  NumericVector hyperPhi(n_community);
  hyperPhi.fill(1);
  NumericMatrix Phi=rdirichletBinomial(nBands,hyperPhi);

  //'Initialize the ThetaGibbs
  NumericMatrix ThetaGibbs(n_gibbs,nLocations*n_community);

  //'Initialize the PhiGibbs
  NumericMatrix PhiGibbs(n_gibbs,nBands*n_community);

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

    //'Initialize the zMAt
    GenericVector zMat = generateZBinomial(matDATA, matPOP, Theta, Phi);

    //Generate gamma MH
    if(bgamma){
      if ((g%50==0) & (g<500)){
        double z = acept/50;
        if ((z>0.4) & (jump<100))   jump=jump*2;
        if ((z<0.1) & (jump>0.001)) jump=jump*0.5;
        gamma = gammaMHBinomial(vMat, gamma, jump,acept);
      }
    }

    //'Generate Theta
    Theta = generateThetaBinomial(zMat,vMat, nLocations, n_community, gamma);

    //'Generate Phi
    Phi = generatePhiBinomial(zMat, nBands, n_community, alpha0, alpha1);

    //'Create the final Theta (n_gibbs,nLocations*n_community) and final Phi (n_gibbs,nBands*n_community)
    updateThetaAndPhiBinomial(ThetaGibbs, Theta, PhiGibbs, Rcpp::transpose(Phi), g);

    //'Initialize the logLikelihood
    double logLikelihood=ll_priorFunctionBinomial(matDATA,matPOP,
                                                  vMat, Theta, Phi,
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



//' @name GibbsSamplingBinomial
//' @title Compute the Gibbs Sampling for LDA Binomial
//' @description Compute the Gibbs Sampling for LDA Binomial
//' @param DATA - DataFrame with Presence and Absecence (Binomial)
//' @param POP - DataFrame with Population Size (Binomial)
//' @param int n_community - Number of communities
//' @param alpha0 - Hyperparameter Beta(alpha0,alpha1)
//' @param alpha1 - Hyperparameter Beta(alpha0,alpha1)
//' @param gamma - Hyperparameter  Beta(1,gamma)
//' @param n_gibbs - Total number of Gibbs Samples
//' @param n_burn - Number of elements to burn-in
//' @param ll_prior - Likelihood compute with Priors ?
//' @param bool display_progress=true - Should I Show the progressBar ?
//' @return List - With Theta(n_gibbs,n_community*nSpecies), Phi(n_gibbs,nLocations*n_community) and logLikelihood
List lda_binomial_burn(DataFrame data,DataFrame pop, int n_community, double alpha0, double alpha1, double gamma, int n_gibbs,int n_burn, bool ll_prior=true, bool display_progress=true) {

  //'Convert to matrix
  NumericMatrix matDATA = internal::convert_using_rfunction(data, "as.matrix");

  //'Convert to matrix
  NumericMatrix matPOP = internal::convert_using_rfunction(pop, "as.matrix");

  //'Total number of locations
  int nLocations = matDATA.nrow();

  //'Total number of bands
  int nBands = matDATA.ncol();

  //'Intialize Theta
  NumericVector hyperTheta(n_community);
  hyperTheta.fill(1);
  NumericMatrix Theta=rdirichletBinomial(nLocations,hyperTheta);
  NumericMatrix vMat(nLocations,n_community);

  //'Intialize Phi
  NumericVector hyperPhi(n_community);
  hyperPhi.fill(1);
  NumericMatrix Phi=rdirichletBinomial(nBands,hyperPhi);

  //'Initialize the ThetaGibbs
  NumericMatrix ThetaGibbs(n_gibbs-n_burn,nLocations*n_community);

  //'Initialize the PhiGibbs
  NumericMatrix PhiGibbs(n_gibbs-n_burn,nBands*n_community);

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

    //'Initialize the zMAt
    GenericVector zMat = generateZBinomial(matDATA, matPOP, Theta, Phi);

    //Generate gamma MH
    if(bgamma){
      if ((g%50==0) & (g<500)){
        double z = acept/50;
        if ((z>0.4) & (jump<100))   jump=jump*2;
        if ((z<0.1) & (jump>0.001)) jump=jump*0.5;
        gamma = gammaMHBinomial(vMat, gamma, jump,acept);
      }
    }

    //'Generate Theta
    Theta = generateThetaBinomial(zMat,vMat, nLocations, n_community, gamma);

    //'Generate Phi
    Phi = generatePhiBinomial(zMat, nBands, n_community, alpha0, alpha1);
    if(g>n_burn){
      //'Create the final Theta (n_gibbs,nLocations*n_community) and final Phi (n_gibbs,nBands*n_community)
      updateThetaAndPhiBinomial(ThetaGibbs, Theta, PhiGibbs, Rcpp::transpose(Phi), cont);

      //'Initialize the logLikelihood
      double logLikelihood=ll_priorFunctionBinomial(matDATA,matPOP,
                                                    vMat, Theta, Phi,
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
                                    Rcpp::Named("Phi")  = PhiGibbs,
                                    Rcpp::Named("logLikelihood")  =logLikelihoodVec);

  return resTemp;
}

