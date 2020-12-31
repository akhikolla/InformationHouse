#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "progress.hpp"
#include <iostream>
#include <ctime>
#include <fstream>
// [[Rcpp::depends(RcppProgress)]]
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;



/***************************************************************************************************************************/
/*********************************                      UTILS          *****************************************************/
/***************************************************************************************************************************/


arma::mat digammaMat(arma::mat x){
  arma::mat xdi(x.n_rows,x.n_cols);
  for(int r=0;r<(int)x.n_rows;r++){
    for(int c=0;c<(int)x.n_cols;c++){
      xdi(r,c)=R::digamma(x(r,c));
    }
  }
  return(xdi);
}

arma::mat lgammaMat(arma::mat x){
  arma::mat xdi(x.n_rows,x.n_cols);
  for(int r=0;r<(int)x.n_rows;r++){
    for(int c=0;c<(int)x.n_cols;c++){
      xdi(r,c)=R::lgamma1p(x(r,c));
    }
  }
  return(xdi);
}

arma::colvec digammaVec(arma::colvec x){
  arma::colvec xdi(x.n_rows);
  for(int r=0;r<(int)x.n_rows;r++){
    xdi(r)=R::digamma(x(r));
  }
  return(xdi);
}

arma::colvec rowSums(arma::mat X){
  int nRows = X.n_rows;
  arma::colvec out(nRows);
  for(int i = 0; i < nRows; i++){
    out(i) = arma::sum(X.row(i));
  }
  return(out);
}

arma::rowvec colSums(arma::mat X){
  int nCols = X.n_cols;
  arma::rowvec out(nCols);
  for(int i = 0; i < nCols; i++){
    out(i) = arma::sum(X.col(i));
  }
  return(out);
}

arma::mat vec2mat(arma::colvec x, int ncol, int nrow){
  arma::mat matTemp(nrow,ncol);
  for(int j=0;j<ncol;j++){
    for(int i=0;i<(int)x.n_elem;i++){
      matTemp(i,j)=x(i);
    }
  }
  return(matTemp);
}

arma::rowvec meltVariational(arma::mat mat){
  //'Initialize the NumericVector
  arma::rowvec vec(mat.n_rows*mat.n_cols);
  //'Initialize the position
  int pos=0;
  for(int col=0;col<(int)mat.n_cols;col++){
    for(int row=0;row<(int)mat.n_rows;row++){
      //'meltAbundance the matrix
      vec(pos)=mat(row,col);
      //'Increment the position
      pos=pos+1;
    }
  }
  return(vec);
}

void updateThetaVariational(arma::mat &ThetaGibbs, arma::mat Theta, arma::mat &PhiGibbs, arma::mat Phi, int gibbs){
  //'meltAbundance the Theta and Phi matrix
  ThetaGibbs.row(gibbs)=meltVariational(Theta);
  PhiGibbs.row(gibbs)=meltVariational(Phi);
}

/***************************************************************************************************************************/
/*********************************            GIBBS SAMPLING FUNCTIONS           *******************************************/
/***************************************************************************************************************************/



double generateELBO(arma::mat a,arma::mat b,arma::mat c,arma::mat d,arma::cube m1,arma::cube m0,arma::mat dat,int n_community,int n_species,int n_obs,int nLocation,double gamma, double a_phi,double b_phi) {
  //Utils
  arma::mat dig_a=digammaMat(a);
  arma::mat dig_b=digammaMat(b);
  arma::mat dig_c=digammaMat(c);
  arma::mat dig_d=digammaMat(d);
  arma::mat dig_soma_ab=digammaMat(a+b);
  arma::mat dig_soma_cd=digammaMat(c+d);
  arma::mat prop=dat/n_obs;
  arma::mat absent=n_obs-dat;

  //Part 1
  arma::mat res1(nLocation, n_species);
  arma::mat res0(nLocation, n_species);

    for (int i=0; i<nLocation;i++){
      for (int j=0; j<n_species;j++){


        arma::rowvec m1sub = m1.subcube( span(i),  span(j), span() );
        arma::colvec m1sub2 = arma::conv_to< arma::colvec >::from(m1sub);

        arma::rowvec m0sub =  m0.subcube( span(i),  span(j), span() );
        arma::colvec m0sub2 = arma::conv_to< arma::colvec >::from(m0sub);

        res1(i,j)=sum(m1sub2%(dig_c.col(j)-dig_soma_cd.col(j)));
        res0(i,j)=sum(m0sub2%(dig_d.col(j)-dig_soma_cd.col(j)));
      }
    }
  arma::mat sumMat = dat%res1+absent%res0;
  double p1=arma::sum(arma::sum(sumMat,0));
  //Part 2

  //get digamma calcs
  arma::mat pos=dig_a-dig_soma_ab;
  arma::mat neg=dig_b-dig_soma_ab;
  arma::mat res(nLocation,n_community);
  res.col(0)=pos.col(0);
  res.col(1)=pos.col(1)+neg.col(0);
  for (int i=2;i< n_community;i++){
    //stick-breaking structure
    arma::mat negTemp = neg.cols( 0, (i-1) );
    res.col(i)=pos.col(i)+rowSums(negTemp);
  }

  double p2=0;
  for (int i=0;i<n_community;i++){
    arma::mat p21 = m1.slice(i)%dat+m0.slice(i)%absent;
    arma::mat p22 = vec2mat(res.col(i),n_species, nLocation);
    arma::mat ptemp = p21%p22;
    p2=p2+arma::sum(arma::sum(ptemp,0));
  }
  //Part 3
  arma::mat tmp = (gamma-1)*(dig_b-dig_soma_ab);
  double p3=arma::sum(arma::sum(tmp,0));
  //Part 4
  double p41=(a_phi-1)*arma::sum(arma::sum(dig_c-dig_soma_cd,0));
  double p42=(b_phi-1)*arma::sum(arma::sum(dig_d-dig_soma_cd,0));
  double p4=p41+p42;
  //Part 5
  double p5=0;
  for (int i=0;i<n_community;i++){
    arma::mat p51=dat%m1.slice(i)%arma::log(m1.slice(i));
    arma::mat p52=absent%m0.slice(i)%arma::log(m0.slice(i));
    p5=p5+arma::sum(arma::sum(p51,0))+arma::sum(arma::sum(p52,0));
  }
  //Part 6
  arma::mat p61=lgammaMat(a+b)-lgammaMat(a)-lgammaMat(b);
  arma::mat p62=(a-1)%(dig_a-dig_soma_ab)+(b-1)%(dig_b-dig_soma_ab);
  double p6=arma::sum(arma::sum(p61,0))+arma::sum(arma::sum(p62,0));
  //Part 7
  arma::mat p71=lgammaMat(c+d)-lgammaMat(c)-lgammaMat(d);
  arma::mat p72=(c-1)%(dig_c-dig_soma_cd)+(d-1)%(dig_d-dig_soma_cd);
  double p7=arma::sum(arma::sum(p71,0))+arma::sum(arma::sum(p72,0));
  return(p1+p2+p3+p4-p5-p6-p7);
}


Rcpp::List generateCDmatrix(int n_community,int n_species,int n_obs, double a_phi, double b_phi, arma::cube m1, arma::cube m0, arma::mat dat) {
  arma::mat matC(n_community,n_species);
  arma::mat matD(n_community,n_species);
  for (int c=0;c<n_community;c++){
    arma::mat m1Temp = m1.slice(c);
    arma::mat m0Temp = m0.slice(c);
    arma::mat tmp = dat%m1Temp;
    matC.row(c)=colSums(tmp)+a_phi;
    arma::mat tmp0 = (n_obs-dat)%m0Temp;
    matD.row(c) = colSums(tmp0)+b_phi;
  }
  //'Store the results
  Rcpp::List resTemp = Rcpp::List::create(Rcpp::Named("matC") = matC,
                                          Rcpp::Named("matD") = matD);
  return(resTemp);
}




Rcpp::List generateABmatrix(int nLocations,int n_community,int n_species,int n_obs, double gamma, arma::cube m1, arma::cube m0, arma::mat dat) {
  arma::mat matA(nLocations,n_community);
  arma::mat matB(nLocations,n_community);
  for (int c=0;c< n_community;c++){
    arma::mat m1Temp = m1.slice(c);
    arma::mat m0Temp = m0.slice(c);
    arma::mat tmp = dat%m1Temp+(n_obs-dat)%m0Temp;
    matA.col(c)=rowSums(tmp)+1.0;
    //Initialize cumsum matrix
    arma::mat cs_m1(nLocations, n_species);
    cs_m1.fill(0.0);
    arma::mat cs_m0(nLocations, n_species);
    cs_m0.fill(0.0);

    //get cumsum of m1 and m0
    if (c==(n_community-1)) {
      m1Temp = m1.slice(c);
      m0Temp = m0.slice(c);
      cs_m1 = m1Temp;
      cs_m0 = m0Temp;
    }
    if (c!=(n_community-1)){
      cs_m1.fill(0.0);
      cs_m0.fill(0.0);
      for (int j=(c+1);j<n_community;j++){
        m1Temp = m1.slice(j);
        m0Temp = m0.slice(j);
        cs_m1=cs_m1+m1Temp;
        cs_m0=cs_m0+m0Temp;
      }
    }
    arma::mat tmp0 = (dat%cs_m1)+((n_obs-dat)%cs_m0);
    matB.col(c)=rowSums(tmp0)+gamma;
  }

  //'Store the results
  Rcpp::List resTemp = Rcpp::List::create(Rcpp::Named("matA") = matA,
                                          Rcpp::Named("matB") = matB);
  return(resTemp);
}



Rcpp::List generateM1M0matrix(int nLocations,int n_community,int n_species, arma::mat matA, arma::mat matB, arma::mat matC, arma::mat matD) {
  //' Calculate digamma
  arma::mat matAdi = digammaMat(matA);
  arma::mat matABdi = digammaMat(matA+matB);
  //' Calculate the difference
  arma::mat stb=matAdi-matABdi;
  //' Normalize last column
  stb.col(n_community-1).fill(0.0);
  arma::colvec res(stb.n_rows);
  res.fill(0.0);
    for (int c=0;c<(int)stb.n_cols-1;c++){
      arma::colvec sumAB = matA.col(c)+matB.col(c);
      res = res + digammaVec(matB.col(c)) - digammaVec(sumAB);
      stb.col(c+1)=stb.col(c+1)+res;
    }

  //' Data hyperparameters
  arma::mat pdat1=digammaMat(matC)-digammaMat(matC+matD);
  arma::mat pdat0=digammaMat(matD)-digammaMat(matC+matD);
  //' Initialize M1 and M0
  arma::cube matM0(nLocations,n_species,n_community);
  arma::cube matM1(nLocations,n_species, n_community);
  for (int i=0;i<nLocations;i++){
    arma::vec stbVec =  arma::conv_to< arma::vec >::from(stb.row(i));
    for (int j=0;j<n_species;j++){
      arma::vec pdat1Vec =  arma::conv_to< arma::vec >::from(pdat1.col(j));
      arma::vec pdat0Vec =  arma::conv_to< arma::vec >::from(pdat0.col(j));
      //'M1
      arma::vec m1vec = arma::exp(pdat1Vec+stbVec);
      double sum1 = arma::sum(m1vec);
      //'M0
      arma::vec m0vec = arma::exp(pdat0Vec+stbVec);
      double sum0 = arma::sum(arma::exp(pdat0Vec+stbVec));
      for(int k=0;k<n_community;k++){
        matM1(i,j,k)=m1vec(k)/sum1;
        matM0(i,j,k)=m0vec(k)/sum0;
      }
    }
  }

  //'Store the results
  Rcpp::List resTemp = Rcpp::List::create(Rcpp::Named("M1") = matM1,
                                    Rcpp::Named("M0")  = matM0);

  return(resTemp);
}

arma::mat generateThetaBinomialVariational(int nLocations,int n_community, arma::mat matA, arma::mat matB) {
  //'Intialize the Theta matrix
  arma::mat Theta(nLocations,n_community);
  //'Create the mean
  arma::mat meanMatDenominator=matA+matB;
  arma::mat meanMat=matA/meanMatDenominator;
  //'Fix in one the last community
  meanMat.col(n_community-1).fill(1.0);
  //' Define the first two communities
  Theta.col(0)=meanMat.col(0);
  Theta.col(1)=meanMat.col(1)%(1.0-meanMat.col(0));
  for(int c=2;c<n_community;c++){
    arma::mat tmpMat = 1.0-meanMat.cols(0,c-1);
    arma::colvec prod = arma::prod(tmpMat, 1);
    Theta.col(c)=meanMat.col(c)%prod;
  }
  return(Theta);
}

/***************************************************************************************************************************/
/*********************************            VARIATIONAL INFERENCE              *******************************************/
/***************************************************************************************************************************/


// [[Rcpp::export]]
Rcpp::List lda_binomial_var(arma::mat data, int n_community, int maxit, int n_obs, double gamma, double a_phi, double b_phi, double thresh, double delta_elbo, arma::cube m1, arma::cube m0) {

 //Elbo Vector
 Rcpp::NumericVector elboVec(maxit);
 //Delta Elbo
 delta_elbo =  std::numeric_limits<double>::max();
 //Number of locations
 int nLocations = data.n_rows;
 //Number of locations
 int n_species = data.n_cols;
 int i=0;

 //'Initialize the ThetaGibbs
 arma::mat ThetaGibbs(maxit,nLocations*n_community);
 arma::mat PhiGibbs(maxit,n_community*n_species);


  while ((i < maxit) & (delta_elbo>thresh)){

    //Get A and B matrices
    Rcpp::List tmp0=generateABmatrix(nLocations, n_community, n_species, n_obs, gamma, m1, m0, data);
    arma::mat matA = tmp0(0);
    arma::mat matB = tmp0(1);

    Rcpp::List tmp1 = generateCDmatrix(n_community, n_species, n_obs, a_phi, b_phi, m1, m0, data);
    arma::mat matC = tmp1(0);
    arma::mat matD = tmp1(1);

    Rcpp::List tmp2 = generateM1M0matrix(nLocations, n_community, n_species, matA, matB, matC, matD);
    arma::cube tm1=tmp2(0);
    m1=tm1;
    arma::cube tm0=tmp2(1);
    m0=tm0;

    arma::mat Theta = generateThetaBinomialVariational(nLocations, n_community, matA, matB);
    arma::mat Phi = matC/(matC+matD);
    updateThetaVariational(ThetaGibbs, Theta, PhiGibbs, Phi, i);

    elboVec(i) = generateELBO(matA, matB, matC, matD, m1, m0, data, n_community, n_species, n_obs, nLocations, gamma, a_phi, b_phi);
    if (i!=0) delta_elbo=std::abs(elboVec(i)-elboVec(i-1));
    i=i+1;
    if(i>maxit) break;
  }

  //'Store the results
  Rcpp::List resTemp = Rcpp::List::create(Rcpp::Named("Theta") = ThetaGibbs,
                                          Rcpp::Named("Phi") = PhiGibbs,
                                          Rcpp::Named("Elbo")  = elboVec);

  return resTemp;
}


