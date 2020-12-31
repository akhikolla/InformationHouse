// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <iostream>
#include <cmath>
using namespace Rcpp;
using namespace std;
using namespace arma;

// [[Rcpp::export]]
int signCPP(double x){
  if(x>0){
    return(1);
  }
  else if(x<0){
    return(-1);
  }
  else return(0);
}


// [[Rcpp::export]]
double SoftThreshCPP(double x,double lambda){
  double t=0;
  t=signCPP(x)*fmax(abs(x)-lambda,0);
  return(t);
}

// [[Rcpp::export]]
arma::Col<double> KernelCPP(arma::Row<double> x, arma::Mat<double> TrainData, double Sigma){
  int p = TrainData.n_cols;
  arma::Mat<double> Diff = TrainData.each_row()-x;
  Diff = square(Diff);
  arma::Col<double> RowSums = Diff*arma::ones(p);
  arma::Col<double> KernelVec = exp(-RowSums / Sigma);
  return(KernelVec);
}

// [[Rcpp::export]]
arma::Col<double> CoordDesCPP(arma::Col<double> w0, arma::Mat<double> Q, arma::Col<double> beta,double Lambda,double Epsilon, int Maxniter=1.0e7){
  int p = w0.size();
  int niter=0; // Declare iteration number.
  double error=1; // Declare error term. This will measure convergence of weigh
  arma::Col<double> w=w0; //Declare updated weight and initialize it to old weights.
  double normDiff=0;
  while((error>Epsilon)&&(niter<Maxniter)){
    normDiff=0;
    for(int i=0;i<p;i++){
      if (w0 [i] != 0){
        if(beta[i] == 0){
          w[i]=0;
        }
        else if(Q(i,i)==0){
          w[i]=signCPP(beta[i]);
         // std::cout<<i<<"-th Q Diagonal Term is zero without previous weight being zero." <<"\n";
        }
        else{
          double a = beta(i)-dot(Q.row(i),w)+Q(i,i)*w(i);
          double b = Lambda/2;
          w[i]=SoftThreshCPP(a,b);
          w[i]=w[i]/Q(i,i);
          if(w[i]>1){
            w[i]=1;
          }
          else if(-w[i]>1){
            w[i]=-1;
          }
        }

      }
      normDiff=normDiff+pow(w[i]-w0[i],2);
    }
    w0=w;
    niter++;
    error=sqrt(normDiff);
  }
  return(w);
}


// [[Rcpp::export]]
arma::Col<double> SolveKOSCPP(arma::Mat<double> YTheta, arma::Mat<double> K,double Gamma){
   //Assumes K is doubly-centered!
   int n = K.n_rows;
   //Solve for ridge solution. 
   arma::Mat<double> Mat = K+ n*Gamma*arma::diagmat(arma::ones(n)); //Add the n*Gamma*I term for ridge penalty
   arma::Col<double> Dvec = solve(Mat,YTheta);
   return(Dvec);
}

// [[Rcpp::export]]
arma::Mat<double> DerivCPP(arma::Row<double> x, arma::Mat<double> Data, arma::Col<double> w0, double sigmaD){
  arma::Mat<double> Diff = Data.each_row()-x;
  Diff = square(Diff);
  arma::Mat<double> DiffPart = (-2/sigmaD)*(Diff*diagmat(w0));
  int p = Data.n_cols;
  arma::Col<double> RowSums = Diff*arma::ones(p);
  arma::Col<double> KernPart = exp(-RowSums*(1/sigmaD));
  arma::Mat<double> result = arma::diagmat(KernPart)*DiffPart;
  return(result);
}

// [[Rcpp::export]]
arma::Mat<double> TMatCPP(arma::Mat<double> Data, arma::Col<double> Dvec, arma::Col<double> w0, double sigmaTm){
  int p = Data.n_cols;
  int n = Dvec.n_rows;
  arma::Mat<double> C = arma::eye(n,n)-(1/n)*arma::ones(n,n);
  Dvec = C*Dvec;
  arma::Mat<double> T = arma::zeros(n,p);
  for(int j=0; j<n; j++){
    T.row(j)=Dvec.t()*DerivCPP(Data.row(j), Data , w0, sigmaTm);
  }
  T = C*T;
  return(T);
}

// [[Rcpp::export]]
arma::Mat<double> compressedTMatCPP(arma::Mat<double> Data, arma::Mat<double> Q, arma::Col<double> compDvec, arma::Col<double> w0, double Sigma){
  int p = Data.n_cols;
  int m = compDvec.n_rows;
  arma::Mat<double> C = arma::eye(m,m)-(1/m)*arma::ones(m,m);
  compDvec = compDvec;
  
  arma::Mat<double> T = arma::zeros(m,p);
  for(int j=0; j<m; j++){
    T.row(j)=compDvec.t()*DerivCPP(Data.row(j), Data , w0, Sigma);
  }
  T = C*T;
  return(T);
}

// [[Rcpp::export]]
double ObjectiveFuncCPP(arma::Col<double> w, arma::Mat<double> Kw, arma::Mat<double> Data, arma::Col<double> DVectors, arma::Mat<double> YTheta, double LambdaOF, double GammaOF, double EpsilonOF=1e-5){
  int n = Kw.n_rows;
  double Regression = (1/n)*accu(square(YTheta-(Kw*DVectors))); //Least Squares component
  double LASSO = LambdaOF*sum(abs(w)); //LASSO Component
  arma::Mat<double> RidgeMat=DVectors.t()*(Kw*DVectors+EpsilonOF*DVectors);
  double Ridge = GammaOF*accu(arma::diagmat(RidgeMat)); //Ridge component
  return(Regression+LASSO+Ridge); //The sum of all three gives the objective function
}

// [[Rcpp::export]]
arma::Col<double> LambdaSeqCpp(double from, double to, double length){
  double Ratio = pow(to/from, 1/(length - 1) );
  arma::Col<double> LambdaSeq = arma::zeros(length);
  for(int i = 0; i < length; i++){
    LambdaSeq[i] = pow(Ratio, i) * from;
  }
  return(LambdaSeq);
}



// [[Rcpp::export]]
arma::Col<double> GetProjectionsCPP(arma::Mat<double> TrainData, arma::Col<int> TrainCat, arma::Mat<double> TestData, arma::Col<double> Dvec, arma::Col<double> w, arma::Mat<double> Kw, double Sigma, double Gamma){
  int n = TrainData.n_rows;
  int p = TrainData.n_cols;
  int m = TestData.n_rows;
  arma::Mat<double> Y = arma::zeros(n,2);
  
  //Fill Categorical Response
  int n1 = 0;
  int n2 = 0;
  for(int i = 0; i < n; i++){
    int Index = TrainCat[i]-1;
    Y(i, Index) = 1;
    if(TrainCat[i] == 1){
      n1++;
    }
    else{
      n2++;
    }
  }
  
  //Generate Transformed Response 
  arma::Col<double> theta = arma::zeros(2,1);
  theta[0] = pow(n2 / n1, 1/2);
  theta[1] = - pow(n1 / n2, 1/2);
  arma::Col<double> YTheta = Y*theta;
  
  // Weight TrainData and TestData
  for(int i = 0; i < p; i++){
    TrainData.col(i) = TrainData.col(i) * w[i];
    TestData.col(i) = TestData.col(i) * w[i];
  }
  
  //Colmeans Vector
  arma::Col<double> Means = Kw*arma::ones(n) / n;
  
  //Center Discriminant Vector
  double DvecMean = mean(Dvec);
  for(int i = 0; i < n; i++){
    Dvec[i] = Dvec[i] - DvecMean;
  }
  
  //Generate Projection Values
  arma::Col<double> ProjValues = arma::zeros(m);
  for(int i = 0; i < m; i++){
    arma::Row<double> x = TestData.row(i);
    arma::Col<double> KernVec = KernelCPP(x, TrainData, Sigma);
    ProjValues[i] = dot(KernVec - Means, Dvec);
  }
  return(ProjValues);
}
