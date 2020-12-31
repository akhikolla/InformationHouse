// [[Rcpp::depends(RcppEigen)]]
#include <Rcpp.h>
#include <RcppEigen.h>
#include "Inversion.h"
#include "RemlMM.h"
#include "RemlMM1Mat.h"
#include "RemlMMHen.h"
#include "RemlMMHenDiag.h"
#include "RemlMM2Mat.h"
#include "PrepMat.h"
#include "MLMM2Mat.h"
#include "MLMM.h"


using namespace Rcpp;
using namespace Eigen;
typedef Map<MatrixXd> Map_MatrixXd;



//[[Rcpp::export]]
List RemlMM(VectorXd Y, MatrixXd X , List VarList , VectorXd Init , int MaxIter , double CritVar , double CritLogLik) {
  List L = MM_RemlRcpp(Y, X, VarList, Init , MaxIter, CritVar, CritLogLik);
  return L;
}

//[[Rcpp::export]]
List RemlMM1Mat(VectorXd Y, MatrixXd X , MatrixXd VarInv, double logdetVar) {
  List L = MM_Reml1MatRcpp(Y, X, VarInv, logdetVar);
  return L;
}

// [[Rcpp::export]]
List sym_inverseRcpp(MatrixXd X) {
  MatrixXd x_inv(X);
  double log_detx,log_detxinv;
  sym_inverse(X,x_inv,log_detx,log_detxinv,0);
  return List::create(Rcpp::Named("inverse")=x_inv,Rcpp::Named("log_det")=log_detx);
}

//[[Rcpp::export]]
List RemlMMHen(VectorXd Y, MatrixXd X , MatrixXd Z , List GList , List GinvList , MatrixXd Rinv , VectorXd logdetV , VectorXd Init , int MaxIter , double CritVar , double CritLogLik) {
  List L = MM_RemlRcppHen(Y, X, Z , GList , GinvList , Rinv , logdetV , Init , MaxIter, CritVar, CritLogLik);
  return L;
}
//[[Rcpp::export]]
List RemlMMHenDiag(VectorXd Y, MatrixXd X , MatrixXd Z , List GList , List GinvList , VectorXd Rinv , VectorXd logdetV , VectorXd Init , int MaxIter , double CritVar , double CritLogLik) {
  List L = MM_RemlRcppHenDiag(Y, X, Z , GList , GinvList , Rinv , logdetV , Init , MaxIter, CritVar, CritLogLik);
  return L;
}

//[[Rcpp::export]]
List MM_Reml2MatRcpp(VectorXd Y, MatrixXd X , MatrixXd U , VectorXd D , VectorXd Init , int MaxIter , double CritVar , double CritLogLik) {
  List L = MM_Reml2Mat(Y, X, U , D, Init, MaxIter, CritVar, CritLogLik);
  return L;
}

//[[Rcpp::export]]
List PrepMat(VectorXd Y , MatrixXd K1 , MatrixXd K2){
  return PrepMatRcpp(Y,K1,K2);
}

//[[Rcpp::export]]
List MM_ML2MatRcpp(VectorXd Y, MatrixXd X , MatrixXd U , VectorXd D, VectorXd Init , int MaxIter , double CritVar , double CritLogLik) {
  List L = MM_ML2Mat(Y, X, U, D , Init , MaxIter, CritVar, CritLogLik);
  return L;
}

//[[Rcpp::export]]
List MLMM(VectorXd Y, MatrixXd X , List VarList , VectorXd Init , int MaxIter , double CritVar , double CritLogLik) {
  List L = MM_MLRcpp(Y, X, VarList, Init , MaxIter, CritVar, CritLogLik);
  return L;
}


template<typename T1, typename T2>
inline void chol_inverse(Eigen::MatrixBase<T1> & x, Eigen::MatrixBase<T2> & xi, double & log_det) {
  LDLT<MatrixXd> ldlt(x);
  
  log_det = ldlt.vectorD().array().log().sum();
  xi.setIdentity();
  ldlt.solveInPlace(xi);
}

// [[Rcpp::export]]
List chol_inverse(NumericMatrix X) {
  Map_MatrixXd x(as<Map<MatrixXd> >(X));
  double log_det;
  NumericMatrix Xi(X.rows(), X.cols());
  Map_MatrixXd xi(as<Map<MatrixXd> >(Xi));
  chol_inverse(x, xi, log_det);
  List L;
  L["inverse"] = Xi;
  L["log_det"] = log_det;
  return L;
}



// cette fonction badibulgue x en l'utilisant pour calculer SD 
// dans le cas oÃ¹ on n'a plus besoin de x une fois qu'on a son inverse Ã§a le fait
// Cette fonction n'utilise que le triangle supÃ©rieur de x...
// et ne remplit que le triangle supÃ©rieur de y !!
// avec eps = 0 calcule l'inverse
// avec eps petit... (1e-6 ?) pseudo inverse
void blocki(Eigen::MatrixXd & x, int x1, int n, Eigen::MatrixXd & y, int y1, double & log_det, double & det, double eps) {
  if(n == 1) {
    double d = (std::abs(x(x1,x1))<eps)?0:x(x1,x1);
    y(y1,y1) = (d==0)?0:1/d;
    det = d;
    log_det = log(d);
    return;
  }
  
  int m1 = n/2;
  int m2 = n-m1;
  
  Block<MatrixXd> A = x.block(x1,x1,m1,m1); 
  //Block<MatrixXd> D = x.block(x1+m1,x1+m1,m2,m2); 
  Block<MatrixXd> B = x.block(x1,x1+m1,m1,m2); 
  
  Block<MatrixXd> TL = y.block(y1,y1,m1,m1);
  Block<MatrixXd> BL = y.block(y1+m1,y1,m2,m1);
  Block<MatrixXd> TR = y.block(y1,y1+m1,m1,m2);
  Block<MatrixXd> BR = y.block(y1+m1,y1+m1,m2,m2);
  
  // BR = inverse(D)
  double log_detD, detD;
  blocki(x,x1+m1,m2,y,y1+m1,log_detD, detD, eps);
  
  // BL = inverse(D)*Bt
  BL.noalias() = BR.selfadjointView<Upper>() * B.transpose();
  
  // le bloc A est Ã©crasÃ© par SD
  A.triangularView<Upper>() -= B*BL; // on ne calcule que le triangle supÃ©rieur car on n'utilise pas l'autre
  
  // TL = inverse(SD)
  double log_detSD, detSD;
  blocki(x,x1,m1,y,y1,log_detSD, detSD, eps);
  
  TR.noalias() = TL.selfadjointView<Upper>()*(-BL.transpose());
  BR.triangularView<Upper>() -= BL*TR;  // on sait que cette matrice doit Ãªtre symmÃ©trique : on ne fait que la moitiÃ© des calculs
  
  log_det = log_detD + log_detSD;
  det = detD * detSD;
}

// la mÃªme en float
void blocki(Eigen::MatrixXf & x, int x1, int n, Eigen::MatrixXf & y, int y1, float & log_det, float & det, float eps) {
  if(n == 1) {
    float d = (std::abs(x(x1,x1))<eps)?0:x(x1,x1);
    y(y1,y1) = (d==0)?0:1/d;
    det = d;
    log_det = log(d);
    return;
  }
  
  int m1 = n/2;
  int m2 = n-m1;
  
  Block<MatrixXf> A = x.block(x1,x1,m1,m1); 
  //Block<MatrixXf> D = x.block(x1+m1,x1+m1,m2,m2); 
  Block<MatrixXf> B = x.block(x1,x1+m1,m1,m2); 
  
  Block<MatrixXf> TL = y.block(y1,y1,m1,m1);
  Block<MatrixXf> BL = y.block(y1+m1,y1,m2,m1);
  Block<MatrixXf> TR = y.block(y1,y1+m1,m1,m2);
  Block<MatrixXf> BR = y.block(y1+m1,y1+m1,m2,m2);
  
  // BR = inverse(D)
  float log_detD, detD;
  blocki(x,x1+m1,m2,y,y1+m1,log_detD, detD, eps);
  
  // BL = inverse(D)*Bt
  BL.noalias() = BR.selfadjointView<Upper>() * B.transpose();
  
  // le bloc A est Ã©crasÃ© par SD
  A.triangularView<Upper>() -= B*BL; // on ne calcule que le triangle supÃ©rieur car on n'utilise pas l'autre
  
  // TL = inverse(SD)
  float log_detSD, detSD;
  blocki(x,x1,m1,y,y1,log_detSD, detSD, eps);
  
  TR.noalias() = TL.selfadjointView<Upper>()*(-BL.transpose());
  BR.triangularView<Upper>() -= BL*TR;  // on sait que cette matrice doit Ãªtre symmÃ©trique : on ne fait que la moitiÃ© des calculs
  
  log_det = log_detD + log_detSD;
  det = detD * detSD;
}
