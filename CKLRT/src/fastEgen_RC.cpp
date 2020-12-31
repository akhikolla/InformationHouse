#include <RcppEigen.h>
#include <Rcpp.h>


// [[Rcpp::depends(RcppEigen)]]
using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::SelfAdjointEigenSolver;
typedef  Map<VectorXd>  MapVecd;
using namespace Rcpp;
//' Eigen_C
//' @param As A sysmetric matrix
//' @keywords internal
//[[Rcpp::export]]

List Eigen_C(NumericMatrix As){
  const Map<MatrixXd> A(as<Map<MatrixXd> >(As));
  SelfAdjointEigenSolver<MatrixXd> es(A);
  return List::create(Named("values") = es.eigenvalues().reverse(),
                      Named("vectors") = es.eigenvectors().rowwise().reverse());

}

//' Eigen_C_value
//' @param As A sysmetric matrix
//' @keywords internal
//[[Rcpp::export]]
Eigen::VectorXd Eigen_C_value(NumericMatrix As){
  const Map<MatrixXd> A(as<Map<MatrixXd> >(As));
  SelfAdjointEigenSolver<MatrixXd> es(A);
  return  es.eigenvalues().reverse();
}

//' MatMult_C
//' @param A first matrix
//' @param B second matrix
//' @keywords internal
// [[Rcpp::export]]
SEXP MatMult_C(Eigen::MatrixXd A, Eigen::MatrixXd B){
  Eigen::MatrixXd C = A * B;
  return Rcpp::wrap(C);
}


//' Sum_C
//' @param AA Vector
//' @keywords internal
// [[Rcpp::export]]
double Sum_C(NumericVector AA){
  const MapVecd A(as<MapVecd>(AA));
  double result= A.sum();
  return result;
}

//' ColSum_C
//' @param AA Matrix
//' @keywords internal
// [[Rcpp::export]]
NumericVector ColSum_C(NumericMatrix AA){
  const Map<MatrixXd> A(as<Map<MatrixXd> >(AA));
  NumericVector result;
  result = A.colwise().sum();
  return result;
}

//' MatrixRowMax_C
//' @param AA Matrix
//' @keywords internal
// [[Rcpp::export]]
NumericVector MatrixRowMax_C(NumericMatrix AA){
  const Map<MatrixXd> A(as<Map<MatrixXd> >(AA));
  NumericVector result;
  result = A.rowwise().maxCoeff();
  return result;
}



//' Elementwisesquare_C
//' @param AA Matrix
//' @keywords internal
// [[Rcpp::export]]
NumericVector Elementwisesquare_C(NumericMatrix AA){
  const Map<MatrixXd> A(as<Map<MatrixXd> >(AA));
  NumericVector result;
  result = A.array().square();
  return result;
}




//' VecMultMat_C
//' @param A Vector
//' @param B Matrix
//' @keywords internal
// [[Rcpp::export]]
NumericVector VecMultMat_C(Eigen::VectorXd A,Eigen::MatrixXd  B){
  Eigen::VectorXd C = A.transpose()*B;
  return Rcpp::wrap(C);
}

//' Vecplus_C
//' @param A Vector
//' @param B Vector
//' @keywords internal
// [[Rcpp::export]]
NumericVector Vecplus_C(Eigen::VectorXd A,Eigen::VectorXd B){
  Eigen::VectorXd C = A+B;
  return Rcpp::wrap(C);
}


//' ColSumtwomatrix_C
//' @param AA Matrix
//' @param BB Matrix
//' @keywords internal
// [[Rcpp::export]]
NumericVector ColSumtwomatrix_C(NumericMatrix AA,NumericMatrix BB){
  NumericVector result;
  result = ColSum_C(AA)+ColSum_C(BB);
  return result;
}

//' ifelsetest_C
//' @param x Vector
//' @keywords internal
//[[Rcpp::export]]
NumericVector ifelsetest_C( NumericVector x){
  return Rcpp::wrap( ifelse( x < 0, 0, x ));
}



//' MatrixPlus_C
//' @param A First Matrix
//' @param B Second Matrix
//' @keywords internal
// [[Rcpp::export]]
SEXP MatrixPlus_C(Eigen::MatrixXd A, Eigen::MatrixXd B){
  Eigen::MatrixXd C = A + B;
  return Rcpp::wrap(C);
}


//' NumxMatrix_C
//' @param A Number
//' @param B Matrix
//' @keywords internal
// [[Rcpp::export]]
SEXP NumxMatrix_C(double A, Eigen::MatrixXd B){
  Eigen::MatrixXd C = A * B;
  return Rcpp::wrap(C);
}

//' LR0_fixRho_C
//' @param LamdasR Lamda Number
//' @param muR mu vector
//' @param w1R w1 vector
//' @param w2R w2 vector
//' @param nminuspx n-px
//' @keywords internal
// [[Rcpp::export]]
NumericMatrix LR0_fixRho_C(NumericVector LamdasR,
                           NumericVector muR,
                           NumericMatrix w1R,
                           NumericMatrix w2R,
                           int nminuspx){
  const Map<MatrixXd> w1(as<Map<MatrixXd> >(w1R));
  const Map<MatrixXd> w2(as<Map<MatrixXd> >(w2R));
  const Map<VectorXd> mu(as<Map<VectorXd> >(muR));
  const Map<VectorXd> Lamdas(as<Map<VectorXd> >(LamdasR));
  int length_lamda = Lamdas.size();
  int N = w2.size();
  double lam;
  Eigen::VectorXd lammu_con;
  Eigen::VectorXd lammu_case;
  Eigen::VectorXd Dn;
  Eigen::VectorXd Nn;
  Eigen::VectorXd DnonNn;
  NumericVector temp;
  NumericMatrix result(N,length_lamda);


  for(int i=0;i<length_lamda;i++){
  lam = Lamdas[i];
  lammu_con = 1/(1+lam*mu.array());
  lammu_case = 1- lammu_con.array();
  Dn = (lammu_con).transpose()*w1+w2;

  Nn = lammu_case.transpose()*w1;
  temp = nminuspx*(1+Nn.array()/Dn.array()).log()-
    (1+(lam*mu).array()).log().sum();

  result(_,i) = ifelsetest_C(temp);
  }
  return Rcpp::wrap(result);

}



//' doubleloop
//' @param K1R K1 matrix
//' @param K2R K2 matrix
//' @param P0R P0 matrix
//' @param AR A matrix
//' @param U1R U1 vector
//' @param wR w matrix
//' @param LamdasR Lamdas vector
//' @param nminuspx n-px
//' @param all_rho the rho vector
//' @param LR0_allRhoR the matrix of likelihood ratio
//' @keywords internal
// [[Rcpp::export]]
NumericMatrix doubleloop(NumericMatrix K1R,
                            NumericMatrix K2R,
                            NumericMatrix P0R,
                            NumericMatrix AR,
                            NumericMatrix U1R,
                            NumericMatrix wR,
                            NumericVector LamdasR,
                            int nminuspx,
                            NumericVector all_rho,
                            NumericMatrix LR0_allRhoR){

  const Map<MatrixXd> K1(as<Map<MatrixXd> >(K1R));
  const Map<MatrixXd> K2(as<Map<MatrixXd> >(K2R));
  const Map<MatrixXd> P0(as<Map<MatrixXd> >(P0R));
  const Map<MatrixXd> A(as<Map<MatrixXd> >(AR));
  const Map<MatrixXd> w(as<Map<MatrixXd> >(wR));
  const Map<MatrixXd> U1(as<Map<MatrixXd> >(U1R));
  const Map<VectorXd> Lamdas(as<Map<VectorXd> >(LamdasR));
  Map<MatrixXd> LR0_allRho(as<Map<MatrixXd> >(LR0_allRhoR));
  double rho;
  Eigen::MatrixXd AKA1 = A.transpose()*K1*A;
  Eigen::MatrixXd AKA2 = A.transpose()*K2*A;
  Eigen::MatrixXd U1W = U1*w;
  for(int j =1;j<all_rho.size();j++){
    rho = all_rho(j);
    Eigen::MatrixXd K = rho*K1+(1-rho)*K2;
    int Knrow = K1.rows();
    SelfAdjointEigenSolver<MatrixXd> es(K);
    Eigen::VectorXd Kevalues = es.eigenvalues().reverse();
    Eigen::MatrixXd Kematrix = es.eigenvectors().rowwise().reverse();
    int k;
    for(k = 0;k<Knrow;k++){
      if((Kevalues(k))<1e-10){
        break;
      }
    }
    Eigen::VectorXd xi = Kevalues.head(k);
    Eigen::MatrixXd xis = xi.array().sqrt().matrix().asDiagonal();
    Eigen::MatrixXd ximat = Kematrix.leftCols(k);
    Eigen::MatrixXd phi=ximat*xis;
    Eigen::MatrixXd phiPphi = phi.transpose()*P0*phi;
    SelfAdjointEigenSolver<MatrixXd> ephiPphi(phiPphi);
    Eigen::VectorXd mu = ephiPphi.eigenvalues().reverse();

    double muximax= mu.maxCoeff();
    if(muximax < xi.maxCoeff()){
      muximax = xi.maxCoeff();
    }
    mu = mu/muximax;
    xi = xi/muximax;
    Eigen::MatrixXd AKA = rho*AKA1+(1-rho)*AKA2;

    SelfAdjointEigenSolver<MatrixXd> eAKA(AKA);
    Eigen::MatrixXd U2 = eAKA.eigenvectors().rowwise().reverse();
    Eigen::MatrixXd ww= U2.transpose()*U1W;
    int n = U2.rows();
    Eigen::MatrixXd ww2 = ww.array().square();
    Eigen::MatrixXd w1 = ww2.topRows(k);

    Eigen::MatrixXd w2 = ww2.bottomRows(n-k).colwise().sum();

    if(mu.size()<k){
      Eigen::VectorXd munew=xi;
      munew.head(mu.size())=mu;
      munew.tail(k-mu.size()).array() = 0 ;
      mu = munew;
    }

    int length_lamda = Lamdas.size();
    int N = w2.size();
    double lam;
    Eigen::VectorXd lammu_con;
    Eigen::VectorXd lammu_case;
    Eigen::VectorXd Dn;
    Eigen::VectorXd Nn;
    Eigen::VectorXd DnonNn;
    NumericVector temp;
    NumericMatrix LR0_fixRho(N,length_lamda);


    for(int i=0;i<length_lamda;i++){
      lam = Lamdas[i];
      lammu_con = 1/(1+lam*mu.array());
      lammu_case = 1- lammu_con.array();
      Dn = (lammu_con).transpose()*w1+w2;

      Nn = lammu_case.transpose()*w1;
      temp = nminuspx*(1+Nn.array()/Dn.array()).log()-
        (1+(lam*mu).array()).log().sum();

      LR0_fixRho(_,i) = ifelsetest_C(temp);
    }
    Map<MatrixXd> LR0_fixRhoE(as<Map<MatrixXd> >(LR0_fixRho));
    LR0_allRho.col(j) = LR0_fixRhoE.rowwise().maxCoeff();
  }


return Rcpp::wrap(LR0_allRho);

                            }


//' LR0_fixRho_LRT_C
//' @param LamdasR Lamda Number
//' @param muR mu vector
//' @param w1R w1 vector
//' @param w2R w2 vector
//' @param nminuspx n-px
//' @param xiR  Vector
//' @keywords internal
// [[Rcpp::export]]
NumericMatrix LR0_fixRho_LRT_C(NumericVector LamdasR,
                           NumericVector muR,
                           NumericMatrix w1R,
                           NumericMatrix w2R,
                           int nminuspx,
                           NumericVector xiR){
  const Map<MatrixXd> w1(as<Map<MatrixXd> >(w1R));
  const Map<MatrixXd> w2(as<Map<MatrixXd> >(w2R));
  const Map<VectorXd> mu(as<Map<VectorXd> >(muR));
  const Map<VectorXd> Lamdas(as<Map<VectorXd> >(LamdasR));
  const Map<VectorXd> xi(as<Map<VectorXd> >(xiR));
  int length_lamda = Lamdas.size();
  int N = w2.size();
  double lam;
  Eigen::VectorXd lammu_con;
  Eigen::VectorXd lammu_case;
  Eigen::VectorXd Dn;
  Eigen::VectorXd Nn;
  Eigen::VectorXd DnonNn;
  NumericVector temp;
  NumericMatrix result(N,length_lamda);


  for(int i=0;i<length_lamda;i++){
    lam = Lamdas[i];
    lammu_con = 1/(1+lam*mu.array());
    lammu_case = 1- lammu_con.array();
    Dn = (lammu_con).transpose()*w1+w2;

    Nn = lammu_case.transpose()*w1;
    temp = nminuspx*(1+Nn.array()/Dn.array()).log()-
      (1+(lam*xi).array()).log().sum();

    result(_,i) = ifelsetest_C(temp);
  }
  return Rcpp::wrap(result);

}

//' doubleloop_LRT
//' @param K1R K1 matrix
//' @param K2R K2 matrix
//' @param P0R P0 matrix
//' @param AR A matrix
//' @param U1R U1 vector
//' @param wR w matrix
//' @param LamdasR Lamdas vector
//' @param nminuspx n-px
//' @param all_rho rho vector
//' @param LR0_allRhoR LR0_allRhomatrix
//' @keywords internal
// [[Rcpp::export]]
NumericMatrix doubleloop_LRT(NumericMatrix K1R,
                         NumericMatrix K2R,
                         NumericMatrix P0R,
                         NumericMatrix AR,
                         NumericMatrix U1R,
                         NumericMatrix wR,
                         NumericVector LamdasR,
                         int nminuspx,
                         NumericVector all_rho,
                         NumericMatrix LR0_allRhoR){

  const Map<MatrixXd> K1(as<Map<MatrixXd> >(K1R));
  const Map<MatrixXd> K2(as<Map<MatrixXd> >(K2R));
  const Map<MatrixXd> P0(as<Map<MatrixXd> >(P0R));
  const Map<MatrixXd> A(as<Map<MatrixXd> >(AR));
  const Map<MatrixXd> w(as<Map<MatrixXd> >(wR));
  const Map<MatrixXd> U1(as<Map<MatrixXd> >(U1R));
  const Map<VectorXd> Lamdas(as<Map<VectorXd> >(LamdasR));
  Map<MatrixXd> LR0_allRho(as<Map<MatrixXd> >(LR0_allRhoR));
  double rho;
  Eigen::MatrixXd AKA1 = A.transpose()*K1*A;
  Eigen::MatrixXd AKA2 = A.transpose()*K2*A;
  Eigen::MatrixXd U1W = U1*w;
  for(int j =1;j<all_rho.size();j++){
    rho = all_rho(j);
    Eigen::MatrixXd K = rho*K1+(1-rho)*K2;
    int Knrow = K1.rows();
    SelfAdjointEigenSolver<MatrixXd> es(K);
    Eigen::VectorXd Kevalues = es.eigenvalues().reverse();
    Eigen::MatrixXd Kematrix = es.eigenvectors().rowwise().reverse();
    int k;
    for(k = 0;k<Knrow;k++){
      if((Kevalues(k))<1e-10){
        break;
      }
    }
    Eigen::VectorXd xi = Kevalues.head(k);
    Eigen::MatrixXd xis = xi.array().sqrt().matrix().asDiagonal();
    Eigen::MatrixXd ximat = Kematrix.leftCols(k);
    Eigen::MatrixXd phi=ximat*xis;
    Eigen::MatrixXd phiPphi = phi.transpose()*P0*phi;
    SelfAdjointEigenSolver<MatrixXd> ephiPphi(phiPphi);
    Eigen::VectorXd mu = ephiPphi.eigenvalues().reverse();

    double muximax= mu.maxCoeff();
    if(muximax < xi.maxCoeff()){
      muximax = xi.maxCoeff();
    }
    mu = mu/muximax;
    xi = xi/muximax;
    Eigen::MatrixXd AKA = rho*AKA1+(1-rho)*AKA2;


    SelfAdjointEigenSolver<MatrixXd> eAKA(AKA);
    Eigen::MatrixXd U2 = eAKA.eigenvectors().rowwise().reverse();
    Eigen::MatrixXd ww= U2.transpose()*U1W;
    int n = U2.rows();
    Eigen::MatrixXd ww2 = ww.array().square();
    Eigen::MatrixXd w1 = ww2.topRows(k);

    Eigen::MatrixXd w2 = ww2.bottomRows(n-k).colwise().sum();

    if(mu.size()<(k)){
      Eigen::VectorXd munew=xi;
      munew.head(mu.size())=mu;
      munew.tail(k-mu.size()).array() = 0 ;
      mu = munew;
    }

    int length_lamda = Lamdas.size();
    int N = w2.size();
    double lam;
    Eigen::VectorXd lammu_con;
    Eigen::VectorXd lammu_case;
    Eigen::VectorXd Dn;
    Eigen::VectorXd Nn;
    Eigen::VectorXd DnonNn;
    NumericVector temp;
    NumericMatrix LR0_fixRho(N,length_lamda);


    for(int i=0;i<length_lamda;i++){
      lam = Lamdas[i];
      lammu_con = 1/(1+lam*mu.array());
      lammu_case = 1- lammu_con.array();
      Dn = (lammu_con).transpose()*w1+w2;

      Nn = lammu_case.transpose()*w1;
      temp = nminuspx*(1+Nn.array()/Dn.array()).log()-
        (1+(lam*xi).array()).log().sum();

      LR0_fixRho(_,i) = ifelsetest_C(temp);
    }
    Map<MatrixXd> LR0_fixRhoE(as<Map<MatrixXd> >(LR0_fixRho));
    LR0_allRho.col(j) = LR0_fixRhoE.rowwise().maxCoeff();
 }


  return Rcpp::wrap(LR0_allRho);

}
