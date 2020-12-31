#include <math.h>
#include <ctime>

#include "s_hull_pro.h"

// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>

using namespace Rcpp;

using Eigen::MatrixXi;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::ArrayXd;
using Eigen::LLT;
using Eigen::Lower;
using Eigen::Map;
using Eigen::Upper;
using Eigen::HouseholderQR;
using Rcpp::as;
typedef Map<MatrixXd> MapMatd;
typedef Map<MatrixXi> MapMati;
typedef Map<VectorXd> MapVecd;
typedef Eigen::ColPivHouseholderQR<MatrixXd> CPivQR;
typedef CPivQR::PermutationType Permutation;

typedef struct triang{
  int nT;
  // indices of points
  std::vector<int> i1;
  std::vector<int> i2;
  std::vector<int> i3;
  // indices of neighbour triangles
  std::vector<int> j1;
  std::vector<int> j2;
  std::vector<int> j3;
  // circumcircle data
  std::vector<double> xc;
  std::vector<double> yc;
  std::vector<double> rc;
  // triangle area and ratio (ir/ccr)
  std::vector<double> ar;
  std::vector<double> rt;
  // convex hull
  std::vector<int> ch;
  int nch;
  // arcs, from to node indices
  std::vector<int> a1;
  std::vector<int> a2;
  // triangles to arcs indices
  std::vector<int> k1;
  std::vector<int> k2;
  std::vector<int> k3;
  int na;
} Triang;


typedef Eigen::Matrix< int , Eigen::Dynamic, 1> VectorXi;

typedef struct edges{
  int nE;
  VectorXi i1;
  VectorXi i2;
  VectorXi t1;
  VectorXi t2;
  MatrixXd xB;
  MatrixXd yB;
  MatrixXd zBl;
  MatrixXd zBr;
} Edges;

typedef struct nn{
  MatrixXi ind;
  MatrixXd dist;
} NN;

typedef struct cc{
  float xc;
  float yc;
  float rc;
  float ar;
} CC;

typedef struct pdest{
  VectorXd betahat;
  VectorXd est;
  VectorXd se;
  double cond;
} PDEst;




#define EIGEN_INITIALIZE_MATRICES_BY_NAN 1
#define EIGEN_USE_BLAS 1

MatrixXd AtA(MatrixXd A);
double threshold();
ArrayXd Dplus(const ArrayXd& d);
double kern2d(double x, double xi, double hx,
              double y, double yi, double hy,
              std::string kernel);

PDEst pD(NumericVector xD, NumericVector yD, NumericVector zD, NN nn,
         double x, double y, CharacterVector kernel, NumericVector h,
         std::string solver, int degree);
PDEst pDsmooth(NumericVector xD, NumericVector yD, NumericVector zD, NN nn,
           double x, double y, CharacterVector kernel, NumericVector h,
               std::string solver, int degree, int n, bool akimaweight);
triang shDt(std::vector<double> x, std::vector<double> y);
NN nN(NumericVector x, NumericVector y);
NN nN(VectorXd x, VectorXd y);
NN extendNN(NN nn, NumericVector X, NumericVector Y,
            NumericVector x, NumericVector y);
CC circum(double r1,double c1, double r2,double c2, double r3,double c3);
VectorXd myDnorm(VectorXd x, double mu, double sd);
