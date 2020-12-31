/*
 Jiaxing Lin
 Aug 25 2016
 */


/*
 main function to perform the constrict optmization using 
 gradient projection method.
 */

#include <iostream>
#include <unordered_set>
#include<fstream>


// include subroutines
#include "initialGuess.h"
#include "gradient.h"
#include "hessian.h"
#include "step1.h"
#include "step2.h"
#include "step3.h"
#include "step4.h"
#include "step5.h"
#include "sumterms.h"
#include "Ispline.h"
#include "Mspline.h"
#include <ctime>
#include <Rcpp.h>
#include <RcppEigen.h>

// global
#include "global.h"

//[[Rcpp::depends(RcppEigen)]]

using Eigen::VectorXd;
using Eigen::MatrixXd;
using namespace std;
using namespace Rcpp;

//[[Rcpp::export]]
NumericVector sieve(NumericVector Utmp, NumericVector Vtmp, NumericVector Markertmp, 
      NumericVector Deltatmp, NumericVector KnotI, NumericVector KnotM, int ki)
{

  VectorXd U(Utmp.size());
  for(int i = 0; i < Utmp.size(); i++)
    U(i) = Utmp(i);
  VectorXd V(Vtmp.size());
  for(int i = 0; i < Vtmp.size(); i++)
    V(i) = Vtmp(i);
  VectorXd Marker(Markertmp.size());
  for(int i = 0; i < Markertmp.size(); i++)
    Marker(i) = Markertmp(i);
  VectorXd Delta( Deltatmp.size());
  for(int i = 0; i < Deltatmp.size(); i++)
    Delta(i) = Deltatmp(i);
  VectorXd knot1(KnotI.size());
  for(int i = 0; i < KnotI.size(); i++)
    knot1(i) = KnotI(i);
  VectorXd knot2(KnotM.size());
  for(int i = 0; i < KnotM.size(); i++)
    knot2(i) = KnotM(i);
 
  //////////////////////////////////
  /* Declare and Alocate variables*/
  //////////////////////////////////
  
  //---------------------------------------------------------------------------
  // dimensions
  int pn = ki-1;
  int qn = ki-1;
  int dim = pn*qn+pn;
  
  // TOL of convergence
  // Yuan suggest 0.003 for epsilon 
  double Epsilon = 0.003;	
  
  // parameter vector
  VectorXd Theta(dim);
  VectorXd ThetaPre(dim);
  ThetaPre.fill(0.2/dim);
  
  // constrain condition matrix
  MatrixXd X(0,dim);
  
  // constrain set
  unordered_set<int> Alpha;
  
  // gradient vector
  VectorXd Gradient(dim);
  
  // hessian matrix
  MatrixXd Hessian(dim, dim);
  MatrixXd HessianPre(dim,dim);
  
  // fesible direction
  VectorXd Fdirection(dim);
  
  // initial guess of parameter Theta
  initialGuess(Theta);
  ThetaPre = Theta;	
  
  // compute spline function values
  MatrixXd ispline_U = Ispline(pn,U,knot1);
  MatrixXd ispline_V = Ispline(pn,V,knot1);
  MatrixXd mspline_m = Mspline(qn, Marker, knot2);
  
  bool converge = false;
  double Gamma  = 0;
  ////////////////////////
  /*Minimization process*/
  ////////////////////////

  int MaxIterCount = 200;
  int iterCount    = 0;
  
  std::clock_t start;
  double duration;

  
  int MemHold = -10;
  min_lamba = &MemHold;
  int MemHold2 = 0;
  cse = &MemHold2;
  
  
  // timing 
  clock_t start0 = clock(); 
  start = std::clock();
  	
  while(!converge && iterCount < MaxIterCount)
  {	
    iterCount++;
    
    //start = std::clock();
    Gradient = gradient(Theta, ispline_U, ispline_V, mspline_m, Delta, pn, qn);
    
    
    // step 0-2: compute the hessian matrix
    Hessian = hessian(Theta, ispline_U, ispline_V, mspline_m, Delta, pn, qn);
    
    // step 1: compute fesiable direction
    Fdirection = feaDirec(Hessian, Gradient, X);	
    
    
    // step 2: compute the gamma value needed in step 3
    Gamma = conGamma(Fdirection, Theta);
    
    
    // step 3: updates Theta.
    updateTheta(Fdirection, Gamma, Theta, ispline_U, ispline_V, mspline_m, Delta, pn, qn);
    
    // step 4: update X matrix.
    updateX(Theta, ThetaPre, X, Alpha);	
    ThetaPre = Theta;
    
    // step 5: checking convergence
    converge = checkConvergeAndUpdateX(Fdirection, Epsilon, Theta, X, Hessian, Gradient, Alpha);		 
  }
  
  NumericVector theta(Theta.size());
  
  if(converge)
  {
  	for(int i = 0; i < Theta.size(); i++)
    theta(i) = Theta(i);
  }
  else
  {	
    Rcpp::Rcout << "Converge Failed!" << endl; 
  	for(int i = 0; i < Theta.size(); i++)
	  theta(i) = 0;
  }
  return theta;
}
