/*****************************/
/* Jiaxing Lin               */
/* 05-10-2016                */
/* Main function to find the */
/* sigm for the maximum like-*/
/* lihood function for multi-*/
/* families.                 */
/*****************************/
//[[Rcpp::depends(BH)]]

#include <iostream>
#include <math.h>
#include "global.h"
#include <boost/bind.hpp>
#include "fminbr.h"
#include <ctime>
#include "fam_LLVar.h"
#include "Rcpp.h"
#include "famSize.h"

using boost::bind;
typedef boost::function<double(double x)> bindtype;

using namespace std;
using namespace Rcpp;

//[[Rcpp::export]]
List varEst( Rcpp::NumericVector fam_group,
             Rcpp::NumericVector alpha,
             Rcpp::NumericVector dtime,
             Rcpp::NumericVector delta,
             Rcpp::NumericVector g,       
             double beta,  // The value of beta
             double lower, // Lower bound of opt regime
             double upper  // Upper bound of opt regime
) 
{
  double* av = new double[alpha.size()];
  int*    dt = new int[dtime.size()];
  int* Delta = new int[delta.size()];
  int*     G = new int[g.size()];
  int      m = max(fam_group);
  double* logat = new double[alpha.size()];
  int*  famsize = new int[m];
  int* famgroup = new int[fam_group.size()]; 
 
  for(int i=0; i<alpha.size(); i++)
    av[i]=alpha[i];
  for(int i=0; i<alpha.size(); i++)
    logat[i]=log(av[i]);
  for(int i=0; i<dtime.size(); i++)
    dt[i]=dtime[i];
  for(int i=0; i<g.size(); i++)
    G[i]=g[i];
  for(int i=0; i<delta.size(); i++)
    Delta[i]=delta[i];  	 
  for(int i=0; i<fam_group.size();i++)
    famgroup[i] = fam_group[i];

  // generate number of subject for each family and put in famsize array.
  famSize(famsize, famgroup, fam_group.size());
   
  global_alpha_v_       = av;       
  global_Dtime_         = dt;   
  global_Delta_         = Delta;
  global_G_             = G;
  global_beta_          = &beta; 
  global_log_alpha_v_   = logat;
  
  bindtype foo = bind(fam_LLVar, _1, famsize, dt, Delta, G, m);
  /* 
   bind the function fam_LL to foo, so 
   that now foo is esentially the same
   as fam_LL, however, will prefilled
   arguments dt, Delta, G, and m. Only
   beta is the argument for function 
   foo.
   */
  double min = fminbr(lower,upper, foo, TOL);
  /*
   use the brent algorithm to search 
   to minimum value from function foo
   regarding to beta.
   */
  
  return List::create(Named("varEst") = (double)min);
  
  /*
   return the max likelihood and 
   corresponding values for beta
   */
}




