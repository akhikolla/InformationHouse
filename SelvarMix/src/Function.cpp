//
//  Function.cpp
//  
//
//  Created by sedki on 09/04/2014.
//
//
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <Rmath.h>
#include <Rdefines.h>

using namespace std;
using namespace Rcpp;
using namespace arma;

#include "Function.h"

long double Quad_Form(colvec x, colvec mu, mat S_Inv){
    colvec a(x.size());
    long double QuadForm;
    a = S_Inv * (x - mu);
    QuadForm = as_scalar(trans(x - mu) * a);
    return(QuadForm);
};

long double ldcppmvt(colvec x, colvec mu, mat SInv, double SLogDet){
    const long double log2pi = log(2.0 * M_PI);
    int xdim = x.size();
    long double constants = -0.5 * xdim * log2pi ;
    long double Qf =  Quad_Form(x, mu, SInv);
    long double lret = constants - (0.5 * SLogDet) - (0.5 * Qf);
    return(lret);
    
};
