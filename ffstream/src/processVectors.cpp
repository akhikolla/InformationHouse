#ifndef GUARD_processvectors_cpp
#define GUARD_processvectors_cpp

#include "processVectors.h"
//#include<Rcpp.h>




//-----------------------------------------------------------------//
//     Example initialisation WITH argument.
//     FFF fff(lambda);
//
//     Example initialisation with no argument.
//     FFF fff2;
//     fff2.print();
//-----------------------------------------------------------------//



//' Compute the FFF mean of a vector
//'
//' Given a vector \code{x} and a value \code{lambda} for a fixed forgetting
//' factor, returns the value of the fixed forgetting factor mean
//' \eqn{\bar{x}_{N, \lambda}}, where \eqn{N} is the length of \code{x}.
//' Algorithm is implemented in C++.
//'
//' @param x Vector of numeric values values.
//'
//' @param lambda Value for the fixed forgetting factor in \eqn{[0,1]}.
//' 
//' @keywords internal
// [[Rcpp::export]]
double cpp_computeFFFMean(Rcpp::NumericVector x, double lambda){
    //initialise with lambda
    FFF fff(lambda);

    //using NumericVector iterator
    for(Rcpp::NumericVector::iterator it = x.begin(); it != x.end(); ++it) {
        fff.update(*it);
    }

    //return Xbar
    return fff.getXbar();
}


//' Compute the AFF mean of a vector
//'
//' Given a vector \code{x} and a value \code{eta} for step size
//' in the stochastic gradient descent for the adaptive forgetting
//' factor, this returns the value of the fixed forgetting factor mean
//' \eqn{\bar{x}_{N, \overrightarrow{\lambda} }}, where \eqn{N} is the 
//' length of \code{x}. Algorithm is implemented in C++.
//'
//' @param x Vector of numeric values values.
//'
//' @param eta Value for the step size in the gradient descent step.
//' 
//' @keywords internal
// [[Rcpp::export]]
double cpp_computeAFFMean(Rcpp::NumericVector x, double eta){
    //initialise with lambda
    AFF aff(eta);

    //using NumericVector iterator
    for(Rcpp::NumericVector::iterator it = x.begin(); it != x.end(); ++it) {
        aff.update(*it);
    }

    //return Xbar
    return aff.getXbar();
}




#endif
