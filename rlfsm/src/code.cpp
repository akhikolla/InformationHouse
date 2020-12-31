#include <Rcpp.h>
#include <iostream>
#include <vector>
#include <cmath>
using namespace std;
using namespace Rcpp;




double PP(int j, int m, double Ha) {

    double NN = pow(j, Ha);
    double DD = pow(m, Ha);
    double RR = NN/DD;
    return(RR);
}


/* This function doesnt work properly and isnt used
 // [[Rcpp::export]]
 double PLUS(int j, int m, double beta) {
 double F=0;
 if (j/m-1 > 0) {
 F=PP(j-m, m, beta);
 } else {

 }
 return(F);
 }
 */


// [[Rcpp::export]]
std::vector<double> a_tilda_cpp(int N, int m, int M, double alpha, double H){

    // it creates 0 in mM to m(M+N) anyway
    std::vector<double> a;
    a.assign(m*(M), 0);

    double X1;
    double X2 = pow(m,(-1/alpha));

    // Computation of a's
    for( int j = 1; j<m+1; j++) {
        //int j = 1;
        X1 = PP(j, m, H-1/alpha);

        a[j-1] = X1*X2;
    }

    for( int j = m+1; j<(m*M)+1; j++) {
        //int j = 1;
        X1 = PP(j, m, H-1/alpha);

        a[j-1] = (X1 - PP(j-m, m, H-1/alpha))*X2;
    }

    return a;
}


