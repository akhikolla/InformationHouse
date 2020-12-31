//// File Name: mnlfa_rcpp_mnlfa.cpp
//// File Version: 0.32



// [[Rcpp::depends(RcppArmadillo)]]

// #include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;


//*** constants
const double pi1 = 3.14159265359;


///********************************************************************
///** mnlfa_rcpp_plogis
// [[Rcpp::export]]
double mnlfa_rcpp_plogis( double x)
{
    double z=0;
    if (x<0){
        z = std::exp(x);
        z = z / ( 1 + z );
    } else {
        z = std::exp(-x);
        z = 1 / ( 1 + z);
    }
    //--- OUTPUT
    return z;
}
///********************************************************************


///********************************************************************
///** mnlfa_rcpp_calc_probs_2pl
// [[Rcpp::export]]
Rcpp::NumericMatrix mnlfa_rcpp_calc_probs_2pl(
    Rcpp::NumericMatrix a, Rcpp::NumericMatrix b, Rcpp::NumericMatrix theta,
    Rcpp::IntegerVector y, Rcpp::LogicalVector y_resp )
{
    int N = a.nrow();
    int TP = theta.nrow();
    Rcpp::NumericMatrix like(N,TP);
    like.fill(1);
    double temp=0;
    for (int nn=0;nn<N; nn++){
        if (y_resp[nn] ){
            for (int tt=0; tt<TP; tt++){
                temp = mnlfa_rcpp_plogis( a(nn,0)*theta(tt,0) - b(nn,0) );
                if (y[nn] == 1){
                    like(nn,tt) = temp;
                } else {
                    like(nn,tt) = 1 - temp;
                }
            }
        }
    }
    //--- OUTPUT
    return like;
}
///********************************************************************


///********************************************************************
///** mnlfa_rcpp_mstep_trait_unidim
// [[Rcpp::export]]
double mnlfa_rcpp_mstep_trait_unidim( Rcpp::NumericMatrix theta, Rcpp::NumericMatrix mu_p,
    Rcpp::NumericMatrix sigma_p, Rcpp::NumericMatrix post)
{
    double theta_tt=0;
    double tmp=0;
    double tmp1=0;
    int N = mu_p.nrow();
    int TP = theta.nrow();
    double val=0;
    double const1 = -std::log(2*pi1)/2;
    double eps=1e-10;;
    for (int tt=0; tt<TP; tt++){
        theta_tt = theta(tt,0);
        for (int nn=0; nn<N; nn++){
            tmp = ( theta_tt - mu_p(nn,0) ) / ( sigma_p(nn,0)+eps);
            tmp1 = const1 - std::log(sigma_p(nn,0)+eps) - tmp*tmp/2;
            val += post(nn,tt)*tmp1;
        }
    }
    //--- OUTPUT
    return val;
}
///********************************************************************

