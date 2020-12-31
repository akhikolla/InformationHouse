//// File Name: immer_rcpp_dnorm_pbivnorm_drezner.cpp
//// File Version: 1.08



// [[Rcpp::depends(RcppArmadillo)]]

// #include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;
// using namespace arma;

double const pi0 = 3.14159265359;


///********************************************************************
///** immer_pnorm
double immer_pnorm(double x)
{
    double val = ::Rf_pnorm5( x, 0.0, 1.0, 1, 0);
    return val;
}
///********************************************************************


///********************************************************************
///** immer_signum
double immer_signum(double x)
{
    double val=-1.0;
    if (x>0){
        val=1.0;
    }
    return val;
}
///********************************************************************

///********************************************************************
///** immer_dnorm2
double immer_dnorm2(double a, double b, double rho)
{
    double val = ( a*a + b*b - 2*rho*a*b ) / ( 2.0 * ( 1 - rho*rho ) );
    val = std::exp(-val) / ( 2.0 * pi0 * std::sqrt( 1 - rho*rho) );
    return val;
}
///********************************************************************

///********************************************************************
///** immer_dnorm
double immer_dnorm(double a)
{
    double val = std::exp( - a*a / 2.0 ) / ( std::sqrt( 2*pi0 ) );
    return val;
}
///********************************************************************

///********************************************************************
///** pbivnorm_drezner_fct_f_cpp
double pbivnorm_drezner_fct_f_cpp(double x, double y, double ap, double bp, double rho)
{
    double minval=-99;
    double res=0;
    if ( ( x > minval ) | ( y > minval ) ){
        res = ap*(2*x-ap) + bp*(2*y-bp) + 2*rho*(x-ap)*(y-bp);
        res = std::exp(res);
    } else {
        res = 0;
    }
    return res;
}
///********************************************************************

///********************************************************************
///** pbivnorm_drezner_all_neg
double pbivnorm_drezner_all_negative_cpp(double a, double b, double rho)
{
    Rcpp::NumericVector A(4);
    Rcpp::NumericVector B(4);
    A[0] = 0.3253030;
    B[0] = 0.1337764;
    A[1] = 0.4211071;
    B[1] = 0.6243247;
    A[2] = 0.1334425;
    B[2] = 1.3425378;
    A[3] = 0.006374323;
    B[3] = 2.2626645;
    double t1 = std::sqrt( 2*(1-rho*rho) );
    double ap = a / t1;
    double bp = b / t1;
    double val=0;
    double res=0;
    for (int ii=0; ii<4; ii++){
        for (int jj=0; jj<4; jj++){
            res = pbivnorm_drezner_fct_f_cpp(B[ii], B[jj], ap, bp, rho);
            val = val + A[ii]*A[jj]*res;
        }
    }
    val = std::sqrt( 1 - rho*rho) / pi0 * val;
    return val;
}
///********************************************************************

///********************************************************************
///** pbivnorm_drezner_product_negative_cpp
double pbivnorm_drezner_product_negative_cpp(double a, double b, double rho)
{
    double val=0;
    double val1=0;
    double a1=-a;
    double b1=-b;
    double rho1=-rho;
    ///------------
    if ( ( a<=0 ) & ( b>=0) & ( rho>=0) ){
        val1 = pbivnorm_drezner_all_negative_cpp(a, b1, rho1);
        val = immer_pnorm(a) - val1;
    }
    ///------------
    if ( ( a>=0 ) & ( b<=0) & ( rho>=0) ){
        val1 = pbivnorm_drezner_all_negative_cpp(a1, b, rho1);
        val = immer_pnorm(b) - val1;
    }
    ///------------
    if ( ( a>=0 ) & ( b>=0) & ( rho<=0) ){
        val1 = pbivnorm_drezner_all_negative_cpp(a1, b1, rho);
        val = immer_pnorm(a) + immer_pnorm(b) - 1 + val1;
    }
    ///------------
    if ( ( a<=0 ) & ( b<=0) & ( rho<=0) ){
        val = pbivnorm_drezner_all_negative_cpp(a, b, rho);
    }

    ///----- output
    return val;
}
///********************************************************************



///********************************************************************
///** pbivnorm_drezner_product_positive_cpp
double pbivnorm_drezner_product_positive_cpp(double a, double b, double rho)
{
    double t1 = std::sqrt( a*a - 2*a*b*rho + b*b);
    double sign_a = immer_signum(a);
    double sign_b = immer_signum(b);
    double rho1 = ( rho * a - b)*sign_a / t1;
    double rho2 = ( rho * b - a)*sign_b / t1    ;
    double delta = ( 1 - sign_a*sign_b )/4;
    double a1 = pbivnorm_drezner_product_negative_cpp(a, 0, rho1);
    double a2 = pbivnorm_drezner_product_negative_cpp(b, 0, rho2);
    double val = a1 + a2 - delta;
    ///----- output
    return val;
}
///********************************************************************

///********************************************************************
///** pbivnorm_drezner_numeric_arguments
double pbivnorm_drezner_numeric_arguments(double a, double b, double rho)
{
    double prod0 = a*b*rho;
    double val=0;
    if (prod0 > 0){
        val = pbivnorm_drezner_product_positive_cpp(a, b, rho);
    }
    if (prod0 <=0){
        val = pbivnorm_drezner_product_negative_cpp(a, b, rho);
    }
    return val;
}
///********************************************************************

///********************************************************************
///** pbivnorm_drezner
// [[Rcpp::export]]
Rcpp::NumericVector pbivnorm_drezner(Rcpp::NumericVector a, Rcpp::NumericVector b,
        Rcpp::NumericVector rho)
{
    int N=a.size();
    Rcpp::NumericVector res(N);
    for (int nn=0; nn<N; nn++){
        res[nn] = pbivnorm_drezner_product_positive_cpp( a[nn], b[nn], rho[nn] );
    }
    return res;
}
///********************************************************************


///********************************************************************
///** pbivnorm_drezner_derivative_rho_numeric
double pbivnorm_drezner_derivative_rho_numeric(double a, double b, double rho)
{
    double val = immer_dnorm2(a,b,rho);
    return val;
}
///********************************************************************

///********************************************************************
///** pbivnorm_drezner_derivative_a_numeric
double pbivnorm_drezner_derivative_a_numeric(double a, double b, double rho)
{
    double t1 = ( b - rho*a ) / std::sqrt( 1 - rho*rho );
    double val = immer_dnorm(a) * immer_pnorm(t1);
    return val;
}
///********************************************************************

///********************************************************************
///** pbivnorm_drezner_derivative_b_numeric
double pbivnorm_drezner_derivative_b_numeric(double a, double b, double rho)
{
    double t1 = ( a - rho*b ) / sqrt( 1 - rho*rho );
    double val = immer_dnorm(b) * immer_pnorm(t1);
    return val;
}
///********************************************************************

///********************************************************************
///** pbivnorm_drezner_derivative
// [[Rcpp::export]]
Rcpp::List pbivnorm_drezner_derivative(Rcpp::NumericVector a,
    Rcpp::NumericVector b, Rcpp::NumericVector rho)
{
    int N = a.size();
    Rcpp::NumericVector a_der(N);
    Rcpp::NumericVector b_der(N);
    Rcpp::NumericVector rho_der(N);
    for (int nn=0; nn<N; nn++){
        a_der[nn] = pbivnorm_drezner_derivative_a_numeric( a[nn], b[nn], rho[nn] );
        b_der[nn] = pbivnorm_drezner_derivative_b_numeric( a[nn], b[nn], rho[nn] );
        rho_der[nn] = pbivnorm_drezner_derivative_rho_numeric( a[nn], b[nn], rho[nn] );
    }
    return Rcpp::List::create(
                Rcpp::Named("a_der") = a_der,
                Rcpp::Named("b_der") = b_der,
                Rcpp::Named("rho_der") = rho_der
            );
}
///********************************************************************

///********************************************************************
///** immer_logdnorm2
double immer_logdnorm2( double x, double y, double mu1, double mu2, double var1,
        double var2, double cov12 )
{
    double val = -std::log(2*pi0);
    double det1 = var1*var2 - cov12*cov12;
    val += - 0.5*std::log(det1);
    double val1 = var2 * std::pow( x - mu1, 2.0) - 2 * cov12 * ( x- mu1)*(y-mu2) + var1*std::pow( y-mu2, 2.0);
    val += - 0.5*val1 / det1;
    return val;
}
///********************************************************************

///********************************************************************
///** immer_logdnorm2_derivative
// [[Rcpp::export]]
Rcpp::List immer_logdnorm2_derivative( double x, double y, double mu1, double mu2, double var1,
        double var2, double cov12 )
{
    // density value
    double val = -std::log(2*pi0);
    double det1 = var1*var2 - cov12*cov12;
    val += - 0.5*std::log(det1);
    double val1 = var2 * std::pow( x - mu1, 2.0) - 2 * cov12 * (x- mu1)*(y-mu2) + var1*std::pow( y-mu2, 2.0);
    val += - 0.5*val1 / det1;
    double logdens = val;
    double det1_squared = std::pow( det1, 2.0);

    // derivative with respect to mu1
    double der_mu1 = ( var2 * ( x - mu1 ) - cov12*(y-mu2) ) / det1;

    // derivative with respect to mu2
    double der_mu2 = ( var1 * ( y - mu2 ) - cov12*(x-mu1) ) / det1;

    // derivative with respect to var1
    double der_var1 = -0.5 / det1 * var2;
    der_var1 +=  0.5*val1 / det1_squared * var2;
    der_var1 += - 0.5 / det1 * std::pow( y - mu2, 2.0);

    // derivative with respect to var2
    double der_var2 = -0.5 / det1 * var1;
    der_var2 += 0.5*val1 / det1_squared * var1;
    der_var2 += -0.5 / det1 * std::pow( x - mu1, 2.0);

    // derivative with respect to cov12
    double der_cov12 = cov12 / det1;
    der_cov12 += - cov12*val1 / det1_squared;
    der_cov12 += ( x - mu1 )*( y-mu2 ) / det1;

    //--- output
    return Rcpp::List::create(
            Rcpp::Named("logdens") = logdens,
            Rcpp::Named("der_mu1") = der_mu1,
            Rcpp::Named("der_mu2") = der_mu2,
            Rcpp::Named("der_var1") = der_var1,
            Rcpp::Named("der_var2") = der_var2,
            Rcpp::Named("der_cov12") = der_cov12
        );
}
///********************************************************************
