#ifndef GUARD_utils_cpp
#define GUARD_utils_cpp

#include "utils.h"
//#include<Rcpp.h>

//cmath not needed
// #include<cmath>


//adapted from code from:
//http://www.johndcook.com/blog/cpp_phi/
//"This code is in the public domain, Do whatever you want with it, no strings attached."
//accessed:23/03/2015
//Changes:
//have removed
//-the need to compute 1/sqrt(2)
//-the need to use fabs() - we were computing sign anyway
//the code depends on the formula given in Abramowitz & Stegun, Equation 7.1.26
//which gives |error| <= 1.5e-7
double stdnormcdf(double x){
    // constants
    double a1 =  0.254829592;
    double a2 = -0.284496736;
    double a3 =  1.421413741;
    double a4 = -1.453152027;
    double a5 =  1.061405429;
    double p  =  0.3275911;
    //double sqrtonehalf = 0.70710678118654757274;
    double sqrtonehalf = 0.70710678118;
 
    // Save the sign of x
    int sign = 1;
    if (x < 0)
        sign = -1;
    //x = fabs(x) * sqrtonehalf;
    //x = sign * x * sqrtonehalf;
    x *= (sign * sqrtonehalf);
 
    //Abramowitz & Stegun formula 7.1.26
    // A&S formula 7.1.26
    double t = 1.0/(1.0 + p*x);
    //y = erf x, using the formula
    double y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-x*x);
 
    return 0.5*(1.0 + sign*y);
}


//when not standard, then standardise
double normcdf(double x, double mu, double sigma){
    if (sigma==0){
        //throw warning? Maybe not
        sigma = 1;
    }
    return stdnormcdf( (x-mu)/sigma );
}





//if p = 0.95 then f(p) = 0.05
//if p = 0.05, then f(p) = 0.05
//if p = 0.975 then f(p) = 0.025
//This DOES NOT convert a two-sided p-value to a one-sided p-value
double convertPvalueToCorrectSide(double p){
    return ( 0.5 - std::abs(0.5 - p)  );
}


//convert two-sided p-value into one-sided pvalue
//actually, this might be correct. The function above might be incorrect
//if p = 0.95 then f(p) = 0.1
//if p = 0.05, then f(p) = 0.1
//if p = 0.975 then f(p) = 0.05
double makeTwoSidedPvalueOneSided(double pval2){
    return  1 - std::abs(1 - 2*pval2) ;
//     return convertPvalueToCorrectSide( 1 - std::abs(1 - 2*pval2) );
}


//compare x in the range (a,b), where a is mapped to 0, and b is mapped to 1.96
double computeOneSidedPvalue(double x, double a, double b){
    double q = 0;
    double p = DEFAULT_PVALUE;
    //check in case of pathological case that a >= b
    if (b > a){
        q = x / (b-a) * N01QUANTILE_0975;
        p = stdnormcdf(q);
        p = convertPvalueToCorrectSide(p);
    }
    return p;
}


//short function for max
// http://en.cppreference.com/w/cpp/algorithm/max
double getMax(double a, double b){
    if (a < b)
        return b;
    return a;
//     return (a < b) ? b : a;
}


//short function for max
// http://en.cppreference.com/w/cpp/algorithm/max
double getMin(double a, double b){
    if (a < b)
        return a;
    return b;
//     return (a < b) ? a : b;
}



double combineTwoOneSidedPvalues(double p1, double p2){
    //combine p1 and p2 into p3 using p3 = 2 * min(p1, p2)
//     return ( 2 * getMin(p1, p2) );
    if (p1 < p2)
        return 2 * p1;
    return 2 * p2;
}


//map a to 0.025, b to 0975
//Formula:
//(x-a) * (q2 - q1) / (b-a) + q1
//this maps [a,b] to [q1, q2]
//x: [a,b]
//x-a : [0, b-a]
//(x-a)/(b-a) : [0, 1]
//(x-a) * (q2-q1)/(b-a) : [0, q2-q1]
//(x-a) * (q2-q1)/(b-a) + q1 : [q1, q2]
//Formula seems to be correct
double computeTwoSidedPvalue(double x, double a, double b){
    double q = 0;
    double p = DEFAULT_PVALUE;
    //check in case of pathological case that a >= b
    if (b > a){
        q = (x-a) * (N01QUANTILE_0975 - N01QUANTILE_0025) / (b-a) + N01QUANTILE_0025;
        p = stdnormcdf(q);
        p = convertPvalueToCorrectSide(p);
    }
    return p;
}

#endif
