#ifndef GUARD_aff_cpp
#define GUARD_aff_cpp

//need to declare fff before aff
#include "fff.h"
//now can declare aff
#include "aff.h"
#include "utils.h"
//#include<Rcpp.h>
//#include<cmath>
//#include<vector>
//using std::vector;


//---------------------------------------------------------------//

//default constructor
//aff1 <- new(AFF)
AFF::AFF() : FFF(),
             eta(INIT_ETA),
             Omega(INIT_OMEGA),
             Delta(INIT_DELTA), 
             xbarDeriv(INIT_XBAR_DERIV),
             Lderiv(INIT_L_DERIV){}

//constructor when eta is specified
//aff2 <- new(AFF, 0.01)
AFF::AFF(double eta_) : FFF(),
                         eta(eta_),
                         Omega(INIT_OMEGA),
                         Delta(INIT_DELTA),
                         xbarDeriv(INIT_XBAR_DERIV),
                         Lderiv(INIT_L_DERIV){}


//print: show the fields
void AFF::print(){
    Rcpp::Rcout << "AFF contents: " << std::endl;
    Rcpp::Rcout << "lambda: " << lambda << std::endl;
    Rcpp::Rcout << "xbar: " << getXbar() << std::endl;
    Rcpp::Rcout << "eta: " << eta << std::endl;
}


//the simple update mechanism
void AFF::update(double obs){
    //Lderiv must be computed first, since it uses values from
    //xbar_old
    computeLderiv(obs);

    updateOmega();
    updateDelta();
    //FFF version of updating - is the same
    //FFF::update(obs);
    FFF::update(obs);

    //now update lambda
    computeXbarDeriv();
    updateLambda(1.0);
}


//the simple update mechanism
void AFF::updateScaled(double obs, double scaleFactor){
    //Lderiv must be computed first, since it uses values from
    //xbar_old
    computeLderiv(obs);

    updateOmega();
    updateDelta();
    //FFF version of updating - is the same
    //FFF::update(obs);
    FFF::update(obs);

    //now update lambda
    computeXbarDeriv();
    //no need to check value of scaleFactor here
    updateLambda(scaleFactor);
}



void AFF::updateOmega(){
    //here w is w_old
    Omega = lambda * Omega + w;
}


void AFF::updateDelta(){
    //here m is m_old
    Delta = lambda * Delta + m;
}


void AFF::computeXbarDeriv(){
    //here all indices are'new'
    xbarDeriv = (Delta - getXbar() * Omega) / w;
    //same as:
    //xbarDeriv = (Delta * w - m * Omega) / (w*w);
}


void AFF::computeLderiv(double obs){
    Lderiv = 2 * (getXbar() - obs) * xbarDeriv;
}


//---------------------------------------------------------------//
//These functions are only used in this class

//make sure lambda is in bounds
void AFF::checkInBoundsLambda(){
// better:?
//     lambda = (lambda < LAMBDA_MAX) ? LAMBDA_MIN : lambda;
//     lambda = (lambda > LAMBDA_MAX) ? LAMBDA_MAX : lambda;
    if (lambda < LAMBDA_MIN)
        lambda = LAMBDA_MIN;
    else if (lambda > LAMBDA_MAX)
        lambda = LAMBDA_MAX;
}


//---------------------------------------------------------------//
//scale factor is 1/sigmahat
void AFF::updateLambda(double scaleFactor){
    lambda = lambda - eta * scaleFactor * Lderiv;
    checkInBoundsLambda();
}


//---------------------------------------------------------------//
//getters just for testing

double AFF::getLambda(){
    return lambda;
}

double AFF::getXbar(){
    return xbar;
}

double AFF::getS2(){
    return s2;
}

double AFF::getOmega(){
    return Omega;
}

double AFF::getDelta(){
    return Delta;
}

double AFF::getXbarDeriv(){
    return xbarDeriv;
}

double AFF::getLderiv(){
    return Lderiv;
}

//---------------------------------------------------------------//

void AFF::processVector(Rcpp::NumericVector vec){
    for (int i=0; i < vec.size(); i++){
        update(vec[i]);
    }
}
 
Rcpp::List AFF::processVectorSave(Rcpp::NumericVector vec){
    std::vector<double> lambdasave(vec.size());
    for (int i=0; i < vec.size(); i++){
        update(vec[i]);
        lambdasave[i] = lambda;
    }
    return Rcpp::List::create(Rcpp::Named(AFF_FIELD_NAME) = Rcpp::wrap(lambdasave) ); 
}
//---------------------------------------------------------------//



RCPP_MODULE(affmodule){
    Rcpp::class_<AFF>( "AFF" )
        .constructor("default constructor")
        .constructor<double>("constructor when eta is specified")
        .field( "eta", &AFF::eta, "documentation for eta")
        .method( "print", &AFF::print, "documentation for print")
         .method( "update", &AFF::update, "documentation for update")
        .property( "lambda", &AFF::getLambda, "documentation for aff lambda")
        .property( "xbar", &AFF::getXbar, "documentation for aff xbar")
        .property( "s2", &AFF::getS2, "documentation for aff s2")
        .property( "Omega", &AFF::getOmega, "documentation for Omega")
        .property( "Delta", &AFF::getDelta, "documentation for Delta")
        .property( "xbarDeriv", &AFF::getXbarDeriv, "documentation for xbarDeriv")
        .property( "Lderiv", &AFF::getLderiv, "documentation for Lderiv")
        .method( "processVector", &AFF::processVector, "documentation for processVector")
        .method( "processVectorSave", &AFF::processVectorSave, "documentation for processVectorSave")
        ;
}


#endif
