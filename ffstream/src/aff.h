#ifndef GUARD_aff_h
#define GUARD_aff_h

#include "fff.h"
//#include<Rcpp.h>
//#include<cmath>
//#include<vector>

//declare constants for AFF
const double LAMBDA_MIN = 0.6;
const double LAMBDA_MAX = 1.0;
const double INIT_ETA = 0.01;
const double INIT_DELTA = 0.0;
const double INIT_OMEGA = 0.0;
const double INIT_XBAR_DERIV = 0.0;
const double INIT_L_DERIV = 0.0;



class AFF : public FFF{
    public:
        AFF();
        AFF(double);
        double eta;
        void print();
        //overrides FFF::update
        void update(double);
        //a second update function that allows a scale factor to be passed
        void updateScaled(double, double);

        //just for testing
        double getLambda();
        double getXbar();
        double getS2();
        double getDelta();
        double getOmega();
        double getXbarDeriv();
        double getLderiv();
        void processVector(Rcpp::NumericVector);

        //used to check update of AFF
        Rcpp::List processVectorSave(Rcpp::NumericVector);
    private:
        double Omega;
        double Delta;
        double xbarDeriv;
        double Lderiv;
        void updateOmega();
        void updateDelta();
        void computeXbarDeriv();
        void computeLderiv(double);
        void updateLambda(double);
        void checkInBoundsLambda();
};


#endif
