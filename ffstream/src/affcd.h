#ifndef GUARD_affcd_h
#define GUARD_affcd_h

#include "detector.h"
#include "fff.h"
#include "aff.h"
#include "utils.h"
#include<Rcpp.h>
#include<vector>

//only needed for change detection version
// NOTE: Change from v0.1.5, now INIT_SCALEFACTOR is 1.
const double INIT_SCALEFACTOR = 1.0;
const double DEFAULT_SCALEFACTOR = 1.0;
const double ZERO_SCALEFACTOR = 0.0;
//Needs to be strictly larger than 1
const int BL_TOO_SMALL = 1;



//-----------------------------------------------------------------------//

//new class
//add BL
class AFFChangeDetector : public Detector{
    public:
        AFFChangeDetector(double); // alpha (eta is default value)
        AFFChangeDetector(double, double, int); //alpha, eta, BL

        double getAlpha();
        double getLambda();
        double getLderiv();

        void print();
        void update(double);
        void processVector(Rcpp::NumericVector);
        Rcpp::List processVectorSave(Rcpp::NumericVector);

        //this function is only here for testing
        void checkIfChange();
        
        //getters and setters for use with Rcpp
        double getAFFxbar();
        void setAFFxbar(double);
    private:
        AFF aff;
        double alpha;
        void stopBurnIn();
        double scaleFactor;  //this is 1/sigmahat, and scales eta
};

#endif
