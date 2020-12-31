#ifndef GUARD_ewmacd_h
#define GUARD_ewmacd_h

#include "detector.h"
//#include "fff.h"
//#include<Rcpp.h>
//#include<vector>

// const double INIT_EWMA_Z = 0.0;
// const double INIT_EWMA_SIGMA_Z = 0.0;
// const double INIT_EWMA_RFACTOR_SIGMA_Z = 1.0;

//-----------------------------------------------------------------------//

class EwmaChangeDetector : public Detector{
    public:
        //r and L and BL
        EwmaChangeDetector(double, double, int);

        void update(double);
        void print();

        double getR();
        double getL();
        double getZ();
        Rcpp::List processVectorSave(Rcpp::NumericVector);
    private:
        double Z;
        double r;
        double L;
        double sigmaZ; //for variance of EWMA
        double rFactorSigmaZ; //the value (1-r)^{2j} - sequentially updated
        void ewmaUpdate(double);
        void stopBurnIn();
        void startBurnIn();
        void computeSigmaZsq();
        void checkIfChange();

        void computePvalue();
};

#endif
