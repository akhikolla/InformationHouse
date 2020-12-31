#ifndef GUARD_cusumcd_h
#define GUARD_cusumcd_h

#include "detector.h"
//#include "fff.h"
//#include "utils.h"
//#include<Rcpp.h>
//#include<vector>

// const double INIT_CUSUM_S = 0.0;
// const double INIT_CUSUM_T = 0.0;


//-----------------------------------------------------------------------//

class CusumChangeDetector : public Detector{
    public: 
        //k and h and BL
        CusumChangeDetector(double, double, int);

        void update(double);
        double getS();
        double getT();
        double getH();
        double getK();
        void print();
        Rcpp::List processVectorSave(Rcpp::NumericVector);

    private:
        double S;
        double T;
        double h;
        double k;
        void cusumUpdate(double);
        double standardiseObs(double, double, double);
        void startBurnIn();
        void checkIfChange();

        void computePvalue();
};

#endif
