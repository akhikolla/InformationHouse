#ifndef GUARD_fffcd_h
#define GUARD_fffcd_h

#include "detector.h"
//#include "fff.h"
//#include<Rcpp.h>
//#include<vector>


//-----------------------------------------------------------------------//

class FFFChangeDetector : public Detector{
    public:
        FFFChangeDetector(double, double); //lambda and alpha
        FFFChangeDetector(double, double, int); //lambda, alpha and BL

        double getAlpha();
        void setAlpha(double);

        void print();
        void update(double);
        void processVector(Rcpp::NumericVector);
        Rcpp::List processVectorSave(Rcpp::NumericVector);
        void checkIfChange();
        
        //getters and setters so can be called by Rcpp...for testing...
        double getFFFxbar();
        void setFFFxbar(double);

    private:
        FFF fff;
        double alpha;
};

#endif
