#ifndef GUARD_fff_h
#define GUARD_fff_h

#include "utils.h"
//#include<Rcpp.h>
//#include<vector>

// const double LAMBDA_NO_FORGETTING = 1.0;
// const double INIT_XBAR = 0.0;
// const double INIT_S2 = 0.0;
// const double INIT_M = 0.0;
// const double INIT_W = 0.0;
// const double INIT_U = 0.0;
// const double INIT_V = 0.0;

//-----------------------------------------------------------------------//

class FFF{
    public:
        //initialize lambda and alpha
        //need default contructor for AFF extension
        FFF();
        //second constructor with lambda specified
        FFF(double);
        void update(double);
        void print();
        void printHeader(int);
        void printAll();
        void reset();
        double getLambda();
        double getXbar();
        double getS2();
        double getW();
        double getU();
        double getV();
        void setXbar(double);
        void processVector(Rcpp::NumericVector);

    protected:
        //made m and w public
        //m is needed for Delta in AFF
       double lambda;
       double m;
       double w;
       double xbar;
       double s2;
    private:
       double u;
       double v;
       void computeXbar();
       void updateM(double);
       void updateW();
       void updateXbar(double);
       void updateU();
       void updateVFromUAndW();
       void updateS2(double);
};

#endif
