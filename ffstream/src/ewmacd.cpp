#ifndef GUARD_ewmacd_cpp
#define GUARD_ewmacd_cpp

#include "ewmacd.h"
//#include "fff.h"
//#include<Rcpp.h>
//#include "utils.h"
//#include<vector>


//constructor
EwmaChangeDetector::EwmaChangeDetector(double r_, double L_, int BL_) :
    Detector(BL_),
    Z(INIT_EWMA_Z),
    r(r_),
    L(L_),
    sigmaZ(INIT_EWMA_SIGMA_Z),
    rFactorSigmaZ(INIT_EWMA_RFACTOR_SIGMA_Z){}


//compute the variance for the EWMA
void EwmaChangeDetector::computeSigmaZsq(){
    //first, sequentially update the rFactor
    //(1-r)^{2j} -> (1-r)^{2j+2}
    rFactorSigmaZ = rFactorSigmaZ * (1-r) * (1-r);

    //now compute sigmaZ
    //long term behaviour is below
    //   sigmaZ = sqrt(  r / (2-r)  ) * streamEstSigma;
    //
    //, but we use sqrt before multiplying by streamEstSigma:
    sigmaZ = sqrt(  r * (1 - rFactorSigmaZ) / (2-r)  ) * getStreamEstSigma();
}


//update Z according to:
//Z_j = (1-r) Z_{j-1} + r x_j
//Z_0 = mu
void EwmaChangeDetector::ewmaUpdate(double obs){
    Z = ( (1-r) * Z )  +  (r * obs);
    computeSigmaZsq();

    //computethe pval
    computePvalue();
}


//compute the p-value
//For EWMA, this will be a two-sided test (different to CUSUM)
void EwmaChangeDetector::computePvalue(){
    //first compute p-value according to location in (a, b)
    //a mapped to 0.025, b mapped to 0.975
    double a = getStreamEstMean() - L * sigmaZ;
    double b = getStreamEstMean() + L * sigmaZ;
    setPval( computeTwoSidedPvalue(Z, a, b) );
}


//stop the burn-in
void EwmaChangeDetector::stopBurnIn(){
    Detector::stopBurnIn();
    //For EWMA need to set Z to streamEstMean
    Z = getStreamEstMean();
    //For EWMA lso need to reset rFactorSigmaZ
    rFactorSigmaZ = 1;
}



//start the burn-in
void EwmaChangeDetector::startBurnIn(){
    Detector::startBurnIn();

    //reset the CUSUM and burn-in counter
    Z = INIT_EWMA_Z;
    sigmaZ = INIT_EWMA_SIGMA_Z;
    rFactorSigmaZ = INIT_EWMA_RFACTOR_SIGMA_Z;

}


//check if there is a change
void EwmaChangeDetector::checkIfChange(){
    double a = getStreamEstMean() - L * sigmaZ;
    double b = getStreamEstMean() + L * sigmaZ;
    if (Z < a)
        changeDetected = true;
    if (Z > b)
        changeDetected = true;
}


//main update for EWMA
void EwmaChangeDetector::update(double obs){

    //if a change occurred on the last update, need to reset to burnIn
    if (changeDetected){
        startBurnIn();
    }

    if (inBurnIn){
        streamEstimator.update(obs);
        BLcount++;
        //do update ewma during burn-in
        //   ewmaUpdate(obs); //actually, do not do the update during burn-in
        if ( !(BLcount < BL) ){
            stopBurnIn();
        }
    } else {
        ewmaUpdate(obs);
        checkIfChange();
    }

}



//getter for r
double EwmaChangeDetector::getR(){
    return r;
}


//getter for L
double EwmaChangeDetector::getL(){
    return L;
}


//getter for z
double EwmaChangeDetector::getZ(){
    return Z;
}


//getter for changeDetected
// bool EwmaChangeDetector::getChangeDetected(){
//     return(changeDetected);
// }


//print the EWMA values
void EwmaChangeDetector::print(){
    Rcpp::Rcout << "r: " << getR() << ", L: " << getL();
    Rcpp::Rcout << ", Z = " << getZ();
    Rcpp::Rcout << ", sigmaZ = " << sigmaZ;
    Rcpp::Rcout << ", rFactorSigmaZ= " << rFactorSigmaZ;
    Rcpp::Rcout << ", Burn in count: " << BLcount;
    Rcpp::Rcout <<  ", changeDetected: " << getChangeDetected();
    Rcpp::Rcout << std::endl;
}



//process a vector
Rcpp::List EwmaChangeDetector::processVectorSave (Rcpp::NumericVector vec){
    std::vector<bool> isChangeVec(vec.size());
    //should only be +1 because of trunc, but plus 2 for safety
    size_t maxNumChangepoints = (int) (vec.size() / BL) + 2;
    std::vector<int> changepoints(maxNumChangepoints);

    //keep a counter of the number of changepoints
    size_t numChangepoints = 0;

    for (int i=0; i < vec.size(); i++){
        update(vec[i]);
        isChangeVec[i] = changeDetected;
        if (changeDetected){
            numChangepoints++;
            //need to pushback i+1; index starts from 0 in cpp, but from 1 in R
            changepoints[numChangepoints-1] = (int) i+1;
        }
    }

    return Rcpp::List::create(Rcpp::Named(IS_CHANGE_DETECTED_FIELD_NAME)=Rcpp::wrap(isChangeVec),
                              Rcpp::Named(CHANGEPOINT_FIELD_NAME)=Rcpp::wrap(std::vector<int>(changepoints.begin(), changepoints.begin()+numChangepoints)));
}



//-----------------------------------------------------------------------//



RCPP_MODULE(ewmacdmodule){

    Rcpp::class_<Detector>( "Detector" )
        .constructor<int>()
        .property( "BL", &Detector::getBL, &Detector::setBL, "documentation for BL")
        .property( "pval", &Detector::getPval, "documentation for pval")
        .property("streamEstMean", &Detector::getStreamEstMean, &Detector::setStreamEstMean, "documentation for streamEstMean")
        .property("streamEstSigma", &Detector::getStreamEstSigma, &Detector::setStreamEstSigma, "documentation for streamEstSigma")
        .property("changeDetected", &Detector::getChangeDetected, "documentation for changeDetected")
        ;

    Rcpp::class_<EwmaChangeDetector>( "EwmaChangeDetector" )
        .derives<Detector>("Detector")
        .constructor<double, double, double>()
        .method( "update" , &EwmaChangeDetector::update, "documentation for update")
        .method( "print" , &EwmaChangeDetector::print, "documentation for print")
        .method( "processVectorSave" , &EwmaChangeDetector::processVectorSave, "documentation for processVectorSave")
//         .property("changeDetected", &EwmaChangeDetector::getChangeDetected, "documentation for changeDetected")
        .property("r", &EwmaChangeDetector::getR, "documentation for r")
        .property("L", &EwmaChangeDetector::getL, "documentation for L")
        .property("Z", &EwmaChangeDetector::getZ, "documentation for z")
        ;


}

#endif
