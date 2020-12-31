#ifndef GUARD_affcd_cpp
#define GUARD_affcd_cpp

#include "affcd.h"
#include "fff.h"
#include "aff.h"
#include "utils.h"
//#include<Rcpp.h>
//#include<vector>
// using std::vector;


//-----------------------------------------------------------------------//

//constructors
AFFChangeDetector::AFFChangeDetector(double alpha_) : 
        Detector(),
        aff(), 
        alpha(alpha_), 
        scaleFactor(INIT_SCALEFACTOR){}

        
AFFChangeDetector::AFFChangeDetector(double alpha_, double eta_, int BL_) : 
        Detector(BL_),
        aff(eta_), 
        alpha(alpha_), 
        scaleFactor(INIT_SCALEFACTOR){}


// print method
void AFFChangeDetector::print(){
    aff.print();
    streamEstimator.printHeader(1);
    Rcpp::Rcout << "Change detector: " << std::endl;
    Rcpp::Rcout << "alpha: " << alpha << std::endl;
    Rcpp::Rcout << "Burn-in count: " << BLcount << std::endl;
    Rcpp::Rcout << "Burn-in length: " << getBL() << std::endl;
    Rcpp::Rcout << "change detected: " << getChangeDetected() << std::endl;
}


//stop the burn-in
void AFFChangeDetector::stopBurnIn(){
    Detector::stopBurnIn();

    // NOTE: Change from v0.1.5, now INIT_SCALEFACTOR is 1.
    // This handles the case when BL = 0, and the streamEstSigma has not been
    // set before starting the run, in which case the update will not work.
    scaleFactor = INIT_SCALEFACTOR;
    if (BL > BL_TOO_SMALL){
        scaleFactor = 1 / streamEstimator.getS2();
    } else {
        //so less than BL_MIN
        if (getStreamEstSigmaSq() > 0){
            //NOTE: Needs to be the square here
            scaleFactor = 1 / getStreamEstSigmaSq();
        } else {
            //case: BL=0, no streamEstSigma
            //INIT_SCALEFACTOR should take care of this case, but it doesn't 
            scaleFactor = DEFAULT_SCALEFACTOR;
        }
    }
}


// //start the burn-in
// void AFFChangeDetector::startBurnIn(){
//     Detector::startBurnIn();
// }

//if the user manually sets the streamEstSigma, it will be taken to mean
//that BL should be zero, and to stop the burn-in
// bool AFFChangeDetector::checkSigmaNotZero(){
//     if (getStreamEstSigma() > 0){
//         Detector::neverBurnIn();
//         double tempSigma = getStreamEstSigma();
//         double tempMean = getStreamEstMean();
//         stopBurnIn();
//         setStreamEstMean(tempMean);
//         setStreamEstSigma(tempSigma);
//         return true;
//     }
//     return false;
// }

//small function just to make code more readable
bool isZero(int x){
    return(x==0);
}

//update the AFF change detector
void AFFChangeDetector::update(double obs){
    //first check if a change was detected after last obs...
    if (changeDetected){
        startBurnIn();
    }

    //need to add this for prechange case
    if (isZero(BL)){
        stopBurnIn();
    }

    //if burn-in do...else
    if(inBurnIn){
        //in burn-in mode - scaleFactor argument is 0
        //  aff.update() will not use scaleFactor
        aff.updateScaled(obs, ZERO_SCALEFACTOR);
        streamEstimator.update(obs);
        BLcount++;

        if ( !(BLcount < BL) ){ 
            stopBurnIn();
        }
    } else {
        //in Detect mode
        aff.updateScaled(obs, scaleFactor);
        checkIfChange();
    }
}



//need to declare normcdf method above this method...
void AFFChangeDetector::checkIfChange(){
    //need to pass adjusted sigma
    double adjustedSigma = std::sqrt(aff.getU()) * getStreamEstSigma();

//     double pval1sided = 2*normcdf(-std::abs(aff.getXbar()), getStreamEstMean(), adjustedSigma);
//     setPval( pval1sided );
    double pval2sided = normcdf(aff.getXbar(), getStreamEstMean(), adjustedSigma);
    setPval( makeTwoSidedPvalueOneSided(pval2sided) );
    changeDetected = (getPval() < alpha);
}



//-----------------------------------------------------------------------//
   
 
//process vector and save changepoints and lambdao (for the moment)
Rcpp::List AFFChangeDetector::processVectorSave(Rcpp::NumericVector vec){
    std::vector<bool> isChangeVec(vec.size());
    int maxNumChangepoints = (int) (vec.size() / BL) + 2;
    std::vector<int> changepoints(maxNumChangepoints);
    int numChangepoints = 0;
    std::vector<double> lambdasave(vec.size());

    for (int i=0; i < vec.size(); i++){
        update(vec[i]);
        //also save the lambda
        lambdasave[i] = aff.getLambda();
        isChangeVec[i] = changeDetected;
        if (changeDetected){
            //Rcpp::Rcout << i << ": change detected" << std::endl;
            numChangepoints++;
            changepoints[numChangepoints-1] = i+1;
        }            
    }

    //return changepoints and lambda
    return Rcpp::List::create(Rcpp::Named(IS_CHANGE_DETECTED_FIELD_NAME)=Rcpp::wrap(isChangeVec),
                  Rcpp::Named(CHANGEPOINT_FIELD_NAME)=
                            Rcpp::wrap(std::vector<int>(changepoints.begin(), 
                                        changepoints.begin()+numChangepoints)),                       Rcpp::Named(AFF_FIELD_NAME) = Rcpp::wrap(lambdasave)); 

}


//just process a vector of observations
void AFFChangeDetector::processVector(Rcpp::NumericVector vec){
    for (int i=0; i < vec.size(); i++){
        update(vec[i]);
    }
}

//-----------------------------------------------------------------------//
//getters and setters for testing
double AFFChangeDetector::getAFFxbar(){
    return aff.getXbar();
}


void AFFChangeDetector::setAFFxbar(double xbar_){
    aff.setXbar(xbar_);
}

// bool AFFChangeDetector::getChangeDetected(){
//     return changeDetected;
// }

double AFFChangeDetector::getLambda(){
    return aff.getLambda();
}

double AFFChangeDetector::getAlpha(){
    return alpha;
}


double AFFChangeDetector::getLderiv(){
    return aff.getLderiv();
}


//-----------------------------------------------------------------------//


RCPP_MODULE(affcdmodule){

    Rcpp::class_<Detector>( "Detector" )
        .constructor<int>()
        .property( "BL", &Detector::getBL, &Detector::setBL, "documentation for BL")
        .property( "pval", &Detector::getPval, "documentation for pval")
        .property("streamEstMean", &Detector::getStreamEstMean, &Detector::setStreamEstMean, "documentation for streamEstMean")
        .property("streamEstSigma", &Detector::getStreamEstSigma, &Detector::setStreamEstSigma, "documentation for streamEstSigma")
        .property("changeDetected", &Detector::getChangeDetected, "documentation for changeDetected")
        ;
    
    Rcpp::class_<AFFChangeDetector>( "AFFChangeDetector" )
        .derives<Detector>("Detector")
        .constructor<double>()
        .constructor<double, double, int>()
        .property( "alpha", &AFFChangeDetector::getAlpha, "documentation for alpha")
        .method( "print", &AFFChangeDetector::print, "documentation for print")
        .method( "update", &AFFChangeDetector::update, "documentation for update")
        .method( "processVector", &AFFChangeDetector::processVector, "documentation for processVector")
        .method( "processVectorSave", &AFFChangeDetector::processVectorSave, "documentation for processVectorSave")
        .property("affxbar", &AFFChangeDetector::getAFFxbar, &AFFChangeDetector::setAFFxbar, "documentation for AFF xbar")
        .method( "checkIfChange", &AFFChangeDetector::checkIfChange, "documentation for checkIfChange")
        .property( "lambda", &AFFChangeDetector::getLambda, "documentation for lambda")
        .property( "Lderiv", &AFFChangeDetector::getLderiv, "documentation for Lderiv")
        ; 
}


#endif
