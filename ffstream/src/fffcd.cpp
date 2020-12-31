#ifndef GUARD_fffcd_cpp
#define GUARD_fffcd_cpp


#include "fffcd.h"
#include "fff.h"
//#include "utils.h"
//#include<Rcpp.h>
//#include<vector>
// using std::vector;


//-----------------------------------------------------------------------//

//constructor - lambda, alpha
FFFChangeDetector::FFFChangeDetector(double lambda_, double alpha_) : 
        Detector(),
        fff(lambda_), 
        alpha(alpha_){}


        
//constructor - lambda, alpha, BL
FFFChangeDetector::FFFChangeDetector(double lambda_, double alpha_, int BL_) : 
        Detector(BL_),
        fff(lambda_), 
        alpha(alpha_){}


// print method
void FFFChangeDetector::print(){
    fff.printHeader(0);
    streamEstimator.printHeader(1);
    Rcpp::Rcout << "Change detector: " << std::endl;
    Rcpp::Rcout << "alpha: " << alpha << std::endl;
    Rcpp::Rcout << "Burn-in count: " << BLcount << std::endl;
    Rcpp::Rcout << "Burn-in length: " << getBL() << std::endl;
    Rcpp::Rcout << "change detected: " << changeDetected << std::endl;
}

//start a new burn-in
//need to redeclare here, because burn-in start process is slightly different
// void FFFChangeDetector::startBurnIn(){
//     Detector::startBurnIn();
// }


//check if a change has occurred
void FFFChangeDetector::checkIfChange(){
    //need to pass adjusted sigma
    double adjustedSigma = std::sqrt(fff.getU()) * getStreamEstSigma();
    double pval2sided = normcdf(fff.getXbar(), getStreamEstMean(), adjustedSigma);
    setPval( makeTwoSidedPvalueOneSided(pval2sided) );
    changeDetected = (getPval() < alpha);
}


void FFFChangeDetector::update(double obs){
    //first check if a change was detected after last obs...
    if (changeDetected){
        startBurnIn();
    }

    //if burn-in do...
    if(inBurnIn){
        //in burn-in mode
        fff.update(obs);
        streamEstimator.update(obs);
        BLcount++;
        if ( !( BLcount < getBL() ) ){ 
            //if ((BLcount < BL) == false){ 
            //if (!(BLcount < BL)){ 
            //if (BLcount >= BL){ 
            stopBurnIn();
            //changeMonitoringState();
            //add extra steps to initialise stream parameters
        }
    } else {
        //in Detect mode
        fff.update(obs);
        //needs to be last thing...
        checkIfChange();
    }
}



    

//-----------------------------------------------------------------------//
   
   
//just process a vector of observations
void FFFChangeDetector::processVector(Rcpp::NumericVector vec){
    for (int i=0; i < vec.size(); i++){
        update(vec[i]);
    }
}

 
//process vector and return changepoints
Rcpp::List FFFChangeDetector::processVectorSave(Rcpp::NumericVector vec){
    std::vector<bool> isChangeVec(vec.size());
    int maxNumChangepoints = (int) ( vec.size() / getBL() ) + 2;
    std::vector<int> changepoints(maxNumChangepoints);
    int numChangepoints = 0;
    
    for (int i=0; i < vec.size(); i++){
        update(vec[i]);
        isChangeVec[i] = changeDetected;
        if (changeDetected){
            //Rcpp::Rcout << i << ": change detected" << std::endl;
            numChangepoints++;
            changepoints[numChangepoints-1] = i+1;
        }            
    }

    return Rcpp::List::create(Rcpp::Named(IS_CHANGE_DETECTED_FIELD_NAME)=Rcpp::wrap(isChangeVec),
                              Rcpp::Named(CHANGEPOINT_FIELD_NAME)=Rcpp::wrap(std::vector<int>(changepoints.begin(), changepoints.begin()+numChangepoints)) ); 
}



//-----------------------------------------------------------------------//
//getters and setters for testing

double FFFChangeDetector::getAlpha(){
    return alpha;
}

void FFFChangeDetector::setAlpha(double alpha_){
    alpha = alpha_;
}

double FFFChangeDetector::getFFFxbar(){
    return fff.getXbar();
}

void FFFChangeDetector::setFFFxbar(double xbar_){
    fff.setXbar(xbar_);
}

// bool FFFChangeDetector::getChangeDetected(){
//     return changeDetected;
// }


//-----------------------------------------------------------------------//



RCPP_MODULE(fffcdmodule){

    Rcpp::class_<Detector>( "Detector" )
        .constructor<int>()
        .property( "BL", &Detector::getBL, &Detector::setBL, "documentation for BL")
        .property( "pval", &Detector::getPval, "documentation for pval")
        .property("streamEstMean", &Detector::getStreamEstMean, &Detector::setStreamEstMean, "documentation for streamEstMean")
        .property("streamEstSigma", &Detector::getStreamEstSigma, &Detector::setStreamEstSigma, "documentation for streamEstSigma")
        .property("changeDetected", &Detector::getChangeDetected, "documentation for changeDetected")
        ;

    Rcpp::class_<FFFChangeDetector>( "FFFChangeDetector" )
        .derives<Detector>("Detector")
        .constructor<double, double>()
        .constructor<double, double, int>()
        .property( "alpha", &FFFChangeDetector::getAlpha, &FFFChangeDetector::setAlpha, "documentation for alpha")
        .method( "print", &FFFChangeDetector::print, "documentation for print")
        .method( "update", &FFFChangeDetector::update, "documentation for update")
        .method( "processVector", &FFFChangeDetector::processVector, "documentation for processVector")
        .method( "processVectorSave", &FFFChangeDetector::processVectorSave, "documentation for processVectorSave")
        .property("fffxbar", &FFFChangeDetector::getFFFxbar, &FFFChangeDetector::setFFFxbar, "documentation for FFF xbar")
        .method( "checkIfChange", &FFFChangeDetector::checkIfChange, "documentation for checkIfChange")
        ;
}

#endif
