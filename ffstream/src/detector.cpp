#ifndef GUARD_detector_cpp
#define GUARD_detector_cpp

#include "detector.h"
//#include "utils.h"


Detector::Detector() :
    BL(DEFAULT_BL),
    BLcount(INIT_BL_COUNT),
    pval(DEFAULT_PVALUE),
    inBurnIn(INIT_IN_BURN_IN), 
    inDetectState(INIT_IN_DETECT_STATE), 
    changeDetected(INIT_CHANGE_DETECTED),
    streamEstimator(LAMBDA_NO_FORGETTING),
    streamEstMean(INIT_STREAM_EST_MEAN),
    streamEstSigma(INIT_STREAM_EST_SIGMA),
    streamEstSigmaSq(INIT_STREAM_EST_SIGMA_SQ){}


Detector::Detector(int BL_) :
    BL(BL_),
    BLcount(INIT_BL_COUNT),
    pval(DEFAULT_PVALUE),
    inBurnIn(INIT_IN_BURN_IN), 
    inDetectState(INIT_IN_DETECT_STATE), 
    changeDetected(INIT_CHANGE_DETECTED),
    streamEstimator(LAMBDA_NO_FORGETTING),
    streamEstMean(INIT_STREAM_EST_MEAN),
    streamEstSigma(INIT_STREAM_EST_SIGMA),
    streamEstSigmaSq(INIT_STREAM_EST_SIGMA_SQ){}


//destructor

//getter
int Detector::getBL(){
    return(BL);
}

//setter
void Detector::setBL(int BL_){
    BL = BL_;
}


double Detector::getStreamEstMean(){
    return streamEstMean;
}

//set the streamEstMean from the streamEstimator's xbar
void Detector::setStreamEstMean(){
    streamEstMean = streamEstimator.getXbar();
}

void Detector::setStreamEstMean(double streamEstMean_){
    streamEstMean = streamEstMean_;
}


double Detector::getStreamEstSigma(){
    return streamEstSigma;
}


void Detector::setStreamEstSigma(double streamEstSigma_){
    streamEstSigma = streamEstSigma_;
    //need to set square:
    streamEstSigmaSq = streamEstSigma * streamEstSigma;
}

//
//set the streamEstSigma from the streamEstimator's s2
double Detector::getStreamEstSigmaSq(){
    return streamEstSigmaSq;
}



//set the streamEstSigma from the streamEstimator's s2
// void Detector::setStreamEstSigma(){
//     streamEstSigma = std::sqrt( streamEstimator.getS2() );
// }


//if set by value, but actually this will never occur
void Detector::setStreamEstSigmaSq(double streamEstSigmaSq_){
    streamEstSigmaSq = streamEstSigmaSq_;
    //need to set square root:
    streamEstSigma = std::sqrt(streamEstSigmaSq);
}


// This will be called after streamEstimator is used
void Detector::setStreamEstSigmaSq(){
    streamEstSigmaSq = streamEstimator.getS2();
    //need to set square root:
    streamEstSigma = std::sqrt(streamEstSigmaSq);
}


//stop the burn in (end of burn in), so initialise streamEstMean and streamEstSigma
void Detector::stopBurnIn(){
    inBurnIn = false;
    //set stream's estimated mean and sigma based on xbar and s2 of streamEstimator
    //NB: only do this if BL unequal to zero
    if (BL != 0){
        Detector::setStreamEstMean();
        //note change here from Sigma to SigmaSq
        Detector::setStreamEstSigmaSq();
    }
    inDetectState = not inBurnIn;
    changeDetected = false;
}

//start a new burn-in
void Detector::startBurnIn(){
    inBurnIn = true;
    inDetectState = not inBurnIn;
    BLcount = INIT_BL_COUNT;
    streamEstimator.reset();
    changeDetected = false;
    pval = DEFAULT_PVALUE;
}


//getter
double Detector::getPval(){
    return(pval);
}

//setter
void Detector::setPval(double pval_){
    pval = pval_;
}


//getter
bool Detector::getChangeDetected(){
    return(changeDetected);
}



//detect multiple changes
Rcpp::List Detector::detectMultiple(Rcpp::NumericVector x){
    int maxNumChangepoints = (int) ( x.size() / BL ) + 2;
    std::vector<int> changepoints(maxNumChangepoints);
    int numChangepoints = 0;
    //counting from 1, because R indexes from 1.
    int index = 1;
    
    //run through the data points
    for(Rcpp::NumericVector::iterator it = x.begin(); it != x.end(); ++it) {
        update(*it);
        if (getChangeDetected()){
            changepoints[numChangepoints] = index;
            numChangepoints++;
        }            
        index++;
    }
    //an attempt to make the extraction of the changepoints more readable
    return Rcpp::List::create( 
            Rcpp::Named(CHANGEPOINT_FIELD_NAME)=
            Rcpp::wrap(std::vector<int>(changepoints.begin(), 
                                        changepoints.begin()+numChangepoints)
                      ) 
            ); 
}


//detect single change
Rcpp::List Detector::detectSingle(Rcpp::NumericVector x){
    //init iterator
    Rcpp::NumericVector::iterator it = x.begin();

    //for index of changepoint
    //counting from 1, because R indexes from 1.
    int index = 1;
    //estimated location of changepoint
    int tauhat = 2;

    //flag to check if change has been found - used in while loop
    bool changeNotFound = true;

    while(  (changeNotFound) and (it != x.end())  ) {
        update(*it);
        if (getChangeDetected()){
            tauhat = index;
            changeNotFound = false;
        }            
        index++;
        ++it;
    }
    //only returning single changepoint
    return Rcpp::List::create(Rcpp::Named(CHANGEPOINT_FIELD_NAME) = tauhat); 
}


//detect single change WITH prechange mean and standard deviation known
Rcpp::List Detector::detectSinglePrechangeKnown(Rcpp::NumericVector x, 
                                                double prechangeMean, 
                                                double prechangeSigma){
    //set BL to be zero (no need for BL)
    setBL(0);
    //no burn-in period
    inBurnIn = false;

    //set prechange mean and sigma
    setStreamEstMean(prechangeMean);
    setStreamEstSigma(prechangeSigma);


    //now just do single change as usual
    return detectSingle(x);
}

//
//detect single change BUT with prechange parameters known
//
//does not really make sense for multiple changes to have prechange parameters known...
//since one does not typically know how many changepoints there are
//
//Do I add reset option for detector?


// RCPP_MODULE(detectormodule){
//     Rcpp::class_<Detector>( "Detector" )
//         .constructor<int>()
//         .field( "BL", &Detector::getBL, &Detector::setBL, "documentation for BL")
//         .method("getStreamEstMean", &Detector::getStreamEstMean, "documentation for streamEstMean")
//         .method("setStreamEstMean", &Detector::setStreamEstMean, "documentation for streamEstMean")
//         ;
// 
// }
// 

#endif
