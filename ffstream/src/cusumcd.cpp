#ifndef GUARD_cusumcd_cpp
#define GUARD_cusumcd_cpp


#include "cusumcd.h"
//#include "fff.h"
//#include<Rcpp.h>
//#include "utils.h"
//#include<vector>


//constructor
CusumChangeDetector::CusumChangeDetector(double k_, double h_, int BL_) :
    Detector(BL_),
    S(INIT_CUSUM_S),
    T(INIT_CUSUM_T),
    h(h_),
    k(k_){}




//print the CUSUM values
void CusumChangeDetector::print(){
    Rcpp::Rcout << "k: " << getK() <<  ", h: " << getH() << std::endl;
    Rcpp::Rcout << "S: " << getS() << ", T: " << getT() << std::endl;
    Rcpp::Rcout << "Burn-in count: " << BLcount << std::endl;
    Rcpp::Rcout << "Burn-in length: " << getBL() << std::endl;
    Rcpp::Rcout <<  "changeDetected: " << getChangeDetected() << std::endl;
}




//standardise the observation for CUSUM
double CusumChangeDetector::standardiseObs(double obs, double mu, double sigma){
    if (sigma==0){
       obs = 0; 
    } else {
        obs = (obs - mu) / sigma;
    }
    return ( obs );
}


//update the S and T CUSUMs
//obshat = standardised obs
//S_{t+1} = max{0, S_t + obshat - k }
void CusumChangeDetector::cusumUpdate(double obs){
    //standardise obs
    double stdObs = standardiseObs( obs, getStreamEstMean(), getStreamEstSigma() );

    //update S
    S = S + stdObs - k;
    if (S < 0)
        S = 0;
    //update T
    T = T - stdObs - k;
    if (T < 0)
        T = 0;

    //compute the pval
    computePvalue();
}



//start the burn-in
void CusumChangeDetector::startBurnIn(){
    Detector::startBurnIn();

    //reset the CUSUM and burn-in counter
    S = INIT_CUSUM_S;
    T = INIT_CUSUM_T;
}


//check if there is a change
void CusumChangeDetector::checkIfChange(){
    if (S > h)
        changeDetected = true;
    if (T > h)
        changeDetected = true;
}


//compute the pvalue
//will combine two one-sided tests into a single value
void CusumChangeDetector::computePvalue(){
    //we already have h and S and T
    //p1 from S
    //p2 from T
    double p1 = computeOneSidedPvalue(S, 0, h);
    double p2 = computeOneSidedPvalue(T, 0, h);

    //reverse bonferonni
    //combine p1 and p2 into p3 using p3 = 2 * min(p1, p2)
    setPval( combineTwoOneSidedPvalues(p1, p2) );
}

//main update function for CUSUM
void CusumChangeDetector::update(double obs){
    //if a change occurred on the last update, need to reset to burnIn
    if (changeDetected){
        startBurnIn();
    }

    //go to burn in
    if (inBurnIn){
        streamEstimator.update(obs);
        BLcount++;
        if ( !(BLcount < BL) ){
            //using Detector::stopBurnIn()
            stopBurnIn();
        }
    } else {
        cusumUpdate(obs); 
        checkIfChange();
    }
}


//getter for S
double CusumChangeDetector::getS(){
    return S;
}


//getter for T
double CusumChangeDetector::getT(){
    return T;
}


//getter for h
double CusumChangeDetector::getH(){
    return h;
}


//getter for k
double CusumChangeDetector::getK(){
    return k;
}




//process a vector
Rcpp::List CusumChangeDetector::processVectorSave (Rcpp::NumericVector vec){
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
            //no more push_back, save to particular location
            numChangepoints++;
            //need to pushback i+1; index starts from 0 in cpp, but from 1 in R
            changepoints[numChangepoints-1] = i+1;
        }

    }
    //now wrap subvector
    return Rcpp::List::create(Rcpp::Named(IS_CHANGE_DETECTED_FIELD_NAME)=Rcpp::wrap(isChangeVec),
                              Rcpp::Named(CHANGEPOINT_FIELD_NAME)=Rcpp::wrap(std::vector<int>(changepoints.begin(), changepoints.begin()+numChangepoints)));
}

//-----------------------------------------------------------------------//



//creation of module
RCPP_MODULE(cusumcdmodule){

    Rcpp::class_<Detector>( "Detector" )
        .constructor<int>()
        .property( "BL", &Detector::getBL, &Detector::setBL, "documentation for BL")
        .property( "pval", &Detector::getPval, "documentation for pval")
        .property("streamEstMean", &Detector::getStreamEstMean, &Detector::setStreamEstMean, "documentation for streamEstMean")
        .property("streamEstSigma", &Detector::getStreamEstSigma, &Detector::setStreamEstSigma, "documentation for streamEstSigma")
        .property("changeDetected", &Detector::getChangeDetected, "documentation for changeDetected")
        ;
    
    Rcpp::class_<CusumChangeDetector>( "CusumChangeDetector" )
        .derives<Detector>("Detector")
        .constructor<double, double, double>()
        .method( "update" , &CusumChangeDetector::update, "documentation for update")
        .method( "print" , &CusumChangeDetector::print, "documentation for print")
        .method( "processVectorSave" , &CusumChangeDetector::processVectorSave, "documentation for processVectorSave")
//         .property("S", &CusumChangeDetector::getS, "documentation for S")
//         .property("T", &CusumChangeDetector::getT, "documentation for T")
        .property("h", &CusumChangeDetector::getH, "documentation for h")
        .property("k", &CusumChangeDetector::getK, "documentation for k")

//         .property("S", &CusumChangeDetector::getS, "documentation for S")
//         .property("T", &CusumChangeDetector::getT, "documentation for T")
        ;
}

#endif
