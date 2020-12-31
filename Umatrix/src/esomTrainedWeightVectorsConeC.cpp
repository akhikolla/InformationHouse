#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
void esomTrainedWeightVectorsConeC(NumericMatrix WeightVectors, NumericVector DataPoint,NumericVector indices,NumericVector DistancesToBm, double Radius, double LearningRate) {
  int ncol = WeightVectors.ncol();

  int l = indices.length();

  for(int j=0; j < l; j++){
    int ind = indices[j] - 1;
    double Scaling = ((Radius + 1) - DistancesToBm[j])/(Radius + 1);

    //Rcout << WeightVectors(ind,0) << ", " << WeightVectors(ind,1) << " mit "<< Scaling << std::endl;

    for(int i=0; i < ncol; i++){ //durch die einzelnen komponenten
      WeightVectors(ind,i) = WeightVectors(ind,i) + (Scaling*LearningRate * (DataPoint[i] - WeightVectors(ind,i)));
    }
  }

}
