#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

// [[Rcpp::export]]
void esomTrainedWeightVectorsGaussC(NumericMatrix WeightVectors, NumericVector DataPoint,NumericVector indices,NumericVector DistancesToBm, double Radius, double LearningRate) {
  int ncol = WeightVectors.ncol();

  int Stddevs = 2;

  double Gaussnorm = (2.0 * std::pow(Radius + 1,2)) / std::pow(Stddevs,2);
  int l = indices.length();

  for(int j=0; j < l; j++){
    int ind = indices[j] - 1;
    double Scaling = std::exp((-std::pow(((double)DistancesToBm[j]),2))/Gaussnorm);

    //Rcout << WeightVectors(ind,0) << ", " << WeightVectors(ind,1) << " mit "<< Scaling << std::endl;

    for(int i=0; i < ncol; i++){ //durch die einzelnen komponenten
      WeightVectors(ind,i) = WeightVectors(ind,i) + (Scaling * LearningRate * (DataPoint[i] - WeightVectors(ind,i)));
    }
  }

}
