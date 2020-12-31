#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

// [[Rcpp::export]]
void esomTrainedWeightVectorsMexicanHatC(NumericMatrix WeightVectors, NumericVector DataPoint,NumericVector indices,
                                         NumericVector DistancesToBm, double Radius,
                                         double LearningRate) {
  
  // source: http://www.codeproject.com/Articles/16273/Self-Organizing-Feature-Maps-Kohonen-maps
  int ncol = WeightVectors.ncol();

  double MexicanHatNorm = 1.0 / ((double)Radius+1);


  int l = indices.length();

  for(int j=0; j < l; j++){
    int ind = indices[j] - 1;
    double square = pow(MexicanHatNorm * (double)DistancesToBm[j],2);
    double Scaling = (1.0 - 2.0*square) * exp(-square); 

 //Rcout << Scaling << std::endl;//WeightVectors(ind,0) << ", " << WeightVectors(ind,1) << " mit "<< Scaling << std::endl;

    for(int i=0; i < ncol; i++){ //durch die einzelnen komponenten
      WeightVectors(ind,i) = WeightVectors(ind,i) + (LearningRate * Scaling * (DataPoint[i] - WeightVectors(ind,i)));
    }
  }

}
