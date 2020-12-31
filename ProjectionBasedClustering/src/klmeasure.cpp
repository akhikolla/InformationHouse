#include <Rcpp.h>
#include "datamatrix.h"
#include "dataset.h"
#include "distancematrix.h"
#include "metric.h"
#include "euclidean.h"
#include "nervprobability.h"

using namespace Rcpp;

using dredviz::DataMatrix;
using dredviz::DataSet;
using dredviz::DistanceMatrix;
using dredviz::Euclidean;
using dredviz::NeRVProbability;

// [[Rcpp::export]]
List c_klmeasure(NumericMatrix Data, NumericMatrix pData, int NeighborhoodSize = 20) {
  DataMatrix mat (Data);
  DataMatrix pmat (pData);
  DataSet    origData (mat);
  DataSet    projData (pmat);
  size_t effectiveNeighbors = NeighborhoodSize;
  DistanceMatrix origDataDist;
  DistanceMatrix projDataDist;
  // TODO: Option einbauen auch direkt eine Distanzmatrix zu ?bergeben
  // dazu dann konvertierung von R-Distanzmatrix zu DistanceMatrix
  
  Euclidean metric;
  DistanceMatrix dmorig(origData, metric);
  origDataDist = dmorig;
  origDataDist.scale (1.0 / origDataDist.getAverage());
  
  DistanceMatrix dmpro(projData, metric);
  projDataDist = dmpro;
  projDataDist.scale (1.0 / projDataDist.getAverage());
  
  NeRVProbability origProb(origDataDist);
  NeRVProbability projProb(projDataDist);

  std::vector<double> sigma(origProb.getRows());
  origProb.findSigma(sigma, effectiveNeighbors);

  origProb.update(sigma);
  projProb.update(sigma);
  
  //Calculate Kullback-Leibler divergences

  double ED_p_q = 0.0;
  double ED_q_p = 0.0;

  for (size_t i = 0; i < origProb.getRows(); i++)
    for (size_t j = 0; j < origProb.getRows(); j++)
    {
      if (j != i)
      {
        ED_p_q += origProb(i,j) * std::log(origProb(i,j) / projProb(i,j));
        ED_q_p += projProb(i,j) * std::log(projProb(i,j) / origProb(i,j));
      }
    }
  
  
  return Rcpp::List::create(
    Rcpp::Named("SmoothedPrecision") = ED_q_p,
    Rcpp::Named("SmoothedRecall")       = ED_p_q
  ) ;
}


