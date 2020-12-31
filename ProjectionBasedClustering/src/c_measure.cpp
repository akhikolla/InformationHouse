#include <Rcpp.h>
#include "measure.h"
#include "dataset.h"
#include "distancematrix.h"
#include "metric.h"
#include "euclidean.h"
#include "randomdatagenerator.h"
#include "conttrust.h"

using namespace Rcpp;

using dredviz::DataMatrix;
using dredviz::DataSet;
using dredviz::DistanceMatrix;
using dredviz::Euclidean;
using dredviz::ContTrust;

// [[Rcpp::export]]
NumericMatrix c_measure(NumericMatrix datamat, NumericMatrix projmat, unsigned int lastNeighbor) {
  DistanceMatrix origDataDist;
  DistanceMatrix projDataDist;
  Euclidean metric;
  
  DataMatrix mat (datamat);
  DataSet    origData (mat);
  
  DistanceMatrix dm (origData, metric);
  origDataDist = dm;
  
  DataMatrix pmat (projmat);
  DataSet    projData (pmat);
    
  DistanceMatrix dmproj (projData, metric);
  projDataDist = dmproj;  
  
  ContTrust ctmeasure (lastNeighbor);
 
  return ctmeasure.measure(origDataDist, projDataDist);
}

