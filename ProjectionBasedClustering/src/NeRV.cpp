#include <Rcpp.h>
using namespace Rcpp;
#include "datamatrix.h"
#include "dataset.h"
#include "distancematrix.h"
#include "metric.h"
#include "euclidean.h"
#include "randomdatagenerator.h"
#include "nervcostfunction.h"
#include "nervoptstrat.h"
#include "goldensectionsearch.h"
#include "conjugategradientopt.h"

using dredviz::DataMatrix;
using dredviz::DataSet;
using dredviz::DistanceMatrix;
using dredviz::Euclidean;
using dredviz::DynamicDouble;
using dredviz::GoldenSectionSearch;
using dredviz::NeRVCostFunction;
using dredviz::ConjugateGradientOpt;
using dredviz::NeRVOptStrat;
using dredviz::RandomDataGenerator;

NumericMatrix DataSet2NumericMatrix(DataSet input){
  NumericMatrix output = NumericMatrix(input.getRows(),input.getCols());
  for(size_t i = 0; i < input.getRows(); i++){
    for(size_t j = 0; j < input.getCols(); j++){
      output(i,j) = input(i,j);
    }
  }
  return output;
}

// (Data, lambda = 0.1, neighbors = 20, iterations = 10, cg_steps = 2, cg_steps_final = 40,randomORPCAinit=T)
// Function pca wird übergeben, damit in RCPP die R eigene PCA funktion aufgerufen werden kann
// Grund: die von NeRV braucht eine eigene externe Library. Portieren dauert dann zu lange.
// [[Rcpp::export]]
NumericMatrix c_NeRV(NumericMatrix data, double lambda, int lastNeighbor, int iterations, int stepsPerRound, int stepsOnLastRound, bool randominit, int outputDimension, Function pca ) {
  DistanceMatrix origDataDist;
  ConjugateGradientOpt *optStep = NULL;
  std::string tmpoutputfile;
  
  // Step 1: Converting to C++ types
  DataMatrix mat (data);
  DataSet    origData (mat);
  Euclidean metric;
  
  DistanceMatrix dm (origData, metric);
  origDataDist = dm;
  
  // Step 2: initializing matrix for projected data
  DataSet    projData;
  if(!randominit){
    NumericMatrix projectedData;
    List pcares = pca(data,outputDimension);
    projectedData = NumericMatrix((SEXP)pcares["ProjectedPoints"] );
    DataMatrix pmat(projectedData);
    projData = DataSet(pmat);
  } else {
    RandomDataGenerator randgen(origDataDist.getRows(), outputDimension, 1);
    randgen.loadData(projData);
  }

  // Step 3: execute the algorithm like in the original NeRV.cc
  std::vector<double> weights(origDataDist.getRows(), 1);

  origDataDist.scale (1.0 / origDataDist.getAverage ());

  double initrad = 0.0;
  initrad = origDataDist.getMax () / 2.0;

  DynamicDouble initialRadius (initrad, 0);

  GoldenSectionSearch linesearch;

  NeRVCostFunction costFunc (origDataDist, projData, linesearch,
                             initialRadius, lambda, lastNeighbor,
                             weights, Rcpp::Rcout);

  optStep = new ConjugateGradientOpt(costFunc, linesearch, Rcpp::Rcerr);

  NeRVOptStrat optStrategy (iterations, stepsPerRound, stepsOnLastRound,
                            tmpoutputfile);

  optStrategy.optimize (origDataDist, projData, *optStep, costFunc,
                        Rcpp::Rcerr);

  // Step 4: Ergebnis zurückkonvertieren

  delete optStep;
  return DataSet2NumericMatrix(projData);
}
