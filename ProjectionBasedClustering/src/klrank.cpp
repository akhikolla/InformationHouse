#include <Rcpp.h>
#include "exception.h"
#include "distancematrix.h"
#include "euclidean.h"
#include "dataset.h"
#include "nervprobability.h"
#include "rankmatrix.h"
using namespace Rcpp;

using dredviz::DataMatrix;
using dredviz::DataSet;
using dredviz::DistanceMatrix;
using dredviz::Euclidean;
using dredviz::NeRVProbability;
using dredviz::RankMatrix;

// [[Rcpp::export]]
List klrank(NumericMatrix Data, NumericMatrix pData, int NeighborhoodSize=20) {
  /* klrank(Data, pData, NeighborhoodSize)
  * Compares the projection in pData with the original data in Data
  * and calculates the rank-based smoothed recall/precision.
  * 
  * INPUT
  * Data		          Matrix of original data
  * pData             Matrix of projected data
  * NeighborhoodSize  Sets the 'effective number of neighbors' used
  * to control the width of the Gaussian.
  *
  * OUTPUT
  * rank-based smoothed recall and smooted precision
  * all results are the results name substracted from one
  * i.e.: precision_best is (1-precision_best) in reality!
  *
  * AUTOR
  * 08/2016 FP
  */
  DataMatrix mat (Data);
  DataMatrix pmat (pData);
  DataSet    origData (mat);
  DataSet    projData (pmat);
  size_t effectiveNeighbors = NeighborhoodSize;
  DistanceMatrix* origDataDist = NULL;
  DistanceMatrix* projDataDist = NULL;
  
  RankMatrix* origRank = NULL;
  RankMatrix* projRankBest = NULL;
  RankMatrix* projRankWorst = NULL;
  RankMatrix* origRankReversed = NULL;
  
  Euclidean metric;
  
  origDataDist = new DistanceMatrix(origData, metric);
  origDataDist->scale (1.0 / origDataDist->getAverage());
  
  origRankReversed = new RankMatrix();
  origRank = new RankMatrix(*origDataDist, *origRankReversed);
  
  origRank->scale (1.0 / origRank->getAverage());
  origRankReversed->scale (1.0 / origRankReversed->getAverage());
  
  projDataDist = new DistanceMatrix(projData, metric);
  projDataDist->scale (1.0 / projDataDist->getAverage());
 
  projRankBest = new RankMatrix(*projDataDist, *origRank, true);
 
  projRankBest->scale(1.0 / projRankBest->getAverage());
 
  projRankWorst = new RankMatrix(*projDataDist, *origRank, false);
 
  projRankWorst->scale(1.0 / projRankWorst->getAverage());
 
  delete projDataDist;
  projDataDist = NULL;
 
 
  NeRVProbability inputProb(*origRank);
 
  NeRVProbability outputProbBest(*projRankBest);
 
  NeRVProbability outputProbWorst(*projRankWorst);
 
  NeRVProbability inputProbReversed(*origRankReversed);
 
 
  std::vector<double> sigma(inputProb.getRows());
 
  effectiveNeighbors = NeighborhoodSize;
 
  double temp = inputProb.findSigma(effectiveNeighbors, 0);
 
  for (size_t i = 0; i < sigma.size(); i++)
  {
    sigma[i] = temp;
  }
  
  inputProb.update(sigma);
  outputProbBest.update(sigma);
  outputProbWorst.update(sigma);
  inputProbReversed.update(sigma);
  
  delete origRankReversed;
  origRankReversed = NULL;
 
 
  //Calculate KL divergences
 
  double ED_p_q_best = 0.0;
  double ED_q_p_best = 0.0;
 
  double ED_p_q_worst = 0.0;
  double ED_q_p_worst = 0.0;
 
  double ED_p_q_max = 0.0;
  double ED_q_p_max = 0.0;
 
  for (size_t i = 0; i < inputProb.getRows(); i++)
   for (size_t j = 0; j < inputProb.getRows(); j++)
   {
     if (j != i)
     {
       // Best case:
       ED_p_q_best += inputProb(i,j)
       * std::log(inputProb(i,j) / outputProbBest(i,j));
       ED_q_p_best += outputProbBest(i,j)
         * std::log(outputProbBest(i,j) / inputProb(i,j));
       
       // Worst case:
       ED_p_q_worst += inputProb(i,j)
         * std::log(inputProb(i,j) / outputProbWorst(i,j));
       ED_q_p_worst += outputProbWorst(i,j)
         * std::log(outputProbWorst(i,j) / inputProb(i,j));
       
       // Maximum error (for normalization):
       
       ED_p_q_max += inputProb(i,j)
         * std::log(inputProb(i,j) / inputProbReversed(i,j));
       ED_q_p_max += inputProbReversed(i,j)
         * std::log(inputProbReversed(i,j) / inputProb(i,j));
       
    }
  }
   
   /* Normalize measures between 0 and 1 */
   
  ED_p_q_best /= ED_p_q_max;
  ED_p_q_worst /= ED_p_q_max;
  ED_q_p_best /= ED_q_p_max;
  ED_q_p_worst /= ED_q_p_max;
 
  delete origRank;
  delete projRankBest;
  delete projRankWorst;
  delete origDataDist;

  
  return Rcpp::List::create(
    Rcpp::Named("precision_best")    = ED_p_q_best,
    Rcpp::Named("precision_worst")       = ED_p_q_worst,
    Rcpp::Named("recall_best")    = ED_q_p_best,
    Rcpp::Named("recall_worst")       = ED_q_p_worst
  ) ; 
}


