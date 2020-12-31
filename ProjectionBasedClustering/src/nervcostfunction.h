/*Copyright (C) Kristian Nybo, Jarkko Venna
 *
 *This software is released under the GNU Lesser General Public
 *License. See the included file LICENSE for details.*/

#ifndef NERVCOSTFUNCTION_HH
#define NERVCOSTFUNCTION_HH

#include "linesearch.h"
#include "distancematrix.h"
#include "dataset.h"
#include "costfunction.h"
#include "dynamicdouble.h"
#include "nervprobability.h"

#define MIN_PROB 1.0e-200       /* The minimimum probability allowed. */

#include <iostream>
#include <vector>

/* NeRVCostFunction implements the cost function used in NeRV, as described
   in "Nonlinear dimensionality reduction as information retrieval" by
   Jarkko Venna and Samuel Kaski. */

namespace dredviz {
class NeRVCostFunction:public CostFunction
{
private:
  const DistanceMatrix & origDist;
  //  DistanceMatrix projDist;

  double minexp;
  /* radius is used as the gaussian when calculating input and output
     probabilities... */

  DynamicDouble radius;

  double lambda;

  /* ...except if radius is smaller than finalNeighborhoods[i], in which case
     the point i will use finalNeighborhoods[i] instead of radius.
     finalNeighborhoods[i] will contain a value such that when it is used
     as the gaussian for the neighborhood probability distribution of point i,
     the entropy of the distribution will be log(effectiveNeighborhoodSize),
     where effectiveNeighborhoodSize is the neighborhoodSize given to
     NeRVCostFunction's constructor. */

  vector < double >finalNeighborhoods;

  /* For convenience, the squares of the current gaussians are always stored
     in this vector, ie., sigma[i] equals
     max(radius.value(), finalNeighborhoods[i]). */

  vector < double >sigmaSqrd;

  /* The line search used in calculateFinalNeighborhoods() */

  LineSearch & linesearch;


  /* inputProb(i,j), a double precision value between 0 and 1, represents
     the probability that point i would pick point j as its neighbor in
     the input data. Analogously, inputProb(i,j) represents the probability
     that point i would pick point j as its neighbor in the output data. */

  NeRVProbability inputProb;
  //  NeRVProbability outputProb;
  DataMatrix outputProb;

  /* The weights for each data point's contribution to the cost function.
     Defaults to 1 for each data point. */
  std::vector<double> weights;

  std::ostream & feedback;

  /* All output distances are shifted by this amount to preserve numerical
     precision (the exponent function becomes numerically unstable for small
     values). */

  double minimumDistance;

  /* The smallest distance in the original data. */

  const double origMinimumDistance;

  /* Temporary values used to calculate the derivative */
  DataMatrix dDval;
  std::vector<double> w;
  //  vector<DataMatrix> ddist;

  /* Initializes finalNeighborhoods. */

  void calculateFinalNeighborhoods (size_t effectiveNeighborhoodSize);

  /* Updates the entire input probability matrix. */

  //  void updateInputProb();

  /* Updates the rowth row of the input probability matrix, using sigma as
     the radius. */

  //  void updateInputProb(size_t row, double sigma);

  /* Updates the entire output probability matrix to represent the data in
     projData. */

  void updateOutputProb (const DataMatrix & projData);
  //  void updateGradientParts(const DataMatrix& projData);
  /* Updates the distance by which all the output distances are shifted
     to preserve numerical precision. */

  void updateMinimumDistance (const DataMatrix & projData);


public:
  NeRVCostFunction (const DistanceMatrix & origDist,
                    DataMatrix & projData, LineSearch & lineSearch,
                    DynamicDouble radius, double lambda,
                    size_t neighborhoodSize,
                    const std::vector<double> &weights,
                    std::ostream & feedback);

  /* See CostFunction for the rest. */

  double getGradient (const DataMatrix & projData, DataMatrix & target);

  double evaluate (const DataMatrix & data);

  void reportParameters(std::string& target);
  void updateDynamicParameters (size_t currentRound, size_t totalRounds,
                                const DataMatrix & projData);

  void updateDataRepresentation (const DataMatrix & newData);

};
}

#endif
