/*Copyright (C) Kristian Nybo, Jarkko Venna
 *
 *This software is released under the GNU Lesser General Public
 *License. See the included file LICENSE for details.*/

#ifndef NERV_PROBABILITY_HH
#define NERV_PROBABILITY_HH

#include <math.h>
#include "datamatrix.h"
#include "distancematrix.h"

namespace dredviz {
class NeRVProbability
{
private:
  DataMatrix prob;
  const DistanceMatrix & data;
  double shift;
  const double minexp;                // Minimum exponent allowed in calculating the probability.
  const double MIN_SIGMA;

public:
  NeRVProbability (const DistanceMatrix & data, double minexp = -450.0);

  void update (size_t row, double sigma2);

  void update (const vector < double >&sigma2);

  double operator () (size_t i, size_t j) const
  {
    return prob (i, j);
  };

  double findSigma (size_t effectiveNeighbors, size_t index);

  void findSigma (vector < double >&sigma2, size_t effectiveNeighbors);

  size_t getRows () const
  {
    return prob.getRows ();
  };

  size_t getCols () const
  {
    return prob.getCols ();
  };
};


}

#endif
