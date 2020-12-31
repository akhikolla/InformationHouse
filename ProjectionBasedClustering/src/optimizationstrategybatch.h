/*Copyright (C) Kristian Nybo, Jarkko Venna
 *
 *This software is released under the GNU Lesser General Public
 *License. See the included file LICENSE for details.*/

#ifndef OPTIMIZATIONSTRATEGYBATCH_HH
#define OPTIMIZATIONSTRATEGYBATCH_HH

#include "dataset.h"
#include "distancematrix.h"
#include "optimizationstepbatch.h"
#include "costfunction.h"

#include <ostream>

namespace dredviz {
class OptimizationStrategyBatch
{
public:
  /* optimize() attempts to find a configuration of points in the desired
     output space that minimizes the cost function.


     const DistanceMatrix& origDist: The original data as a distance matrix,
     that is, origDist(i,j) is the distance between points i and j in the
     original data.

     Matrix& initialProjData: a matrix that represents the original data
     projected into the output space: initialProjData(i,j) is the jth
     coordinate of the ith data point in the output space.

     OptimizationStepBatch& optStep: The algorithm that performs a batch of
     optimization steps.

     CostFunction& costFunc: The cost function to be minimized.

     std::ostream feedback: optimize() will report its progress into this
     stream.
   */

  virtual void optimize (const DistanceMatrix & origDist,
                         DataMatrix & initialProjData,
                         OptimizationStepBatch & optStep,
                         CostFunction & costFunc, std::ostream & feedback) =
                           0;

  virtual ~ OptimizationStrategyBatch ()
  {
  }
};

}
#endif
