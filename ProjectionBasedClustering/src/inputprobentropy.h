/*Copyright (C) Kristian Nybo, Jarkko Venna
 *
 *This software is released under the GNU Lesser General Public
 *License. See the included file LICENSE for details.*/

#include "nervprobability.h"
#include "costfunction.h"

namespace dredviz {
class InputProbEntropy:public CostFunction
{
private:
  double desiredEntropy;
  size_t index;
  NeRVProbability & prob;

public:
  /* effectiveNeighborhoodSize -- the number of points that the final
     radius of each point should ('effectively') contain.

     const size_t* index -- this should be a pointer to the variable that
     indexes the rows of the probability matrix for which the entropy is
     being calculated. */

  InputProbEntropy (size_t effectiveNeighborhoodSize, size_t index,
                    NeRVProbability & prob);

  /* Calculates the entropy for the indexth row of targetCostFunc. */
  double evaluate (const DataMatrix & sigma);

  void setRow (size_t i)
  {
    index = i;
  }

  /* Dummy functions; we don't need the gradient, and there are no
     dynamic parameters to be updated. */

  double getGradient (const DataMatrix &, DataMatrix &)
  {
    return 0.0;
  }

  void reportParameters(std::string& target)
  {
    target = "";
  }


  void updateDynamicParameters (size_t, size_t, const DataMatrix &)
  {
  }
  void updateDataRepresentation (const DataMatrix &)
  {
  }
};
}
