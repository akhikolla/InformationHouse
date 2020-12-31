/*Copyright (C) Kristian Nybo, Jarkko Venna
 *
 *This software is released under the GNU Lesser General Public
 *License. See the included file LICENSE for details.*/

#include <iostream>
#include <limits>
#include <cmath>
#include "inputprobentropy.h"

namespace dredviz {
double
InputProbEntropy::evaluate (const DataMatrix & sigma)
{
  if (sigma (0, 0) <= 0.0)
    return std::numeric_limits < double >::max ();

  prob.update (index, sigma (0, 0));

  double entropy = 0.0;

  for (size_t j = 0; j < prob.getCols (); j++)
  {
    if (j != index)
    {
      entropy -= prob (index, j) * std::log (prob (index, j)) / std::log (2.0);
    }
  }

  //std::cout << "Entropy error " << fabs(entropy - desiredEntropy) << "\n";

  return fabs (entropy - desiredEntropy);
}


InputProbEntropy::InputProbEntropy (size_t effectiveNeighborhoodSize, const size_t index, NeRVProbability & prob):desiredEntropy ((effectiveNeighborhoodSize <
        2) ? 1.0 : std::log (effectiveNeighborhoodSize) / std::log (2.0)),
    index (index), prob (prob)
{
}
}
