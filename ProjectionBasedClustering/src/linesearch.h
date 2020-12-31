/*Copyright (C) Kristian Nybo, Jarkko Venna
 *
 *This software is released under the GNU Lesser General Public
 *License. See the included file LICENSE for details.*/

#ifndef LINESEARCH_HH
#define LINESEARCH_HH

#include "costfunction.h"

namespace dredviz {
class LineSearch
{
public:
  /* Calculates the minimum of costFunction along the line that passes through
     pointA and pointB. Returns the minimum in pointA. */

  virtual ~ LineSearch ()
  {
  }
  virtual double operator () (CostFunction & costFunction, DataMatrix & point,
                              const DataMatrix & negativeGradient,
                              double initialStepSize, double &finalCost) = 0;
};
}

#endif
