/*Copyright (C) Kristian Nybo, Jarkko Venna
 *
 *This software is released under the GNU Lesser General Public
 *License. See the included file LICENSE for details.*/

#ifndef GOLDENSECTIONSEARCH_HH
#define GOLDENSECTIONSEARCH_HH

#include "linesearch.h"
#include "datamatrix.h"

#include <iostream>

namespace dredviz {
class GoldenSectionSearch:public LineSearch
{
private:
  /* The magnification ratio used in the initial bracketing. */

  const double INITIAL_BRACKETING_MAGNIFICATION;

  /* The ratio with which the length of an interval is multiplied when
     picking a new point. */

  const double RATIO;

  /* The line search will terminate once the minimum is known with a margin
     of error equal to or less than TOLERANCE. */

  const double TOLERANCE;

  const size_t MAX_ITER;

  /* referencePoint is chosen as the origin on the line along which we
     minimize, ie., all points are expressed as referencePoint + c * direction,
     where c is a scalar and direction is a vector (matrix) that points along
     the line. */

  DataMatrix referencePoint;

  double alpha;
  double beta;
  double gamma;
  double xi;

  double costA;                 /* Cost at the 'left' end of the bracket */
  double currentMinimumCost;
  double costC;                 /* Cost at the 'right' end of the bracket */
  double costX;

  DataMatrix x;


  //  DataMatrix direction;
  bool findInitialBracket (CostFunction & costFunction,
                           DataMatrix & pointA,
                           const DataMatrix & negativeGradient,
                           double initialStepSize);

public:
  GoldenSectionSearch ();

  double operator () (CostFunction & costFunction, DataMatrix & pointA,
                      const DataMatrix & negativeGradient,
                      double initialStepSize, double &finalCost);
};
}

#endif
