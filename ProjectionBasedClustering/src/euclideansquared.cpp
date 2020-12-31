/*Copyright (C) Kristian Nybo, Jarkko Venna
 *
 *This software is released under the GNU Lesser General Public
 *License. See the included file LICENSE for details.*/

#include "euclideansquared.h"

namespace dredviz {
double
EuclideanSquared::operator () (const DataMatrix & data, size_t row1,
                               size_t row2)
{
  double distanceSquared = 0.0;

  for (size_t col = 0; col < data.getCols (); col++)
  {
    distanceSquared += (data (row1, col) - data (row2, col)) *
                       (data (row1, col) - data (row2, col));
  }

  return distanceSquared;
}
}
