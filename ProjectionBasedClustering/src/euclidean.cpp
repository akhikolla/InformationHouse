/*Copyright (C) Kristian Nybo, Jarkko Venna
 *
 *This software is released under the GNU Lesser General Public
 *License. See the included file LICENSE for details.*/

#include "euclidean.h"
#include <math.h>

namespace dredviz {
double
Euclidean::operator () (const DataMatrix & data, size_t row1, size_t row2)
{
  double distanceSquared = 0;

  for (size_t i = 0; i < data.getCols (); i++)
    distanceSquared += pow (data (row1, i) - data (row2, i), 2);

  return sqrt (distanceSquared);
}
}
