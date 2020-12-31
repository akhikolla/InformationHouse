/*Copyright (C) Kristian Nybo, Jarkko Venna
 *
 *This software is released under the GNU Lesser General Public
 *License. See the included file LICENSE for details.*/

#ifndef EUCLIDEANSQUARED_HH
#define EUCLIDEANSQUARED_HH

#include "metric.h"

namespace dredviz {
class EuclideanSquared:public Metric
{
private:
  bool isSymmetric (void) const
  {
    return true;
  }

public:
  double operator () (const DataMatrix & data, size_t row1, size_t row2);
};
}

#endif
