/*Copyright (C) Kristian Nybo, Jarkko Venna
 *
 *This software is released under the GNU Lesser General Public
 *License. See the included file LICENSE for details.*/

#ifndef METRIC_HH
#define METRIC_HH

#include "dataset.h"

namespace dredviz {
class Metric
{
public:
  virtual ~ Metric ()
  {
  }

  //Returns the distance between the points data(row1, col1) and
  //data(row2, col2).
  virtual double operator ()
  (const DataMatrix & data, size_t row1, size_t row2) = 0;

  virtual bool isSymmetric (void) const = 0;
};
}

#endif
