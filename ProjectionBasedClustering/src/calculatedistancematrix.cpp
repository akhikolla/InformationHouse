/*Copyright (C) Kristian Nybo, Jarkko Venna
 *
 *This software is released under the GNU Lesser General Public
 *License. See the included file LICENSE for details.*/

#include "calculatedistancematrix.h"

namespace dredviz {

void
CalculateDistanceMatrix::operator () (const DataMatrix & data,
                                      Metric & metric, DataMatrix & target)
{
  //  Matrix distmat(data.getRows(), data.getRows());


  if (metric.isSymmetric ())
  {
    for (size_t i = 0; i < data.getRows (); i++)
      for (size_t j = i; j < data.getRows (); j++)
        if (i == j)
          target (i, j) = 0;
        else
          target (i, j) = target (j, i) = metric (data, i, j);
  }
  else
  {
    for (size_t i = 0; i < data.getRows (); i++)
      for (size_t j = 0; j < data.getRows (); j++)
        if (i == j)
          target (i, j) = 0;
        else
          target (i, j) = target (j, i) = metric (data, i, j);
  }

  //  target = distmat;
}
}
