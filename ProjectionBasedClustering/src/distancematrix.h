/*Copyright (C) Kristian Nybo, Jarkko Venna
 *
 *This software is released under the GNU Lesser General Public
 *License. See the included file LICENSE for details.*/

#ifndef DISTANCEMATRIX_HH
#define DISTANCEMATRIX_HH

#include "metric.h"
#include "calculatedistancematrix.h"
#include "datamatrix.h"

namespace dredviz {

class DistanceMatrix:public DataMatrix
{

public:

  /* Creates a dim * dim UNINITIALIZED distance matrix. */

  DistanceMatrix (size_t dim = 1) : DataMatrix(dim,dim) {}

  DistanceMatrix (const DataMatrix & dataPoints, Metric & metric);

  /* If no metric is specified, euclidean distance will be used. */

  DistanceMatrix (const DataMatrix & dataPoints);

  /* These work like the corresponding DataMatrix functions, except that they
     ignore the diagonal elements, because we're not interested in a point's
     distance from itself. */

  double getMin () const;
  double getMax () const;
  double getAverage () const;
  void scale (double scalar);
};
}

#endif
