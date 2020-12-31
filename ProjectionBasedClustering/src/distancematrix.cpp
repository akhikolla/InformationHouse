/*Copyright (C) Kristian Nybo, Jarkko Venna
 *
 *This software is released under the GNU Lesser General Public
 *License. See the included file LICENSE for details.*/

#include "distancematrix.h"
#include "euclidean.h"

namespace dredviz {
DistanceMatrix::DistanceMatrix (const DataMatrix & dataPoints,
                                Metric & metric):
    DataMatrix (dataPoints.getRows (), dataPoints.getRows ())
{
  CalculateDistanceMatrix calc;

  calc (dataPoints, metric, *this);
}


DistanceMatrix::DistanceMatrix (const DataMatrix & dataPoints):
    DataMatrix (dataPoints.getRows (), dataPoints.getRows ())
{
  CalculateDistanceMatrix calc;
  Euclidean metric;

  calc (dataPoints, metric, *this);
}



double
DistanceMatrix::getAverage () const
{
  double average = 0.0;

  for (size_t i = 0; i < getRows (); i++)
    for (size_t j = 0; j < getCols (); j++)
      if (i != j)
        average += (*this) (i, j);

  return average / ((getRows () - 1) * getCols ());
}


double
DistanceMatrix::getMax () const
{
  double max = (*this) (0, 0);

  for (size_t i = 0; i < getRows (); i++)
    for (size_t j = 0; j < getCols (); j++)
      if (i != j && (*this) (i, j) > max)
        max = (*this) (i, j);

  return max;
}

double
DistanceMatrix::getMin () const
{
  double min = (*this) (0, 0);

  for (size_t i = 0; i < getRows (); i++)
    for (size_t j = 0; j < getCols (); j++)
      if (i != j && (*this) (i, j) < min)
        min = (*this) (i, j);

  return min;
}

void
DistanceMatrix::scale (double factor)
{
  for (size_t i = 0; i < getRows (); i++)
    for (size_t j = 0; j < getCols (); j++)
      if (i != j)
        (*this) (i, j) *= factor;
}
}
