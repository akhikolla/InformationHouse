/*Copyright (C) Kristian Nybo, Jarkko Venna
 *
 *This software is released under the GNU Lesser General Public
 *License. See the included file LICENSE for details.*/

#include "dataset.h"

namespace dredviz {
DataSet::DataSet ():DataMatrix (1, 1)
{
}

DataSet::DataSet (const DataMatrix & data):
    DataMatrix (data)
{
}


void
DataSet::addLabels (vector < vector < string > >&labels)
{
  this->labels = labels;
}


void
DataSet::addAxisLabels (vector < string > &labels)
{
  this->axisLabels = labels;
}


const vector < vector < string > >&
DataSet::getLabels () const
{
  return this->labels;
}


const vector < string > &
DataSet::getAxisLabels () const
{
  return this->axisLabels;
}

DataSet & DataSet::operator= (const DataSet & dataSet)
{
  if (this == &dataSet)
    return *this;

  /* Copy the DataMatrix part of dataSet */

  DataMatrix::operator= (dataSet);

  if (dataSet.hasLabels ())
    labels = dataSet.labels;

  if (dataSet.hasAxisLabels ())
    axisLabels = dataSet.axisLabels;

  return *this;
}
}
