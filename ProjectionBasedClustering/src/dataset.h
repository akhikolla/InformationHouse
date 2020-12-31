/*Copyright (C) Kristian Nybo, Jarkko Venna
 *
 *This software is released under the GNU Lesser General Public
 *License. See the included file LICENSE for details.*/

#ifndef DATASET_HH
#define DATASET_HH

#include "datamatrix.h"

#include <vector>
#include <string>

using
std::vector;
using
std::string;

namespace dredviz {
class
  DataSet: public DataMatrix
{
protected:
  vector<vector<string> > labels;
  vector <string> axisLabels;

public:

  /* Methods specific to DataSet */

  DataSet ();
  DataSet (const DataMatrix & data);
  DataSet (size_t row, size_t col) : DataMatrix(row, col) {}

  virtual ~DataSet () {};

  void addLabels (vector < vector < string > >&labels);
  void addAxisLabels (vector < string > &labels);

  const vector<vector<string> >& getLabels () const;

  const vector<string>& getAxisLabels () const;


  DataSet & operator= (const DataSet & data);



  //Returns the number of labels per data point.

  bool hasLabels () const
  {
    return false;
  }                             //placeholder

  bool hasAxisLabels () const
  {
    return false;
  }
};
}

#endif
