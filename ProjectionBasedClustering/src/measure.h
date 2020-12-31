/*Copyright (C) Kristian Nybo, Jarkko Venna
 *
 *This software is released under the GNU Lesser General Public
 *License. See the included file LICENSE for details.*/

#ifndef MEASURE_HH
#define MEASURE_HH

#include "distancematrix.h"
#include "exception.h"
#include <iostream>
#include <Rcpp.h>

using Rcpp::NumericMatrix;

namespace dredviz {

class Measure
{
public:
  virtual ~ Measure ()
  {
  }
  virtual NumericMatrix measure (const DistanceMatrix & origData,
                        const DistanceMatrix & projData) = 0;

class RowMismatch:public Exception
  {
  public:
    RowMismatch (std::string errMsg):Exception (errMsg)
    {
    }
  };
};

}

#endif
