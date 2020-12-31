/*Copyright (C) Kristian Nybo, Jarkko Venna
 *
 *This software is released under the GNU Lesser General Public
 *License. See the included file LICENSE for details.*/

#ifndef CONTTRUST_HH
#define CONTTRUST_HH

#include "measure.h"
#include "distancematrix.h"
#include "datamatrix.h"
#include <Rcpp.h>

using Rcpp::NumericMatrix;

namespace dredviz {
class ContTrust:public Measure
{
private:
  size_t maxNeighborhood;
  DataMatrix measures; //best case & worst case trustworthiness & continuity
//  DataMatrix precrec; //best case precision & recall

public:
  virtual ~ ContTrust ()
  {
  }
  ContTrust (size_t maxNeighborhood);

  NumericMatrix measure(const DistanceMatrix & origDist,
                        const DistanceMatrix & projDist);
  // wcase =0 average, wcase=1 best, wcase=2 worst
  double getContinuity(size_t numNeighbors, int wcase=0);
  double getTrustworthiness(size_t numNeighbors, int wcase=0);

};
}

#endif
