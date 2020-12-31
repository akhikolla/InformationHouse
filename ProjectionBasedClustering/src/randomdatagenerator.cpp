/*Copyright (C) Kristian Nybo, Jarkko Venna
 *
 *This software is released under the GNU Lesser General Public
 *License. See the included file LICENSE for details.*/

#include <stdlib.h>
#include <cmath>
#include <Rmath.h>
#include <Rcpp.h>
#include "randomdatagenerator.h"

namespace dredviz {

RandomDataGenerator::RandomDataGenerator (size_t numberOfPoints,
    size_t dimension, double range):
    range (range),
    matrix (numberOfPoints, dimension)
{
  Rcpp::Rcout << "Using current .Random.seed as RNG seed.\n";
}


void
RandomDataGenerator::loadData (DataSet & target)
{
  GetRNGstate();
  
  for (size_t i = 0; i < matrix.getRows (); i++)
    for (size_t j = 0; j < matrix.getCols (); j++)
      matrix (i, j) = R::unif_rand() * range;
  
  PutRNGstate();
  DataSet data (matrix);

  target = data;
}
}
