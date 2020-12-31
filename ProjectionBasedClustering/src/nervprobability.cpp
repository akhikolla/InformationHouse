/*Copyright (C) Kristian Nybo, Jarkko Venna
 *
 *This software is released under the GNU Lesser General Public
 *License. See the included file LICENSE for details.*/

#include <limits>
#include "nervprobability.h"
#include "goldensectionsearch.h"
#include "inputprobentropy.h"

namespace dredviz {
NeRVProbability::NeRVProbability (const DistanceMatrix & origData,
                                  double minp):
    prob (origData.getRows (), origData.getRows ()),
    data (origData),
    shift (origData.getMin ()),
    minexp (minp),
    MIN_SIGMA (1E-10)
{
  shift = shift * shift;
}

void
NeRVProbability::update (size_t row, double sigma2)
{

  double sum = 0.0;
  for (size_t col = 0; col < data.getCols (); col++)
  {
    if (col != row)
    {

      double sqdist = (-data (row, col) * data (row, col) + shift) / sigma2;

      if (sqdist > minexp)
        prob (row, col) = exp (sqdist);
      else
        prob (row, col) = exp (minexp);


      sum += prob (row, col);
    }
    else
    {
      prob (row, col) = 0.0;
    }
  }


  for (size_t col = 0; col < getCols (); col++)
  {
    if (col != row)
    {
      prob (row, col) /= sum;
    }
  }
}

void
NeRVProbability::update (const vector < double >&sigma2)
{
  for (size_t row = 0; row < getRows (); ++row)
    update (row, sigma2[row]);
}


double NeRVProbability::findSigma (size_t effectiveNeighbors, size_t index)
{
  DataMatrix sigmaA (1, 1);
  DataMatrix gradient (1, 1);

  GoldenSectionSearch linesearch;
  InputProbEntropy entropy (effectiveNeighbors, 0, (*this));

  gradient (0, 0) = 1.0;

  sigmaA (0, 0) = std::numeric_limits < double >::min ();

  entropy.setRow (index);
  double dummy;
  linesearch (entropy, sigmaA, gradient, 1.0, dummy);

  return sigmaA(0, 0) > MIN_SIGMA ? sigmaA(0,0) : MIN_SIGMA;
}


void
NeRVProbability::findSigma (vector < double >&sigma2,
                            size_t effectiveNeighbors)
{
  DataMatrix sigmaA (1, 1);
  DataMatrix gradient (1, 1);

  GoldenSectionSearch linesearch;
  InputProbEntropy entropy (effectiveNeighbors, 0, (*this));

  gradient (0, 0) = 1.0;

  for (size_t i = 0; i < prob.getRows (); i++)
  {
    sigmaA (0, 0) = std::numeric_limits < double >::min ();

    entropy.setRow (i);
    double dummy;
    linesearch (entropy, sigmaA, gradient, 1.0, dummy);
    sigma2[i] = sigmaA (0, 0) > MIN_SIGMA ? sigmaA(0,0) : MIN_SIGMA;
    //std::cerr <<"sigma2 "<<i<<": "<<sigma2[i]<<std::endl;
  }

}
}
