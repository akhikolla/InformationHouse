/*Copyright (C) Kristian Nybo, Jarkko Venna
 *
 *This software is released under the GNU Lesser General Public
 *License. See the included file LICENSE for details.*/

#include "goldensectionsearch.h"
#include <limits>

namespace dredviz {
GoldenSectionSearch::GoldenSectionSearch ():INITIAL_BRACKETING_MAGNIFICATION (2.0), RATIO (0.61803399),
    TOLERANCE (1.0e-2), MAX_ITER (40), referencePoint (1, 1),
    alpha (0.0), beta (0.0), gamma (0.0), xi (0.0),
    costA (0.0), currentMinimumCost (0.0), costC (0.0), costX (0.0), x (1, 1)
{
}


bool
GoldenSectionSearch::

findInitialBracket (CostFunction & costFunction,
                    DataMatrix & pointA,
                    const DataMatrix & negativeGradient,
                    double initialStepSize)
{

  size_t iter = 0;
  costA = costFunction.evaluate (pointA);

  //  std::cout<<"initial step: "<<initialStepSize<<std::endl;

  for (size_t i = 0; i < pointA.getRows (); i++)
    for (size_t j = 0; j < pointA.getCols (); j++)
      x (i, j) = pointA (i, j) + initialStepSize * negativeGradient (i, j);

  double costB = costFunction.evaluate (x);

  alpha = 0.0;
  beta = initialStepSize;

  if (costA < costB)
  {
    // minimum is somewhere between the two
    gamma = beta;
    currentMinimumCost = costA;
    costC = costB;

    double stepx = initialStepSize / INITIAL_BRACKETING_MAGNIFICATION;

    do
    {
      ++iter;
      for (size_t i = 0; i < x.getRows (); i++)
        for (size_t j = 0; j < x.getCols (); j++)
          x (i, j) = referencePoint (i, j) + stepx * negativeGradient (i, j);
      beta = stepx;
      stepx /= INITIAL_BRACKETING_MAGNIFICATION;
      costX = costFunction.evaluate (x);
    }
    while (costX > costA && iter < MAX_ITER);
    if (costX < currentMinimumCost)
    {
      currentMinimumCost = costX;
    }
    else
    {
      beta = std::numeric_limits < double >::epsilon ();
    }
    costB = costX;

    //           std::cerr<<"lineserach bracketing cost down:"<< beta <<"( " <<costA<<" , "<<costX << " , "<<costB<<" ) diff ax:"<<costA-costX<<" diff bx:"<<costC-costX<<std::endl;

  }
  else
  {
    // find a larger cost beond beta to bracket the minimum

    currentMinimumCost = costB;

    double stepx = initialStepSize * INITIAL_BRACKETING_MAGNIFICATION;

    do
    {
      ++iter;
      for (size_t i = 0; i < x.getRows (); i++)
        for (size_t j = 0; j < x.getCols (); j++)
          x (i, j) = referencePoint (i, j) + stepx * negativeGradient (i, j);
      gamma = stepx;
      stepx *= INITIAL_BRACKETING_MAGNIFICATION;
      costX = costFunction.evaluate (x);
    }
    while (costX < currentMinimumCost && iter < MAX_ITER);

    //            std::cerr<<"linesearch bracketing cost:"<< stepx <<"( " <<costA<<" , "<<currentMinimumCost << " , "<<costX<<" ) diff:"<<costX-currentMinimumCost<<std::endl;


    if (iter == MAX_ITER && costX < costB)
    {
      // If we reach here either the gradient is very close to zero or we
      // take a huge jump towards an undesirable local optima.
      // Only take a small step and try again.
      costC = costB;
      gamma = beta;
      beta = (alpha + gamma) / 2.0;
    }
    else
    {
      costC = costX;
    }
  }
  if (costA < costC && costA < costB)
  {
    // starting position is the best found
    return false;
  }
  else
    return true;
}



double
GoldenSectionSearch::operator ()
(CostFunction & costFunction, DataMatrix & point,
 const DataMatrix & negativeGradient,
 double initialStepSize, double &finalCost)
{

  x = point;
  referencePoint = point;

  if (!findInitialBracket
      (costFunction, point, negativeGradient, initialStepSize))
    return 0.0;
  double bracketLength = gamma;

  size_t iter = 0;


  double tolerance = TOLERANCE;

  // avoid problems if the bracketLength is smaller than TOLERANCE

  if (bracketLength < TOLERANCE)
    tolerance = TOLERANCE * bracketLength;

  while (gamma - alpha > tolerance && iter < MAX_ITER)
  {
    /* Choose a point x from the larger of the segments (a,b) and (b,c). */
    iter++;

    /* Is (b,c) larger than (a,b)? */
    if (gamma - beta > beta - alpha)
    {
      xi = beta + RATIO * (gamma - beta);

      for (size_t i = 0; i < x.getRows (); i++)
        for (size_t j = 0; j < x.getCols (); j++)
          x (i, j) = referencePoint (i, j) + xi * negativeGradient (i, j);

      costX = costFunction.evaluate (x);

      /* If f(x) < f(b), the new bracket is (b, x, c), otherwise
         it's (a, b, x) */

      if (costX < currentMinimumCost)
      {
        alpha = beta;
        costA = currentMinimumCost;

        beta = xi;
        currentMinimumCost = costX;
      }
      else
      {
        gamma = xi;
        costC = costX;
      }
    }
    else
    {
      xi = alpha + RATIO * (beta - alpha);

      for (size_t i = 0; i < x.getRows (); i++)
        for (size_t j = 0; j < x.getCols (); j++)
          x (i, j) = referencePoint (i, j) + xi * negativeGradient (i, j);

      costX = costFunction.evaluate (x);

      /* If f(x) < f(b), the new bracket is (a, x, b), otherwise
         it's (x, b, c) */

      if (costX < currentMinimumCost)
      {
        gamma = beta;
        costC = currentMinimumCost;

        beta = xi;
        currentMinimumCost = costX;
      }
      else
      {
        alpha = xi;
        costA = costC;
      }
    }
  }

  /* Let B be the midpoint, and return whichever of A, B and C has the lowest
     cost. */


  beta = (alpha + gamma) / 2.0;

  double nextInitialStep = beta;

  for (size_t i = 0; i < point.getRows (); i++)
    for (size_t j = 0; j < point.getCols (); j++)
      point (i, j) = referencePoint (i, j) + beta * negativeGradient (i, j);


  currentMinimumCost = costFunction.evaluate (point);



  // verify that the middle point has a cost smaller than the ends
  if (costA < costC)
  {
    if (currentMinimumCost > costA)
    {
      finalCost = costA;
      {
        for (size_t i = 0; i < point.getRows (); i++)
          for (size_t j = 0; j < point.getCols (); j++)
            point (i, j) =
              referencePoint (i, j) + alpha * negativeGradient (i, j);

        nextInitialStep = alpha + TOLERANCE;
      }
    }
    else
    {
      finalCost = currentMinimumCost;
      //      std:cerr<<"lineserach end cost:  " <<currentMinimumCost<<"\n";
    }
  }
  else
  {
    if (currentMinimumCost > costC)
    {
      finalCost = costC;
      for (size_t i = 0; i < point.getRows (); i++)
        for (size_t j = 0; j < point.getCols (); j++)
          point (i, j) =
            referencePoint (i, j) + gamma * negativeGradient (i, j);

      nextInitialStep = gamma;
    }
    else
    {
      finalCost = currentMinimumCost;
    }
  }

  return nextInitialStep;
}
}
