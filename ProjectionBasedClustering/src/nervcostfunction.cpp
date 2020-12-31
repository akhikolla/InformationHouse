/*Copyright (C) Kristian Nybo, Jarkko Venna
 *
 *This software is released under the GNU Lesser General Public
 *License. See the included file LICENSE for details.*/

#include "nervcostfunction.h"

#include "euclidean.h"
#include "euclideansquared.h"

#include <string>

#include <cmath>
#include <limits>

#include <Rcpp.h>
#include <sstream>

namespace dredviz {
void
NeRVCostFunction::calculateFinalNeighborhoods (size_t effectiveNeighbors)
{
  inputProb.findSigma(finalNeighborhoods, effectiveNeighbors);
  for (size_t i = 0; i < sigmaSqrd.size(); i++)
  {
    if (2.0 * radius.value() * radius.value() > finalNeighborhoods[i])
    {
      sigmaSqrd[i] = 2.0 * radius.value() * radius.value();
    }
    else
    {
      sigmaSqrd[i] = finalNeighborhoods[i];
    }
  }
}



void
NeRVCostFunction::updateOutputProb (const DataMatrix & projData)
{
  double minexp = -200.0;
  //  CalculateDistanceMatrix calc;
  //  Euclidean metric;

//  calc(projData, metric, projDist);

// outputProb.update(sigmaSqrd);

  for (size_t i = 0; i < projData.getRows(); ++i)
  {
    double sum = 0.0;
    for (size_t j = 0; j < projData.getRows(); ++j)
    {
      if (i != j)
      {
        double sqdist = 0.0;

        for (size_t d = 0; d < projData.getCols(); ++d)
        {
          double diff = (projData(i, d) - projData(j, d));
          sqdist += diff * diff;
        }
        sqdist = (-sqdist + minimumDistance) / sigmaSqrd[i];

        if (sqdist > minexp)
          outputProb(i, j) = exp(sqdist);
        else
          outputProb(i, j) = exp(minexp);

        sum += outputProb(i, j);
      }
      else
      {
        outputProb(i, j) = 0.0;
      }
    }

    for (size_t d = 0; d < outputProb.getCols(); ++d)
    {
      outputProb(i, d) /= sum;
    }
  }
}

void
NeRVCostFunction::updateDynamicParameters(size_t currentRound,
    size_t totalRounds,
    const DataMatrix & projData)
{
  radius.update(currentRound, totalRounds);

  for (size_t i = 0; i < sigmaSqrd.size(); i++)
  {
    if (2.0 * radius.value() * radius.value() > finalNeighborhoods[i])
    {
      sigmaSqrd[i] = 2.0 * radius.value() * radius.value();
    }
    else
    {
      sigmaSqrd[i] = finalNeighborhoods[i];
    }
    //      std::cerr<<"sigma update, rad:"<<radius.value()<<" final:"<<finalNeighborhoods[i]<<" Sigma:"<<sigmaSqrd[i]<<std::endl;
  }

  updateMinimumDistance(projData);
  inputProb.update(sigmaSqrd);
}


double
NeRVCostFunction::evaluate(const DataMatrix & projData)
{
  updateOutputProb(projData);

  double cost = 0.0;


  for (size_t row = 0; row < inputProb.getRows(); row++)
  {
    for (size_t col = 0; col < inputProb.getCols(); col++)
    {
      if (row != col)
      {
        //              std::cerr << "p: "<< inputProb(row,col)<<" q: "<<outputProb(row,col)<< std::endl;
        cost += weights[row] * lambda * inputProb(row, col) *
                (std::log(inputProb(row, col)) - std::log(outputProb(row, col)));

        cost += weights[row] * (1.0 - lambda) * outputProb(row, col) *
                (std::log(outputProb(row, col)) - std::log(inputProb(row, col)));
      }
    }
  }
  //  feedback << "Current cost: " << cost << std::endl;
  //feedback.flush();

  return cost / projData.getRows();
}


void
NeRVCostFunction::updateDataRepresentation(const DataMatrix & newData)
{
  updateOutputProb(newData);
}



double
NeRVCostFunction::getGradient(const DataMatrix & projData,
                              DataMatrix & gradient)
{

  double sqsum=0.0;
  size_t n=gradient.getRows();
  size_t dim=gradient.getCols();

  updateOutputProb (projData);

  for (size_t i = 0; i < n; ++i)
  {
    w[i] = 0.0;
    for (size_t j = 0; j < n; ++j)
    {
      if (i != j)
      {
        dDval(i, j) =
          (1.0 - lambda) * (std::log (outputProb(i, j)) - std::log (inputProb(i, j)) +
                            1.0)
          -lambda * inputProb(i, j) / outputProb(i, j);

        w[i] += dDval(i,j) * outputProb(i,j) / sigmaSqrd[i];
      }
    }
  }

  for (size_t g_row = 0; g_row < n; g_row++)
  {
    for (size_t g_col = 0; g_col < dim; g_col++)
      gradient(g_row,g_col) = 0.0;

    for (size_t g_col = 0; g_col < dim; g_col++)
    {
      //first term
      double tmp1 = 0.0;
      double tmp2 = 0.0;

      for (size_t j = 0; j < n; j++)
      {
        if(j != g_row)
        {
          tmp1 += dDval(g_row,j) * outputProb(g_row,j) / sigmaSqrd[g_row];
        }
      }

      for (size_t k = 0; k < n; k++)
      {
        if(k != g_row)
        {
          tmp2 += outputProb(g_row,k)
            * (projData(g_row,g_col) - projData(k,g_col));
        }
      }

      gradient(g_row,g_col) += tmp1 * tmp2;

      //second term

      for (size_t i = 0; i < n; i++)
      {
        if(i != g_row)
        {
          gradient(g_row,g_col) += w[i] * outputProb(i,g_row)
            * (projData(g_row,g_col) - projData(i, g_col));
        }
      }

      //third term

      for(size_t j = 0; j < n; j++)
      {
        if(j != g_row)
        {
          gradient(g_row, g_col)
            -= dDval(g_row, j) * outputProb(g_row, j)
            * (projData(g_row, g_col) - projData(j, g_col))
            / sigmaSqrd[g_row];
        }
      }

      //fourth term

      for(size_t i = 0; i < n; i++)
      {
        if(i != g_row)
        {
          gradient(g_row,g_col)
            -= dDval(i, g_row) * outputProb(i,g_row)
            * (projData(g_row, g_col) - projData(i, g_col))
            / sigmaSqrd[i];
        }
      }
      gradient(g_row, g_col) *= weights[g_row];
    }
  }

  /*
  DataMatrix jgradient(gradient.getRows(), gradient.getCols());
  for (size_t t = 0; t < n; ++t)
  {
    for (size_t d = 0; d < dim; ++d)
      jgradient(t, d) = 0.0;

    for (size_t j = 0; j < jgradient.getRows(); ++j)
    {
      if (j != t)
      {

        for (size_t d = 0; d < dim; ++d)
        {
          jgradient(t, d) -=
            dDval(j, t) * outputProb(j,t) * (outputProb(j,t) - 1.0) *
            (projData(j,d) - projData(t,d)) / sigmaSqrd[j];

        }


        for (size_t d = 0; d < dim; ++d)
        {
          jgradient(t, d) -= dDval(t, j) * outputProb(t, j) *
                            (projData(t, d) - projData(j, d)) / sigmaSqrd[t];
        }
        for (size_t k = 0; k < n; ++k)
        {
          if (k != t)
          {
            for (size_t d = 0; d < dim; ++d)
            {
              jgradient(t, d) +=
                dDval(t, j) * outputProb(t, k) * outputProb(t, j) *
                (projData(t, d) - projData(k, d)) / sigmaSqrd[t];
            }
            if (k != j)
            {
              for (size_t d = 0; d < dim; ++d)
              {
                jgradient(t, d) -=
                  dDval(j, k) * outputProb(j, k) * outputProb(j, t) *
                  (projData(j, d) - projData(t, d)) / sigmaSqrd[j];
              }
            }

          }
        }
      }

    }
  }

  double diff = 0;
  */

  for (size_t t = 0; t < n; t++)
  {
    for (size_t d = 0; d < dim; ++d)
    {
      sqsum += (gradient (t, d) * gradient(t,d));
      /*diff += (gradient(t,d) - jgradient(t,d))
        * (gradient(t,d) - jgradient(t,d));*/
    }
  }

  Rcpp::Rcout << "gradient " << sqsum << std::endl;

  return sqsum;
}

/*
double
NeRVCostFunction::getGradient(const DataMatrix & projData,
                              DataMatrix & gradient)
{

  double sqsum=0.0;
  size_t n=gradient.getRows();
  size_t dim=gradient.getCols();

  updateOutputProb (projData);

  for (size_t i = 0; i < n; ++i)
  {
    for (size_t j = 0; j < n; ++j)
    {
      if (i != j)
      {
        dDval(i, j) =
          (1.0 - lambda) * (std::log (outputProb(i, j)) - std::log (inputProb(i, j)) +
                            1.0) - lambda * inputProb(i, j) / outputProb(i, j);
      }
    }
  }

  for (size_t t = 0; t < n; ++t)
  {
    for (size_t d = 0; d < dim; ++d)
      gradient(t, d) = 0.0;

    for (size_t j = 0; j < gradient.getRows(); ++j)
    {
      if (j != t)
      {

        for (size_t d = 0; d < dim; ++d)
        {
          gradient(t, d) -=
            dDval(j, t) * outputProb(j,t) * (outputProb(j,t) - 1.0) *
            (projData(j,d) - projData(t,d)) / sigmaSqrd[j];

        }


        for (size_t d = 0; d < dim; ++d)
        {
          gradient(t, d) -= dDval(t, j) * outputProb(t, j) *
                            (projData(t, d) - projData(j, d)) / sigmaSqrd[t];
        }
        for (size_t k = 0; k < n; ++k)
        {
          if (k != t)
          {
            for (size_t d = 0; d < dim; ++d)
            {
              gradient(t, d) +=
                dDval(t, j) * outputProb(t, k) * outputProb(t, j) *
                (projData(t, d) - projData(k, d)) / sigmaSqrd[t];
            }
            if (k != j)
            {
              for (size_t d = 0; d < dim; ++d)
              {
                gradient(t, d) -=
                  dDval(j, k) * outputProb(j, k) * outputProb(j, t) *
                  (projData(j, d) - projData(t, d)) / sigmaSqrd[j];
              }
            }

          }
        }
      }

    }
    for (size_t d = 0; d < dim; ++d)
      sqsum += gradient (t, d) * gradient (t, d);
  }
  return sqsum;
}
*/

NeRVCostFunction::NeRVCostFunction(const DistanceMatrix & origDist,
                                   DataMatrix & projData,
                                   LineSearch & linesearch,
                                   DynamicDouble radius,
                                   double lambda, size_t neighborhoodSize,
                                   const std::vector<double>& weights,
                                   std::ostream & feedback):
    origDist(origDist),
    minexp(-450),
    radius(radius),
    lambda(lambda),
    finalNeighborhoods(origDist.getRows(),0),
    sigmaSqrd(origDist.getRows()),
    linesearch(linesearch),
    inputProb(origDist),
    outputProb(projData.getRows(),projData.getRows()),
    weights(weights),
    feedback(feedback),
    minimumDistance(origDist.getMax()),
    origMinimumDistance(origDist.getMin()),
    dDval(projData.getRows(),projData.getRows()),
    w(projData.getRows())

{
  calculateFinalNeighborhoods(neighborhoodSize);
  inputProb.update(sigmaSqrd);

  //Printing out probabilities for inspection
  //std::string file = "nervprob.prb";
  //SOMPackExporter ex(file);
  //ex.exportData(inputProb);

  updateMinimumDistance(projData);
}


void
NeRVCostFunction::updateMinimumDistance(const DataMatrix & projData)
{
  EuclideanSquared distSqrd;

  minimumDistance = std::numeric_limits < double >::max ();

  double distance = 0.0;

  for (size_t i = 0; i < projData.getRows (); i++)
  {
    for (size_t j = 0; j < i; j++)
    {
      distance = distSqrd(projData, i, j);

      if (distance < minimumDistance)
        minimumDistance = distance;
    }
  }
}


void NeRVCostFunction::reportParameters(std::string& target)
{
  std::ostringstream params;
  params << "Lambda: " << lambda << "\nCurrent radius: " << radius.value()
    << "\n";

  target = params.str();
}
}
