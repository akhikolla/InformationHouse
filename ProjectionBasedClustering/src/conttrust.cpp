/*Copyright (C) Kristian Nybo, Jarkko Venna
 *
 *This software is released under the GNU Lesser General Public
 *License. See the included file LICENSE for details.*/

#include "conttrust.h"
#include "datamatrix.h"
#include "distcomp.h"
#include <vector>
#include <algorithm>
#include <Rcpp.h>

using Rcpp::NumericMatrix;

namespace dredviz {
ContTrust::ContTrust (size_t maxNeighborhood):
    maxNeighborhood (maxNeighborhood), measures (maxNeighborhood, 6)
    //precrec(maxNeighborhood, 2)
{
  for (size_t i = 0; i < maxNeighborhood; ++i)
  {
    measures (i, 0) = -1.0;
    measures (i, 1) = -1.0;
    measures (i, 2) = -1.0;
    measures (i, 3) = -1.0;
    //precrec(i,0) = precrec(i,1) = 0;
  }
}


NumericMatrix
ContTrust::measure (const DistanceMatrix & origDist,
                    const DistanceMatrix & projDist)
{

  //  double trustErrBest, contErrBest, trustErrWorst, contErrWorst, normalizer;
  //  trustErrBest = contErrBest = trustErrWorst = contErrWorst = normalizer =
  //  0.0;

  /*projDist and origDist are both square matrices, so projDist.getRows(),
     projDist.getCols(), origDist.getRows() and origDist.getCols() are all the
     same. */

  const size_t N = projDist.getRows ();

  if (projDist.getRows () != origDist.getRows ())
  {
    std::string errMsg = "The original data and the projected data must have"
                         " the same number of points.";
    throw (RowMismatch (errMsg));
  }


  for (size_t i = 0; i < maxNeighborhood; ++i)
  {
    measures (i, 0) = 0.0;      //trust best
    measures (i, 1) = 0.0;      //trust worst
    measures (i, 2) = 0.0;      //cont best
    measures (i, 3) = 0.0;      //cont worst
  }

  std::vector < size_t > origRank (N, 0);
  std::vector < size_t > projRank (N, 0);


  CompareIndicesDist origCmp (origDist);
  CompareIndicesProjDist projCmp (projDist, origRank, true);



  /* A temporary vector. Once sorted, indicesSortedByRank[i] == the index of the point with rank i */

  std::vector < size_t > indicesSortedByRank (N - 1);

  /* Calculate continuity and trustworthiness. With 'best case' comparison
     points in the output (input) space that are equally distant from point
     i are sorted so that they are in the same order as they were in the input
     (output) space. */


  for (size_t i = 0; i < N; i++)
  {
    origCmp.setIndex (i);
    // fill the vector. Leave out the index of the current data point
    //to avoid messing out the sort
    size_t ind = 0;
    for (size_t j = 0; j < N; ++j)
    {
      if (j != i)
      {
        indicesSortedByRank[ind] = j;
        ++ind;
      }
    }
    sort (indicesSortedByRank.begin (), indicesSortedByRank.end (), origCmp);
    // fill the ranks in the order of the data points
    for (size_t j = 0; j < indicesSortedByRank.size (); ++j)
    {
      origRank[indicesSortedByRank[j]] = j + 1;

    }

    origRank[i] = 0;
    /*
       std::cerr << "orig sorted "<<i<<":";
       for(size_t j=0;j<N-1;++j)
       std::cerr <<indicesSortedByRank[j]<<" ";
       std::cerr << std::endl;
     */

    // Best case sorting
    projCmp.setIndex (i);
    projCmp.setBest (true);

    sort (indicesSortedByRank.begin (), indicesSortedByRank.end (), projCmp);
    // rearrange the ranks in the order of the data points
    for (size_t j = 0; j < indicesSortedByRank.size (); ++j)
    {
      projRank[indicesSortedByRank[j]] = j + 1;
    }
    projRank[i] = 0;

    /*      std::cerr << "proj sorted best "<<i<<":";
       for(size_t j=0;j<N-1;++j)
       std::cerr <<indicesSortedByRank[j]<<" ";
       std::cerr << std::endl;
     */
    for (size_t j = 0; j < N; j++)
    {
      if (i != j)
      {

        for (size_t neighborhoodSize = 1;
             neighborhoodSize <= maxNeighborhood; neighborhoodSize++)
        {
          if (projRank[j] <= neighborhoodSize
              && origRank[j] > neighborhoodSize)
          {
            measures (neighborhoodSize - 1, 0) += origRank[j] -
                                                  neighborhoodSize;
            //precrec(neighborhoodSize - 1, 0) += 1; //one false positive
          }

          if (origRank[j] <= neighborhoodSize
              && projRank[j] > neighborhoodSize)
          {
            measures (neighborhoodSize - 1, 2) += projRank[j] -
                                                  neighborhoodSize;
            //precrec(neighborhoodSize - 1, 1) += 1; //one miss
          }

        }
      }

    }
    //Worst case sorting

    projCmp.setBest (false);

    sort (indicesSortedByRank.begin (), indicesSortedByRank.end (), projCmp);
    // rearrange the ranks in the order of the data points
    for (size_t j = 0; j < indicesSortedByRank.size (); ++j)
    {
      projRank[indicesSortedByRank[j]] = j + 1;
    }
    projRank[i] = 0;
    /*
       std::cerr << "proj sorted worst"<<i<<":";
       for(size_t j=0;j<N-1;++j)
       std::cerr <<indicesSortedByRank[j]<<" ";
       std::cerr << std::endl;
     */
    for (size_t j = 0; j < N; j++)
    {
      if (i != j)
      {

        for (size_t neighborhoodSize = 1;
             neighborhoodSize <= maxNeighborhood; neighborhoodSize++)
        {
          if (projRank[j] <= neighborhoodSize
              && origRank[j] > neighborhoodSize)
          {
            measures (neighborhoodSize - 1, 1) += origRank[j] -
                                                  neighborhoodSize;
          }

          if (origRank[j] <= neighborhoodSize
              && projRank[j] > neighborhoodSize)
          {
            measures (neighborhoodSize - 1, 3) += projRank[j] -
                                                  neighborhoodSize;
          }

        }
      }

    }

  }
  NumericMatrix output(maxNeighborhood, 7);
  //ostr << "Neighborhood size; worst-case trustworthiness; average trustworthiness; best-case trustworthiness; worst-case continuity; average continuity; best-case continuity\n";
  for (size_t neighborhoodSize = 1;
       neighborhoodSize <= maxNeighborhood; neighborhoodSize++)
  {
    double normalizer = -1.0;

    if (neighborhoodSize < (double) N / 2.0)
    {
      normalizer = 2.0 / (N * neighborhoodSize *
                          (2.0 * N - 3.0 * neighborhoodSize - 1.0));
    }
    else
    {
      normalizer =
        2.0 / (N * (N - neighborhoodSize) * (N - neighborhoodSize - 1.0));
    }

    measures (neighborhoodSize - 1, 0) *= normalizer;
    measures (neighborhoodSize - 1, 1) *= normalizer;
    measures (neighborhoodSize - 1, 2) *= normalizer;
    measures (neighborhoodSize - 1, 3) *= normalizer;

    output(neighborhoodSize - 1, 0) = neighborhoodSize;
    output(neighborhoodSize - 1, 1) = measures(neighborhoodSize - 1, 1);
    output(neighborhoodSize - 1, 2) = (measures(neighborhoodSize - 1, 0) + measures(neighborhoodSize - 1, 1)) / 2.0;
    output(neighborhoodSize - 1, 3) = measures(neighborhoodSize - 1, 0);
    output(neighborhoodSize - 1, 4) = measures(neighborhoodSize - 1, 3);
    output(neighborhoodSize - 1, 5) = (measures(neighborhoodSize - 1, 2) + measures(neighborhoodSize - 1, 3)) / 2.0;
    output(neighborhoodSize - 1, 6) = measures(neighborhoodSize - 1, 2);
    /*
    ostr << neighborhoodSize << " " 
      << measures (neighborhoodSize - 1, 1) << "  " //worst TW
      << (measures (neighborhoodSize - 1, 0) +
         measures (neighborhoodSize - 1, 1)) / 2.0 
      << "  " << measures (neighborhoodSize - 1, 0) << "  "  //best TW
      << measures (neighborhoodSize - 1, 3)  //worst CONT
      << "  " << (measures (neighborhoodSize - 1, 2) + 
          measures (neighborhoodSize - 1, 3)) / 2.0 
      << "  " << measures (neighborhoodSize - 1, 2) << "  "    //best CONT
      << std::endl;
      */
  }
  return output;
}

double
ContTrust::getContinuity (size_t neighborhoodSize, int wcase)
{
  if (wcase == 0)
    return 1.0 - (measures (neighborhoodSize - 1, 2) +
                  measures (neighborhoodSize - 1, 3)) / 2.0;
  if (wcase == 1)
    return 1.0 - measures (neighborhoodSize - 1, 2);
  if (wcase == 2)
    return 1.0 - measures (neighborhoodSize - 1, 3);
  return -1.0;
}


double
ContTrust::getTrustworthiness (size_t neighborhoodSize, int wcase)
{
  if (wcase == 0)
    return 1.0 - (measures (neighborhoodSize - 1, 0) +
                  measures (neighborhoodSize - 1, 1)) / 2.0;
  if (wcase == 1)
    return 1.0 - measures (neighborhoodSize - 1, 0);
  if (wcase == 2)
    return 1.0 - measures (neighborhoodSize - 1, 1);
  return -1.0;
}
}
