#include <algorithm>

#include "rankmatrix.h"

#include "distcomp.h"

namespace dredviz {

RankMatrix::RankMatrix(const DistanceMatrix& dist)
    : DistanceMatrix(dist.getRows())
{
  CompareIndicesDist cmp(dist);
  calculateRanks(cmp, dist);
}


RankMatrix::RankMatrix(const DistanceMatrix& dist, RankMatrix& reverseRanks)
    : DistanceMatrix(dist.getRows())
{
  CompareIndicesDist cmp(dist);
  calculateRanks(cmp, dist, reverseRanks);
}


RankMatrix::RankMatrix(const DistanceMatrix& dist,
                       const RankMatrix& inputRanks, bool assumeBestCase)
    : DistanceMatrix(dist.getRows())
{
  CompareIndicesUsingRankMatrix cmp(dist, inputRanks, assumeBestCase);
  calculateRanks(cmp, dist);
}


void RankMatrix::calculateRanks(CompareIndicesDist& cmpIndices,
                                const DistanceMatrix& dist)
{
  const size_t N = dist.getRows();

  /* Once sorted, indicesSortedByRank[i] == the index of the point with rank
     i */

  std::vector < size_t > indicesSortedByRank(N - 1);

  for (size_t i = 0; i < N; i++)
  {
    cmpIndices.setIndex(i);

    size_t ind = 0;

    for (size_t j = 0; j < N; ++j)
    {
      if (j != i)
      {
        indicesSortedByRank[ind] = j;
        ++ind;
      }
    }

    sort (indicesSortedByRank.begin(), indicesSortedByRank.end(), cmpIndices);

    for (size_t j = 0; j < indicesSortedByRank.size(); ++j)
    {
      (*this) (i, indicesSortedByRank[j]) = j + 1;
    }


    (*this) (i, i) = 0;
  }
}


void RankMatrix::calculateRanks(CompareIndicesDist& cmpIndices,
                                const DistanceMatrix& dist, RankMatrix& reverseRanks)
{
  const size_t N = dist.getRows();

  reverseRanks = RankMatrix(N); // Make sure reverseRanks is N*N

  /* Once sorted, indicesSortedByRank[i] == the index of the point with rank
     i */

  std::vector < size_t > indicesSortedByRank(N - 1);

  for (size_t i = 0; i < N; i++)
  {
    cmpIndices.setIndex(i);

    size_t ind = 0;

    for (size_t j = 0; j < N; ++j)
    {
      if (j != i)
      {
        indicesSortedByRank[ind] = j;
        ++ind;
      }
    }

    sort (indicesSortedByRank.begin(), indicesSortedByRank.end(), cmpIndices);

    for (size_t j = 0; j < indicesSortedByRank.size(); ++j)
    {
      (*this) (i, indicesSortedByRank[j]) = j + 1;

      reverseRanks(i, indicesSortedByRank[j]) = N - 1 - j;
    }


    (*this) (i, i) = 0;

    reverseRanks(i, i) = 0;
  }
}


void RankMatrix::calculateRanks(CompareIndicesUsingRankMatrix& cmpIndices,
                                const DistanceMatrix& dist)
{
  const size_t N = dist.getRows();

  /* Once sorted, indicesSortedByRank[i] == the index of the point with rank
     i */

  std::vector < size_t > indicesSortedByRank(N - 1);

  for (size_t i = 0; i < N; i++)
  {
    cmpIndices.setIndex(i);

    size_t ind = 0;

    for (size_t j = 0; j < N; ++j)
    {
      if (j != i)
      {
        indicesSortedByRank[ind] = j;
        ++ind;
      }
    }

    sort (indicesSortedByRank.begin(), indicesSortedByRank.end(), cmpIndices);

    for (size_t j = 0; j < indicesSortedByRank.size(); ++j)
    {
      (*this) (i, indicesSortedByRank[j]) = j + 1;
    }


    (*this) (i, i) = 0;
  }
}
}
