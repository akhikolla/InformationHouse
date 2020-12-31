#ifndef RANKMATRIX_HH
#define RANKMATRIX_HH

#include "distancematrix.h"
#include "distcomp.h"

namespace dredviz {
class CompareIndicesUsingRankMatrix;

class RankMatrix : public DistanceMatrix
{
private:
  void calculateRanks(CompareIndicesDist& cmp, const DistanceMatrix& dist);

  void calculateRanks(CompareIndicesDist& cmp, const DistanceMatrix& dist,
                      RankMatrix& reverseRanks);

  void calculateRanks(CompareIndicesUsingRankMatrix& cmp,
                      const DistanceMatrix& dist);

public:

  /* Creates a matrix r such that r(i,j) is the rank of point j when all
     the points have been sorted according to their distance (as reported
     by dist) from point i. For example, if j is the point second closest to
     point i, r(i,j) = 2.

     Note that this constructor assumes that all the ranks are uniquely
     determined, ie., that for each point i, no two points j and k are at
     exactly the same distance from it. If this is not the case, use
     RankMatrix(DistanceMatrix&, RankMatrix&, ...) instead. */

  RankMatrix(const DistanceMatrix& dist);


  /* Otherwise identical to the previous constructor, but stores the ranks
   * in the constructed RankMatrix in reverse order in reverseRanks. That is,
   * if r(i,j) = 1 for the constructed RankMatrix, reverseRanks(i,j) = N-1,
   * where N = dist.getRows(). (The reverse ranks are used to calculate the
   * maximal error in klrank.cc) */

  RankMatrix(const DistanceMatrix& dist, RankMatrix& reverseRanks);


  /* When the rank r(i,j) is not uniquely determined because points j and k
     are both at the same distance from point i, this constructor uses
     another RankMatrix, inputRanks, to decide the tie.
     If assumeBestCase == true, j and k will be ranked in the same order as
     they are ranked in inputRanks; otherwise the order will be opposite. */

  RankMatrix(const DistanceMatrix& dist, const RankMatrix& inputRanks,
             bool assumeBestCase);


  /* Create an uninitialized rank matrix of dimension dim*dim. */

  RankMatrix(size_t dim = 1) : DistanceMatrix(dim) { }
};


class CompareIndicesUsingRankMatrix
{
private:
  const DistanceMatrix & dist;
  const RankMatrix& origRanks;
  size_t index;
  bool best;
public:
  CompareIndicesUsingRankMatrix(const DistanceMatrix & distdat,
                                const RankMatrix& origRanks, bool best):
      dist(distdat),
      origRanks(origRanks),
      index(0),
      best(best)
  {};

  void setIndex(size_t in)
  {
    index=in;
  };
  void setBest(bool dir)
  {
    best=dir;
  };

  bool operator () (size_t j1, size_t j2)
  {
    // if(j1==j2)
    //  return true;
    if (dist(index, j1) == dist(index, j2))
    {
      //	std::cerr << index<<best<<" equal dist:("<<j1<<","<<j2<<"): "<<origRanks[j1]<<" "<<origRanks[j2]<< " "<<(origRanks[j1] < origRanks[j2])<<std::endl;
      if (best)
        return (origRanks(index,j1) < origRanks(index,j2));
      else
        return (origRanks(index,j1) > origRanks(index,j2));
    }
    //   std::cerr << index<<best<<" dist:("<<j1<<","<<j2<<"): "<<dist(index,j1)<<" "<<dist(index,j2)<< " "<<(dist(index,j1) < dist(index,j2))<<std::endl;
    return (dist(index,j1)<dist(index,j2));
  };

};
}
#endif //RANKMATRIX_HH
