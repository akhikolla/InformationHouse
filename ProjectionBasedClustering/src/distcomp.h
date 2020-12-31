/*Copyright (C) Kristian Nybo, Jarkko Venna
 *
 *This software is released under the GNU Lesser General Public
 *License. See the included file LICENSE for details.*/
#ifndef DISTCOMP_HH
#define DISTCOMP_HH

#include "distancematrix.h"

namespace dredviz {
class CompareIndicesDist
{
private:
  const DistanceMatrix & dist;
  size_t index;
public:
  CompareIndicesDist(const DistanceMatrix & distdat):
      dist(distdat),
      index(0){};

  void setIndex(size_t in)
  {
    index=in;
  };

  bool operator () (size_t j1, size_t j2)
  {
    return (dist(index,j1)<dist(index,j2));
  };

};



class CompareIndicesProjDist
{
private:
  const DistanceMatrix & dist;
  const std::vector<size_t> & origRanks;
  size_t index;
  bool best;
public:
  CompareIndicesProjDist(const DistanceMatrix & distdat,
                         const std::vector<size_t> & origRanks, bool best):
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
        return (origRanks[j1] < origRanks[j2]);
      else
        return (origRanks[j1] > origRanks[j2]);
    }
    //   std::cerr << index<<best<<" dist:("<<j1<<","<<j2<<"): "<<dist(index,j1)<<" "<<dist(index,j2)<< " "<<(dist(index,j1) < dist(index,j2))<<std::endl;
    return (dist(index,j1)<dist(index,j2));
  };

};


}



#endif //DISTCOMP_HH
