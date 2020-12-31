
#include "auxFun.h"

/*************
* Fast compute multiscale statistic via quasiconvex funtion on sum of two sets
* Input:
*   data:    numeric vector specifying sample points
*   type:    integer specifying penalized normal (0), Poisson (1)
* Output:
*   vms:     value of multiscale statistic
* Note:
*   A tailorized implementation of 'convHullConstrSum' with
*      stIndexP = 0
*      vSizeP   = n
*      sizeP    = n
*      sizeQ    = n
*      JQ: n-1, ... 0
*      setP = {(i, csi) | i = 1, ..., n}
*      setQ = {(-i,-csi) | i = (n-1), ..., 0}
****************/
// [[Rcpp::export(".mStat")]]
double mStat(NumericVector data, int type = 0) {
  int n = data.size();
  NumericVector cs(n+1); // cumulative sum
  cs[0] = 0.;
  for (int i = 0; i < n; ++i)
    cs[i+1] = cs[i] + data[i];
  double vms = -INFINITY; // value of multiscale statistics
  double (*vs) (double,int,int); // compute value of local statistics
  
  switch (type) {
  case 0: // penalized normal case
    vs = pNorm;
    break;
  case 1: // Poisson
    vs = Pois;
    break;
  default:
    REprintf("Unsupported type %d!", type);
  break;
  }
  
  // Find neighbours of each pi in convHull(pi,...,pn) for valid part of P
  // via Graham scan algorithm
  //      Upper neighbour point
  IntegerVector indexConvP(n);   // index for (upper/lower) convex hull of P
  IntegerVector indexUpNbP(n-1); // index for upper neighbour point of pi in convHull(pi,...,pn) for i < n
  int vSizeCP = 1;                    // valid size of 'indexConvP'
  indexConvP[0] = n-1;            // the last point is always in the convex hull
  for (int i = n-2; i >= 0; --i) {
    while (vSizeCP > 1 &&
           acw(i+1, cs[i+1],
               indexConvP[vSizeCP-1]+1, cs[indexConvP[vSizeCP-1]+1],
                                          indexConvP[vSizeCP-2]+1, cs[indexConvP[vSizeCP-2]+1])) {
      --vSizeCP;
    }
    indexConvP[vSizeCP++]  = i;
    indexUpNbP[i] = indexConvP[vSizeCP-2];
  }
  //      Lower neighbour point
  IntegerVector indexLoNbP(n-1);
  vSizeCP = 1; // note that indexConvP[0] == n-1;
  for (int i = n-2; i >= 0; --i) {
    while (vSizeCP > 1 &&
           cw(i+1, cs[i+1],
              indexConvP[vSizeCP-1]+1, cs[indexConvP[vSizeCP-1]+1],
                                         indexConvP[vSizeCP-2]+1, cs[indexConvP[vSizeCP-2]+1]))
      --vSizeCP;
    indexConvP[vSizeCP++]  = i;
    indexLoNbP[i] = indexConvP[vSizeCP-2];
  }
  
  // Iteratively update R
  IntegerVector indexUpConvQ(n);   // index for upper convex hull of Q
  int vSizeUCQ = 0;                    // valid size of 'indexUpConvQ'
  IntegerVector indexLoConvQ(n);   // index for lower convex hull of Q
  int vSizeLCQ = 0;                    // valid size of 'indexLoConvQ'
  int auxSizeUKQ, auxSizeLKQ;          // # points in K in upper and lower parts severally
  int auxIndexQ = n-1;
  int auxIndexQcs;
  for (int vIndexP = 0; vIndexP < n; ++vIndexP) {
    //      Find convex hull of Q
    for (int i = auxIndexQ; i >= n-vIndexP-1; --i) {
      while (vSizeUCQ > 1 &&             // along upper side
             acw(i+1-n, -cs[n-i-1],
                 indexUpConvQ[vSizeUCQ-1]+1-n, -cs[n-indexUpConvQ[vSizeUCQ-1]-1],
                                                  indexUpConvQ[vSizeUCQ-2]+1-n, -cs[n-indexUpConvQ[vSizeUCQ-2]-1]))
        --vSizeUCQ;
      indexUpConvQ[vSizeUCQ++] = i;
      while (vSizeLCQ > 1 &&              // along lower side
             cw(i+1-n, -cs[n-i-1],
                indexLoConvQ[vSizeLCQ-1]+1-n, -cs[n-indexLoConvQ[vSizeLCQ-1]-1],
                                                 indexLoConvQ[vSizeLCQ-2]+1-n, -cs[n-indexLoConvQ[vSizeLCQ-2]-1]))
        --vSizeLCQ;
      indexLoConvQ[vSizeLCQ++] = i;
    }
    auxIndexQ = n-vIndexP-2;
    //      Add points to 'R'
    auxIndexQcs = n-indexUpConvQ[vSizeUCQ-1]-1;
    vms = std::max(vms, vs(cs[vIndexP+1]-cs[auxIndexQcs],vIndexP+1-auxIndexQcs,n));
    
    if (indexLoConvQ[vSizeLCQ-1] != indexUpConvQ[vSizeUCQ-1]) {
      auxIndexQcs = n-indexLoConvQ[vSizeLCQ-1]-1;
      vms = std::max(vms, vs(cs[vIndexP+1]-cs[auxIndexQcs],vIndexP+1-auxIndexQcs,n));
    }
    if (vIndexP == n-1) { // last point in P: include all remaining points in Q
      for (int i = vSizeUCQ-1; i > 0; --i) { // along upper side
        auxIndexQcs = n-indexUpConvQ[i-1]-1;
        vms = std::max(vms, vs(cs[vIndexP+1]-cs[auxIndexQcs],vIndexP+1-auxIndexQcs,n));
      }
      for (int i = vSizeLCQ-1; i > 0; --i) { // along lower side
        auxIndexQcs = n-indexLoConvQ[i-1]-1;
        vms = std::max(vms, vs(cs[vIndexP+1]-cs[auxIndexQcs],vIndexP+1-auxIndexQcs,n));
      }
    } else {
      auxSizeUKQ = 0;
      for (int i = vSizeUCQ-1; i > 0; --i) { // along upper side
        if (cw(0, 0, indexUpNbP[vIndexP]-vIndexP,
               cs[indexUpNbP[vIndexP]+1]-cs[vIndexP+1],
                                           indexUpConvQ[i-1]-indexUpConvQ[i],
                                                                         cs[n-indexUpConvQ[i]-1]-cs[n-indexUpConvQ[i-1]-1])) {
          break;
        } else {
          auxIndexQcs = n-indexUpConvQ[i-1]-1;
          vms = std::max(vms, vs(cs[vIndexP+1]-cs[auxIndexQcs],vIndexP+1-auxIndexQcs,n));
          ++auxSizeUKQ;
        }
      }
      auxSizeLKQ = 0;
      for (int i = vSizeLCQ-1; i > 0; --i) { // along lower side
        if (acw(0, 0, indexLoNbP[vIndexP]-vIndexP,
                cs[indexLoNbP[vIndexP]+1]-cs[vIndexP+1],
                                            indexLoConvQ[i-1]-indexLoConvQ[i],
                                                                          cs[n-indexLoConvQ[i]-1]-cs[n-indexLoConvQ[i-1]-1])) {
          break;
        } else {
          auxIndexQcs = n-indexLoConvQ[i-1]-1;
          vms = std::max(vms, vs(cs[vIndexP+1]-cs[auxIndexQcs],vIndexP+1-auxIndexQcs,n));
          ++auxSizeLKQ;
        }
      }
      //      Remove interior points in K
      vSizeUCQ -= auxSizeUKQ;
      vSizeLCQ -= auxSizeLKQ;
      if (auxSizeUKQ > 0 && auxSizeLKQ == 0)
        indexUpConvQ[vSizeUCQ++] = indexLoConvQ[vSizeLCQ-1];
      else if (auxSizeUKQ == 0 && auxSizeLKQ > 0)
        indexLoConvQ[vSizeLCQ++] = indexUpConvQ[vSizeUCQ-1];
    }
  }
  
  return vms;
}
