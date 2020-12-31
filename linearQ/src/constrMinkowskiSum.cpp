//Coded by Chengcheng Huang and Housen Li

//----------------------------------------------------------------------------------------------//
// Compute (P + Q)_{x1 >= thd} in R^2
// For detailed algorithm see
//      Bernholt, T., Eisenbrand, F., & Hofmeister, T. (2009). Constrained Minkowski
//          Sums: A Geometric Framework for Solving Interval Problems in Computational
//          Biology Efficiently. Discrete & Computational Geometry, 42(1), 22–36.
//----------------------------------------------------------------------------------------------//

#include "auxFun.h"

// satisfy constraint (1) or otherwise (0) for (p1, p2)
bool fitConstraint(double p1x1, double p2x1, double thd) {
  return (p1x1 + p2x1) >= thd;
}

/*************
* Convex hull of constrained (Minkowski) sum of two sets in R^2
*      R = convHull( (P + Q)_{x1 >= thd} )
* Input: (|P| <= |Q|, solution is non-empty, P & Q ordered increasingly wrt x1)
*   Px1,x2:    numeric vectors specifying point set P
*   Qx1,x2:    numeric vectors specifying point set Q
*   thd:       theshold in the constraint
* Output:
*   R:         data frame (index of P, index of Q)
* Note:
*   Feasibility of input is checked in R codes
****************/
// [[Rcpp::export(".convHullConstrSum")]]
DataFrame convHullConstrSum(const NumericVector &Px1, const NumericVector &Px2,
                            const NumericVector &Qx1, const NumericVector &Qx2, double thd) {
  int sizeP = Px1.size();
  int sizeQ = Qx1.size();
  
  // Determine the valid part of P
  int stIndexP = 0;
  while (!fitConstraint(Px1[stIndexP], Qx1[sizeQ-1], thd) && stIndexP < sizeP)
    ++stIndexP;
  const int vSizeP = sizeP - stIndexP; // valid size of P
  
  int sizeR = 3*vSizeP + sizeQ - 2; // additional 'vSizeP' for possible repetition
  IntegerVector RindexP(sizeR);
  IntegerVector RindexQ(sizeR);
  
  // For each p in P, find the smallest index of Q
  // such that the constraint is satisfied
  IntegerVector JQ(vSizeP); // it is decreasing
  int auxIndexQ = sizeQ - 1;
  for (int i = 0; i < JQ.size(); ++i) {
    while (fitConstraint(Px1[stIndexP+i], Qx1[auxIndexQ], thd) && auxIndexQ >= 0)
      JQ[i] = auxIndexQ--;
    auxIndexQ = JQ[i];
  }
  
  // Find neighbours of each pi in convHull(pi,...,pn) for valid part of P
  // via Graham scan algorithm
  //      Upper neighbour point
  IntegerVector indexConvP(vSizeP);   // index for (upper/lower) convex hull of P
  IntegerVector indexUpNbP(vSizeP-1); // index for upper neighbour point of pi in convHull(pi,...,pn) for i < n
  int vSizeCP = 1;                    // valid size of 'indexConvP'
  indexConvP[0] = sizeP-1;            // the last point is always in the convex hull
  for (int i = sizeP-2; i >= stIndexP; --i) {
    while (vSizeCP > 1 &&
           acw(Px1[i], Px2[i],
               Px1[indexConvP[vSizeCP-1]], Px2[indexConvP[vSizeCP-1]],
                                              Px1[indexConvP[vSizeCP-2]], Px2[indexConvP[vSizeCP-2]])) {
      --vSizeCP;
    }
    indexConvP[vSizeCP++]  = i;
    indexUpNbP[i-stIndexP] = indexConvP[vSizeCP-2];
  }
  //      Lower neighbour point
  IntegerVector indexLoNbP(vSizeP-1);
  vSizeCP = 1; // note that indexConvP[0] == P.size()-1;
  for (int i = sizeP-2; i >= stIndexP; --i) {
    while (vSizeCP > 1 &&
           cw(Px1[i], Px2[i],
              Px1[indexConvP[vSizeCP-1]], Px2[indexConvP[vSizeCP-1]],
                                             Px1[indexConvP[vSizeCP-2]], Px2[indexConvP[vSizeCP-2]]))
      --vSizeCP;
    indexConvP[vSizeCP++]  = i;
    indexLoNbP[i-stIndexP] = indexConvP[vSizeCP-2];
  }
  
  // Iteratively update R
  int vSizeR   = 0;                    // valid size of R
  IntegerVector indexUpConvQ(sizeQ);   // index for upper convex hull of Q
  int vSizeUCQ = 0;                    // valid size of 'indexUpConvQ'
  IntegerVector indexLoConvQ(sizeQ);   // index for lower convex hull of Q
  int vSizeLCQ = 0;                    // valid size of 'indexLoConvQ'
  int auxSizeUKQ, auxSizeLKQ;          // # points in K in upper and lower parts severally
  auxIndexQ = sizeQ-1;
  for (int vIndexP = 0; vIndexP < vSizeP; ++vIndexP) {
    //      Find convex hull of Q
    for (int i = auxIndexQ; i >= JQ[vIndexP]; --i) {
      while (vSizeUCQ > 1 &&             // along upper side
             acw(Qx1[i], Qx2[i],
                 Qx1[indexUpConvQ[vSizeUCQ-1]], Qx2[indexUpConvQ[vSizeUCQ-1]],
                                                   Qx1[indexUpConvQ[vSizeUCQ-2]], Qx2[indexUpConvQ[vSizeUCQ-2]]))
        --vSizeUCQ;
      indexUpConvQ[vSizeUCQ++] = i;
      while (vSizeLCQ > 1 &&              // along lower side
             cw(Qx1[i], Qx2[i],
                Qx1[indexLoConvQ[vSizeLCQ-1]], Qx2[indexLoConvQ[vSizeLCQ-1]],
                                                  Qx1[indexLoConvQ[vSizeLCQ-2]], Qx2[indexLoConvQ[vSizeLCQ-2]]))
        --vSizeLCQ;
      indexLoConvQ[vSizeLCQ++] = i;
    }
    auxIndexQ = JQ[vIndexP] - 1;
    //      Add points to 'R'
    RindexP[vSizeR]   = vIndexP + stIndexP;    // include leftmost point of Q
    RindexQ[vSizeR++] = indexUpConvQ[vSizeUCQ-1];
    if (indexLoConvQ[vSizeLCQ-1] != indexUpConvQ[vSizeUCQ-1]) {
      RindexP[vSizeR]   = vIndexP + stIndexP; // include leftmost point of Q
      RindexQ[vSizeR++] = indexLoConvQ[vSizeLCQ-1];
    }
    if (vIndexP == vSizeP-1) { // last point in P: include all remaining points in Q
      for (int i = vSizeUCQ-1; i > 0; --i) { // along upper side
        RindexP[vSizeR]   = vIndexP + stIndexP;
        RindexQ[vSizeR++] = indexUpConvQ[i-1];
      }
      for (int i = vSizeLCQ-1; i > 0; --i) { // along lower side
        RindexP[vSizeR]   = vIndexP + stIndexP;
        RindexQ[vSizeR++] = indexLoConvQ[i-1];
      }
    } else {
      auxSizeUKQ = 0;
      for (int i = vSizeUCQ-1; i > 0; --i) { // along upper side
        if (cw(0, 0, Px1[indexUpNbP[vIndexP]]-Px1[vIndexP+stIndexP],
               Px2[indexUpNbP[vIndexP]]-Px2[vIndexP+stIndexP],
                                           Qx1[indexUpConvQ[i-1]]-Qx1[indexUpConvQ[i]],
                                                                     Qx2[indexUpConvQ[i-1]]-Qx2[indexUpConvQ[i]])) {
          break;
        } else {
          RindexP[vSizeR]   = vIndexP + stIndexP;
          RindexQ[vSizeR++] = indexUpConvQ[i-1];
          ++auxSizeUKQ;
        }
      }
      auxSizeLKQ = 0;
      for (int i = vSizeLCQ-1; i > 0; --i) { // along lower side
        if (acw(0, 0, Px1[indexLoNbP[vIndexP]]-Px1[vIndexP+stIndexP],
                Px2[indexLoNbP[vIndexP]]-Px2[vIndexP+stIndexP],
                                            Qx1[indexLoConvQ[i-1]]-Qx1[indexLoConvQ[i]],
                                                                      Qx2[indexLoConvQ[i-1]]-Qx2[indexLoConvQ[i]])) {
          break;
        } else {
          RindexP[vSizeR]   = vIndexP + stIndexP;
          RindexQ[vSizeR++] = indexLoConvQ[i-1];
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
  
  // Return value
  IntegerVector vRindexP(vSizeR);
  IntegerVector vRindexQ(vSizeR);
  for (int i = 0; i < vSizeR; ++i) {
    vRindexP[i] = RindexP[i];
    vRindexQ[i] = RindexQ[i];
  }
  return DataFrame::create(Named("indexP") = vRindexP, Named("indexQ") = vRindexQ);
}



