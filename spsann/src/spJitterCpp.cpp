/*******************************************************************************
Candidates for random perturbation of spatial coordinates
x: coordinates of the sample points
y: coordinates of the candidate locations
xmax, xmin, ymax, ymin: manimum and minimum shift in the x and y coordinates
idx: point to be jittered
*******************************************************************************/

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export(".spJitterCpp")]]

IntegerVector spJitterCpp(NumericMatrix x, NumericMatrix y, double xmax,
                          double xmin, double ymax, double ymin, int idx) {
   int ncol = 2, nrow = y.nrow(), i;
   NumericVector pt0(ncol, 0.0000), dx(ncol, 0.0000), dy(ncol, 0.0000);
   IntegerVector pt1(nrow, 0);
   
   /* Get the coordinates of the point to be jittered */
   pt0[0] = x(idx - 1, 0);
   pt0[1] = x(idx - 1, 1);
   
   /* Calculate the maximum shift in the x and y coodinates - search window */
   dx[0] = pt0[0] - xmax + xmin;
   dx[1] = pt0[0] + xmax + xmin;
   dy[0] = pt0[1] - ymax + ymin;
   dy[1] = pt0[1] + ymax + ymin;
   
   /* Get the row index of the candidate points that fall within the search
      window */
   for (i = 0; i < nrow; i++) {
     if (y(i, 0) >= dx[0] && y(i, 0) <= dx[1] && 
         y(i, 1) >= dy[0] && y(i, 1) <= dy[1]) {
       pt1[i] = i + 1;
     }
   }
   
//   /* Sort vector of row indexes and count the number of zeros */
//   std::sort(pt1.begin(), pt1.end(), std::greater<int>());
//   int nz = std::count(pt1.begin(), pt1.end(), 0);
//   
//   /* Copy nonzero values */
//   IntegerVector pt2(nrow - nz);
//   for (i = 0; i < nrow - nz; i++) {
//     pt2[i] = pt1[i];
//   }
   
//   return (pt2);
   return (pt1);
}
