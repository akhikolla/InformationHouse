#include <Rcpp.h>
using namespace Rcpp;

// Get contiguous cells (rook case only)

// [[Rcpp::export(name = ".contigCells")]]
IntegerVector contigCells_cpp(int pt, int bgr, NumericMatrix mtx) {
  int rr;
  int cc;
  int dim1 = mtx.nrow();
  int dim2 = mtx.ncol();
  IntegerVector r(4);
  IntegerVector c(4);
  IntegerVector ad;
  int idx;
  if (pt % dim1 == 0) {
    rr = dim1;
    cc = pt / dim1;
  } else {
    cc = trunc(pt / dim1) + 1;
    rr = pt - (cc-1) * dim1;
  }
  r[0] = rr-1;
  r[1] = rr+1;
  r[2] = rr;
  r[3] = rr;
  c[0] = cc;
  c[1] = cc;
  c[2] = cc-1;
  c[3] = cc+1;
  for (int i = 0; i < 4; i++){
    if(r[i] > 0 && r[i] <= dim1 && c[i] > 0 && c[i] <= dim2){
      idx = r[i] + (c[i] - 1) * dim1;
      if(mtx[idx-1] == bgr){
        ad.push_back(idx);
      }
    }
  }
  return(ad);
}

// Assign values
// [[Rcpp::export(name = ".assignValues")]]
NumericMatrix assignValues_cpp(int val, IntegerVector ad, NumericMatrix mtx) {
  for (int i = 0; i < ad.length(); i++){
    mtx[ad[i]-1] = val;
  }
  return(mtx);
}

// Transpose index of input cells

// [[Rcpp::export(name = ".indexTranspose")]]
IntegerVector indexTranspose_cpp(IntegerVector id, int dim1, int dim2) {
  int n = id.size();
  int rr = 0;
  int cc = 0;
  IntegerVector out(n);
  for (int i = 0; i < n; i++){
    if(id[i] % dim1 == 0){
      rr = dim1;
      cc = id[i] / dim1;
    } else {
      cc = trunc(id[i] / dim1) + 1;
      rr = id[i] - (cc - 1) * dim1;
    }
    out[i] = cc + (rr-1) * dim2;
  }
  return(out);
}

// Remove single tones

// // [[Rcpp::export(name = "rmSingle")]]
// IntegerVector rmSingle_cpp(NumericMatrix mtx, bool rm) {
//   int dim1 = mtx.ncol();
//   int dim2 = mtx.nrow();
//   int singles;
//   int v;
//   int vval;
// //
//   int n = mtx.length();
//   for (int i = 0; i < n; i++){
//     v[i] = mtx[i];
//     vval[i] = mtx[i];
//   v <- which(is.finite(v))
//   for (pt in v){  ## Faster than sapply or vapply!
//     if(pt %% dim2 == 0){
//       cc <- dim2
//       rr <- pt / dim2
//     } else {
//       rr <- trunc(pt / dim2) + 1
//       cc <- pt - (rr-1) * dim2
//     }
//     ad <- c(rr-1,rr+1,rr,rr,cc,cc,cc-1,cc+1)
//       ad[ad	 <= 0 | c(ad[1:4] > dim1, ad[5:8] > dim2)] <- NA
//       ad <- ad[5:8] + (ad[1:4]-1)*dim2
//       ad <- ad[is.finite(ad)]
//     if(all(.subset(vval, ad) != .subset(vval, pt))){
//       if(rm == TRUE){
//         vval[pt] <- sample(.subset(vval, ad), 1)
//       } else {
//         singles <- c(singles, pt)
//       }
//     }
//   }
//   if(rm == TRUE){
//     rst[] <- vval
//     return(rst)
//   } else {
//     return(singles)
//   }
// }
