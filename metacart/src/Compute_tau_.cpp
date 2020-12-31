#include <Rcpp.h>
using namespace Rcpp;

//' Compute tau2 
//' 
//' @param x1 the effect size g in the unsplit leaves
//' @param x2 the sampling variance vi in the unsplit leaves
//' @param x3 the lable of the unsplit leaves
//' @param x4 the sorted effect size in the parent leaf
//' @param x5 the sorted sampling variance in the parent leaf
//' @param xuni the unique labels in the unsplit leaves
//' @keywords internal
// [[Rcpp::export(".compute_tau_")]]
NumericVector compute_tau_(NumericVector x1,NumericVector x2, 
                             NumericVector x3, NumericVector xuni,
                             NumericVector x4, NumericVector x5){
  // x1 is the effect size in the unsplit leaves
  // x2 is the sampling variance in the unsplit leaves
  // x3 is the lable of the unsplit leaves
  // xuni is the unique lables of the unsplit leaves
  // x4 is the sorted effect size in the parent leaf
  // x5 is the sorted sampling variance in the parent leaf
  
  //===============The estimate in the other leaves=================//
  int i;
  double sumQoleaf = 0;
  double sumColeaf = 0;
  double sWoleaf = 0;
  for (i = 0; i < xuni.length(); i++) {
    double tempWY2 = 0;
    double tempWY = 0;
    double tempW = 0;
    double tempW2 = 0;
    int ii;
    for (ii = 0; ii < x3.length(); ii++) {
      if (x3[ii] == xuni[i]) {
        tempWY2 = tempWY2 + pow(x1[ii], 2)/x2[ii];
        tempWY = tempWY + x1[ii]/x2[ii];
        tempW = tempW + 1/x2[ii];
        tempW2 = tempW2 + pow(x2[ii], -2);
      }
    }
    sumQoleaf = sumQoleaf + tempWY2 - pow(tempWY,2)/tempW;
    sumColeaf = sumColeaf + tempW2/tempW;
    sWoleaf = sWoleaf + tempW;
  }
  //===============The estimate in the child leaves=================//
  NumericVector wy;
  NumericVector wy2;
  NumericVector wts;
  NumericVector w2;
  int j;
  for (j = 0; j < x4.length(); j++) {
    wy.push_back(x4[j]/x5[j]);
    wy2.push_back(pow(x4[j],2)/x5[j]);
    wts.push_back(1/x5[j]);
    w2.push_back(pow(x5[j], -2));
  }
  NumericVector cwy = cumsum(wy);
  NumericVector cwy2 = cumsum(wy2);
  NumericVector cwts = cumsum(wts);
  NumericVector cw2 = cumsum(w2);
  NumericVector tau2;
  double swy = sum(wy);
  double swy2 = sum(wy2);
  double sw = sum(wts);
  double sw2 = sum(w2);
  for (j = 0; j < x4.length() - 1; j++){
    double tempQ = swy2-pow(cwy[j],2)/cwts[j]-pow(swy-cwy[j],2)/(sw-cwts[j]);
    double tempC = sw+sWoleaf-cw2[j]/cwts[j]-(sw2 - cw2[j])/(sw - cwts[j])-sumColeaf;
    double temptau2 = (tempQ+sumQoleaf-(x4.length()+x1.length() -xuni.length() -2))/tempC;
    if (temptau2 < 0) {
      temptau2 = 0;
    }
    tau2.push_back(temptau2);
  }
  
  return tau2;
}


