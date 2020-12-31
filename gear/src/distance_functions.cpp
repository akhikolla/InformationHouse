#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

// great circle distance for a single set of
// longitude/latitude coordinates. A port
// of the C code by Roger Bivand in the sp::spDists
// function
// [[Rcpp::export]]
NumericMatrix gcdist1(NumericVector lon, NumericVector lat, double eps) {
  int n = lon.length();
  NumericMatrix distmat(n, n);

  double d2r = M_PI/180;
  NumericVector lonr = lon * d2r;
  NumericVector latr = lat * d2r;

  double a = 6378.137;
  double f = 1.0/298.257223563;
  double F, G, L, sinG2, cosG2, sinF2, cosF2, sinL2, cosL2, S, C, w, R, D, H1, H2;

  for (int i = 0; i < n; i++) {
    for(int j = 0; j < i; j++) {
      if (abs(lonr[i] - lonr[j]) < eps && abs(latr[i] - latr[j]) < eps) {
        distmat(i, j) = 0;
      } else {
        F = (latr[i] + latr[j])/2;
        G = (latr[i] - latr[j])/2;
        L = (lonr[i] - lonr[j])/2;

        sinG2 = pow(sin(G), 2);
        cosG2 = pow(cos(G), 2);
        sinF2 = pow(sin(F), 2);
        cosF2 = pow(cos(F), 2);
        sinL2 = pow(sin(L), 2);
        cosL2 = pow(cos(L), 2);

        S = sinG2 * cosL2 + cosF2 * sinL2;
        C = cosG2 * cosL2 + sinF2 * sinL2;

        w = atan(sqrt(S/C));
        R = sqrt(S*C)/w;

        D = 2 * w * a;
        H1 = (3*R - 1) / (2 * C);
        H2 = (3*R + 1) / (2 * S);

        distmat(i, j) = D * (1 + f * H1 * sinF2 * cosG2 - f * H2 * cosF2 * sinG2);
        distmat(j, i) = distmat(i, j);
      }
    }
  }
  return distmat;
}

// great circle distance for two sets of
// longitude/latitude coordinates. A port
// of the C code by Roger Bivand in the sp::spDists
// function
// [[Rcpp::export]]
NumericMatrix gcdist2(NumericVector lon1, NumericVector lat1,
                      NumericVector lon2, NumericVector lat2,
                      double eps) {
  int n = lon1.length();
  int m = lon2.length();
  NumericMatrix distmat(n, m);

  double d2r = M_PI/180;
  NumericVector lon1r = lon1 * d2r;
  NumericVector lon2r = lon2 * d2r;
  NumericVector lat1r = lat1 * d2r;
  NumericVector lat2r = lat2 * d2r;

  double a = 6378.137;
  double f = 1.0/298.257223563;
  double F, G, L, sinG2, cosG2, sinF2, cosF2, sinL2, cosL2, S, C, w, R, D, H1, H2;

  for (int i = 0; i < n; i++) {
    for(int j = 0; j < m; j++) {
      if (abs(lon1r[i] - lon2r[j]) < eps && abs(lat1r[i] - lat2r[j]) < eps) {
        distmat(i, j) = 0;
      } else {
        F = (lat1r[i] + lat2r[j])/2;
        G = (lat1r[i] - lat2r[j])/2;
        L = (lon1r[i] - lon2r[j])/2;

        sinG2 = pow(sin(G), 2);
        cosG2 = pow(cos(G), 2);
        sinF2 = pow(sin(F), 2);
        cosF2 = pow(cos(F), 2);
        sinL2 = pow(sin(L), 2);
        cosL2 = pow(cos(L), 2);

        S = sinG2 * cosL2 + cosF2 * sinL2;
        C = cosG2 * cosL2 + sinF2 * sinL2;

        w = atan(sqrt(S/C));
        R = sqrt(S*C)/w;

        D = 2 * w * a;
        H1 = (3*R - 1) / (2 * C);
        H2 = (3*R + 1) / (2 * S);

        distmat(i, j) = D * (1 + f * H1 * sinF2 * cosG2 - f * H2 * cosF2 * sinG2);
      }
    }
  }
  return distmat;
}

// compute Euclidean distance for a single set of
// coordinates
// [[Rcpp::export]]
NumericMatrix eucdist1(NumericVector x, NumericVector y, double eps) {
  int n = x.length();
  NumericMatrix distmat(n, n);

  for (int i = 0; i < n; i++) {
    for(int j = 0; j < i; j++) {
      if ( abs(x[i] - x[j]) < eps && abs(y[i] - y[j]) < eps) {
        distmat(i, j) = 0;
      } else {
        distmat(i, j) = sqrt(pow(x[i] - x[j], 2) + pow(y[i] - y[j], 2));
        distmat(j, i) = distmat(i, j);
      }
    }
  }
  return distmat;
}

// compute Euclidean distance for two sets of
// coordinates
// [[Rcpp::export]]
NumericMatrix eucdist2(NumericVector x1, NumericVector y1,
                       NumericVector x2, NumericVector y2,
                       double eps) {
  int n = x1.length();
  int m = x2.length();
  NumericMatrix distmat(n, m);

  for (int i = 0; i < n; i++) {
    for(int j = 0; j < m; j++) {
      if (abs(x1[i] - x2[j]) < eps && abs(y1[i] - y2[j]) < eps) {
        distmat(i, j) = 0;
      } else {
        distmat(i, j) = sqrt(pow(x1[i] - x2[j], 2) + pow(y1[i] - y2[j], 2));
      }
    }
  }
  return distmat;
}
