
#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;
using namespace std;

double sign(double a);

double scalar_quantile(double p, int n, int r);
NumericVector vector_quantile(double p, int n, int r);

double pens(int n, int lens, double q); 
double penfs(int n, int len, NumericVector q);
double eta(int llen, int ulen, int n, double sd, NumericVector q);






