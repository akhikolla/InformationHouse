#include <Rcpp.h>
using namespace Rcpp;

//' Nash-Sutcliffe Efficiency
//'
//' @param yhat Model outputs
//' @param y Observations
//' @return NSE
//' @examples
//' NSE(rnorm(100), rnorm(100))
//' @export
//[[Rcpp::export]]
double NSE(NumericVector yhat, NumericVector y) {
  // Rcpp comes with a lot of syntactic sugars so you can write C++ code almost like
  // R code. Unfortunately it does not have the sugar '^' element-wise exponent yet, so
  // we use x * x instead of x^2.
  double ybar = mean(y);
  double rss = sum((y - yhat) * (y - yhat));
  double tss = sum((y - ybar) * (y - ybar));
  return 1 - rss/tss;
}

//' Normalized root-mean-square error
//'
//' RMSE is normalized by the normalization constant
//' @param yhat Model outputs
//' @param y Observations
//' @param normConst The normalization constant
//' @return normalized RMSE
//' @examples
//' x <- rnorm(100)
//' y <- rnorm(100)
//' nRMSE(x, y, sd(y))
//' @export
//[[Rcpp::export]]
double nRMSE(NumericVector yhat, NumericVector y, double normConst) {
  double rmse = sqrt(mean((y - yhat) * (y - yhat)));
  return rmse / normConst;
}

//' Pearson's correlation
//'
//' Calculate the Pearson's correlation using the numerically stable formulation (see References). Internal function.
//' @param x First variable
//' @param y Second variable
//' @return Pearson's correlation
//' @section Reference John D. Cook's article at https://www.johndcook.com/blog/2008/11/05/how-to-calculate-pearson-correlation-accurately/
//[[Rcpp::export]]
double corr(NumericVector x, NumericVector y) {
  double xbar = mean(x);
  double ybar = mean(y);
  double sx = sd(x);
  double sy = sd(y);
  NumericVector part1 = (x - xbar) / sx;
  NumericVector part2 = (y - ybar) / sy;
  return sum(part1 * part2) / ((x.size() - 1));
}

//' Kling-Gupta Efficiency
//'
//' @param yhat Model outputs
//' @param y Observations
//' @return KGE
//' @examples
//' KGE(rnorm(100), rnorm(100))
//' @export
//[[Rcpp::export]]
double KGE(NumericVector yhat, NumericVector y) {
  double mu = mean(y);
  double mu_hat = mean(yhat);
  double sigma = sd(y);
  double sigma_hat = sd(yhat);
  double r = corr(yhat, y);
  double alpha = sigma_hat / sigma;
  double beta = mu_hat / mu;
  double EDsq = (r - 1) * (r - 1) + (alpha - 1) * (alpha - 1) + (beta - 1) * (beta - 1);
  return 1 - sqrt(EDsq);
}

//' Reduction of Error
//'
//' @param yhat Model outputs in the validation set
//' @param y Observations in the validation set
//' @param yc_bar Mean observations in the calibration set
//' @return RE
//' @examples
//' x <- rnorm(100)
//' y <- rnorm(100)
//' yc_bar <- mean(x[1:50])
//' RE(x[51:100], y[51:100], yc_bar)
//' @export
//[[Rcpp::export]]
double RE(NumericVector yhat, NumericVector y, double yc_bar) {
  double rss = sum((y - yhat) * (y - yhat));
  double tss = sum((y - yc_bar) * (y - yc_bar));
  return 1 - rss/tss;
}
