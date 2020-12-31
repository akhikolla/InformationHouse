#include <Rcpp.h>
using namespace Rcpp;


//#############################################################################@
//########################### avg_rank #########################################
//#############################################################################@

// I took this function from https://stackoverflow.com/a/39170596/2275286
class Comparator {
private:
  const Rcpp::NumericVector& ref;

  bool is_na(double x) const
  {
    return Rcpp::traits::is_na<REALSXP>(x);
  }

public:
  Comparator(const Rcpp::NumericVector& ref_)
    : ref(ref_)
  {}

  bool operator()(const int ilhs, const int irhs) const
  {
    double lhs = ref[ilhs], rhs = ref[irhs];
    if (is_na(lhs)) return false;
    if (is_na(rhs)) return true;
    return lhs < rhs;
  }
};

// [[Rcpp::export]]
Rcpp::NumericVector avg_rank(Rcpp::NumericVector x)
{
  R_xlen_t sz = x.size();
  Rcpp::IntegerVector w = Rcpp::seq(0, sz - 1);
  std::sort(w.begin(), w.end(), Comparator(x));

  Rcpp::NumericVector r = Rcpp::no_init_vector(sz);
  for (R_xlen_t n, i = 0; i < sz; i += n) {
    n = 1;
    while (i + n < sz && x[w[i]] == x[w[i + n]]) ++n;
    for (R_xlen_t k = 0; k < n; k++) {
      r[w[i + k]] = i + (n + 1) / 2.;
    }
  }
  return r;
}


//#############################################################################@
//########################### cat_sim_single_cpp ###############################
//#############################################################################@
// [[Rcpp::export]]
double biserial_cpp(Rcpp::NumericVector score,
                                 Rcpp::NumericVector total_score,
                                 std::string type = "default") {
  
  // Remove the NAs from the vectors
  // It is assumed that if score is NA total_score is automatically NA
  total_score = total_score[!is_na(score)];
  score = score[!is_na(score)];

  double n = score.size();
  Rcpp::NumericVector dev(n); // Deviation scores

  if (type == "rank") {
    // Kraemer pointed out that Cureton's (1958, 1964) rank biserial correlation
    // coefficient "essentially replaces observations with their ranks and then
    // applies Brogden's approach. " (p.280)
    total_score = avg_rank(total_score);
    type = "brogden";
  }
  if ((type == "clemans-lord") | (type == "brogden")) {
    // Calculate deviation
    dev = total_score - mean(total_score);
    NumericVector sorted_score = clone(score);
    std::sort(sorted_score.begin(), sorted_score.end(), std::greater<int>());
    double num = sum(score * dev);
    if ((type == "clemans-lord") & (num < 0))
      return -(num / sum((1 - sorted_score) * dev));
    // For positive values Lord's modification is the same as Brogden.
    return (sum(score * dev) /  sum(sorted_score * dev));
  }

  // Mean of total scores of  examinees who correctly answered the item
  NumericVector ts1 = total_score[score == 1];
  NumericVector ts0 = total_score[score == 0];
  double mu0 = mean(ts0);
  double mu1 = mean(ts1);
  double n0 = ts0.size();
  double n1 = ts1.size();
  double sigma = sd(total_score);
  double u = R::dnorm(R::qnorm(n1/n, 0, 1.0, 1, 0),  0, 1.0, 0);

  // Rcout << "rnorm " << R::qnorm(.75, 0, 1.0, 1, 0) << std::endl;
  // Rcout << "n: " << n << std::endl;
  // Rcout << "n0: " << n0 << std::endl;
  // Rcout << "mu0: " << mu0 << std::endl;
  // Rcout << "n1: " << n1 << std::endl;
  // Rcout << "mu1: " << mu1 << std::endl;
  // Rcout << "u[0] " << u << std::endl;
  // Rcout << "u.size " << u.size() << std::endl;
  // u = Rcpp::wrap(Rcpp::qnorm(u, 1.0, 0.0));
  // u = dnorm(u[0]);

  return ((mu1 - mu0) / sigma) * (n1 * n0 / (u * n * n));
}
