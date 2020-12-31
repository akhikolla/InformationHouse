#include <Rcpp.h>
#include "itempool_class_methods.h"
#include "prob.h"
using namespace Rcpp;

//##############################################################################
//##############################################################################
//########################### sim_resp_cpp #####################################
//##############################################################################
//##############################################################################

// [[Rcpp::export]]
int sim_resp_4pm_bare_cpp(double theta, Rcpp::S4 item) {
  double p = prob_4pm_bare_cpp(theta, item);
  double u = as<double>(runif(1, 0, 1));
  if (p > u) {
    return 1;
  } else return 0;
}

//##############################################################################
//########################### sim_resp_poly_bare_cpp ###########################
//##############################################################################

// [[Rcpp::export]]
int sim_resp_poly_bare_cpp(double theta, Rcpp::S4 item) {
  // This method is based on De Ayala (1994) The Influence of
  // Multidimensionality on the, p. 158-159: "The generation of an examinee’s
  // polytomous response string was accomplished by calculating the
  // probability of responding to each item alternative according to the
  // MGRM; the scaling factor D was set to 1.0. Based on the probability for
  // each alternative, cumulative probabilities were obtained for each
  // alternative. A random error component was incorporated into each
  // response by selecting a random number from a uniform distribution [0, 1]
  // and comparing it to the cumulative probabilities. The ordinal position
  // of the first cumulative probability that was greater than the random
  // number was taken as the examinee’s response to the item."
  double u = as<double>(runif(1, 0, 1));
  // Cumulative Probability
  Rcpp::NumericVector cp = cumsum(prob_poly_bare_cpp(theta, item));
  int cp_size = cp.size();
  for (int i=cp_size-2; i > -1; i--)
    // Rcout << i << " - u = " << u << "  , cp = " << cp[i] << "\n";
    if (u > cp[i]) return(i+1);
  return 0;
}


//##############################################################################
//########################### sim_resp_bare_cpp ################################
//##############################################################################

// [[Rcpp::export]]
int sim_resp_bare_cpp(double theta, Rcpp::S4 item) {
  std::string model = as<std::string>(item.slot("model"));
  if (model == "GRM" || model == "GPCM") {
    return sim_resp_poly_bare_cpp(theta, item);
    // The default value is 4pm.
  } else return sim_resp_4pm_bare_cpp(theta, item);
}

