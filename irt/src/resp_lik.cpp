#include <Rcpp.h>
#include "itempool_class_methods.h"
#include "prob.h"
using namespace Rcpp;

//#############################################################################@
//#############################################################################@
//########################### resp_lik #########################################
//#############################################################################@
//#############################################################################@


//#############################################################################@
//########################### resp_lik_bare_item_cpp ###########################
//#############################################################################@

// [[Rcpp::export]]
double resp_lik_bare_item_cpp(double resp, double theta, Rcpp::S4 item) {
  // This function calculates the response likelihood of an item for a
  // given response.

  // Deal with missing responses, return NA directly
  if (NumericVector::is_na(resp))
    return NA_REAL;

  // Get the Psychometric Model name
  std::string model = as<std::string>(item.slot("model"));

  if (model == "GPCM" || model == "GPCM2" || model == "PCM" || model == "GRM") {
    Rcpp::NumericVector P(resp);
    if (model == "GPCM" || model == "PCM" || model == "GPCM2") {
      P = prob_gpcm_bare_cpp(theta, item);
    } else if (model == "GRM") {
      P = prob_grm_bare_cpp(theta, item);
    }
    // The following line assumes that the resp goes from 0 to maximum number
    // of categories.
    return P[resp];
  } else if (model == "Rasch" || model == "1PL" || model == "2PL" ||
    model == "3PL" || model == "4PL") {
    double P = prob_4pm_bare_cpp(theta, item);
    // The following line is important (instead of second line) because 
    // it accounts for resp values that are not 0 or 1. 
    return pow(P, resp) * pow(1.0-P, 1.0-resp);
    return resp * P + (1 - resp) * (1 - P);
  }
  return NA_REAL;
}


//##############################################################################
//########################### resp_lik_item_cpp ################################
//##############################################################################
// [[Rcpp::export]]
Rcpp::NumericVector resp_lik_item_cpp(Rcpp::NumericVector resp,
                                      Rcpp::NumericVector theta, Rcpp::S4 item)
{
  // Calculate response likelihood for one item and multiple theta's (and
  // responses)
  unsigned int num_of_theta = theta.size();
  Rcpp::NumericVector output(num_of_theta);
  for(unsigned int i = 0; i < num_of_theta; i++)
    output[i] = resp_lik_bare_item_cpp(resp[i], theta[i], item);
  return output;
}


//##############################################################################
//########################### resp_lik_bare_testlet_cpp ########################
//##############################################################################
// [[Rcpp::export]]
double resp_lik_bare_testlet_cpp(Rcpp::NumericVector resp, double theta,
                                 Rcpp::S4 testlet)
{
  // Calculate response log-likelihood for a testlet and single examinee
  Rcpp::List item_list = as<List>(testlet.slot("item_list"));
  unsigned int num_of_items = item_list.size();
  double output = 0;
  Rcpp::S4 item; // This will be item
  for(unsigned int i = 0; i < num_of_items; i++) {
    item = as<Rcpp::S4>(item_list(i));
    if (!NumericVector::is_na(resp[i]))
        output = output + resp_lik_bare_item_cpp(
          resp(i), theta, as<Rcpp::S4>(item_list(i)));
  }
  return output;
}


//##############################################################################
//########################### resp_lik_testlet_cpp #############################
//##############################################################################
// [[Rcpp::export]]
Rcpp::NumericVector resp_lik_testlet_cpp(Rcpp::NumericMatrix resp,
                                         Rcpp::NumericVector theta,
                                         Rcpp::S4 testlet)
{
  // Calculate response log-likelihood for an Itempool and multiple
  // theta's (and response strings)
  unsigned int num_of_theta = theta.size();
  NumericVector output(num_of_theta);
  for(unsigned int i = 0; i < num_of_theta; i++) {
    // Get the row belong to the examinee. It is assumed that each row represents
    // an examinee.
    NumericVector resp_vector = resp(i, _);
    output[i] = resp_lik_bare_testlet_cpp(resp_vector, theta[i], testlet);
  }
  return output;
}


//#############################################################################@
//########################### resp_lik_bare_itempool_cpp ######################
//#############################################################################@

// [[Rcpp::export]]
double resp_lik_bare_itempool_cpp(Rcpp::NumericVector resp, double theta,
                                   Rcpp::S4 ip) {
  // This function calculates the response likelihood of an item pool for a
  // given response string and one theta.

  // Assuming that theta.size() ==  resp.size(), though it may not be the case
  // always. There might be an instance where more responses can happen if
  // item pool has testlets.
  int no_items = resp.size();
  double result = 1;
  S4 item;
  // Indicator variable for whether all responses are missing (true) or at
  // least there is one non-missing response (false).
  bool resp_all_na = true;
  List item_list = ip.slot("item_list");
  for (int i = 0; i < no_items; i++) {
    // iterate over non-missing responses
    if (!R_IsNA(resp[i])) {
      resp_all_na = false; // one non-na observed
      item = as<Rcpp::S4>(item_list[i]);
      result = result * resp_lik_bare_item_cpp(resp[i], theta, item);
    }
  }
  if (resp_all_na) result = NA_REAL;   // should it return 0 or NA?
  return result;
}


//##############################################################################
//########################### resp_lik_itempool_cpp ###########################
//##############################################################################
// [[Rcpp::export]]
Rcpp::NumericVector resp_lik_itempool_cpp(Rcpp::NumericMatrix resp,
                                           Rcpp::NumericVector theta,
                                           Rcpp::S4 ip)
{
  // Calculate response log-likelihood for an Itempool and multiple
  // theta's (and response strings)
  unsigned int num_of_theta = theta.size();
  NumericVector output(num_of_theta);
  for(unsigned int i = 0; i < num_of_theta; i++) {
    // Get the row belong to the examinee. It is assumed that each row represents
    // an examinee.
    NumericVector resp_vector = resp(i, _);
    output[i] = resp_lik_bare_itempool_cpp(resp_vector, theta[i], ip);
  }
  return output;
}


