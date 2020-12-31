#include <Rcpp.h>
#include "itempool_class_methods.h"
#include "sim_resp.h"
#include "ability_estimation.h"
#include "cat_select_next_item.h"
#include "info.h"
using namespace Rcpp;


//#############################################################################@
//########################### next_step_cat_cpp ################################
//#############################################################################@
// Initiate CAT or move CAT to the next step
//
// @description NA
//
Rcpp::List next_step_cat_cpp(Rcpp::List true_ability,
                             Rcpp::List cd,
                             Rcpp::Nullable<Rcpp::List> est_history,
                             Rcpp::Nullable<Rcpp::List> additional_args) {

  // Create est_history  (estimate history) object. est_history_step is just
  // one step of est_history list.
  // est_history is a list that holds each step of the adaptive test. The
  // first element is "1" which represents the beginning of the test.
  //   "est_before" is the estimated ability before the item's administration
  //   "se_before" is the estimated standard error before the item's administration
  //   "item" is the item that will be administered
  //   "testlet" is the testlet object that administered item belongs to
  //   "resp" is the response to the item that will be administered
  //   "est_after" is the estimated ability after the item's administration
  //   "se_after" is the estimated standard error after the item's administration
  Rcpp::List est_history_step = Rcpp::List::create(
    Rcpp::Named("est_before") = NA_REAL,
    // TODO: Ideally initiate it as the square root of the prior variance
    Rcpp::Named("se_before") = NA_REAL,
    Rcpp::Named("testlet") = R_NilValue,
    Rcpp::Named("item") = R_NilValue,
    Rcpp::Named("resp") = NA_INTEGER,
    Rcpp::Named("est_after") = NA_REAL,
    Rcpp::Named("se_after") = NA_REAL
    );
  Rcpp::List eh;
  // 'additional_args' is a list that will be passed to intermediate functions
  // such as "select_next_item_cpp", "est_ability_cat_cpp", "terminate_cat_cpp"
  // and hold arguments that should be persistent across all steps of adaptive
  // test. Currently only implemented in "select_next_item_cpp" but will be
  // implemented on other functions as needed.
  Rcpp::List aa; // Surragote for 'additional_args', initiate as an empty list

  // Decide whether this is initial item or next item
  if (Rf_isNull(est_history)) { // Initial Item

    std::string ability_type = cd["ability_type"];
    //////////  Initial Ability Estimate //////////
    // Rcout << "(cat_sim_single_cpp) - Initial Ability Estimate" << std::endl;
    if (cd.containsElementNamed("first_item_rule") &&
        !Rf_isNull(cd("first_item_rule"))) {
        // Get first_item_rule and first_item_par
        std::string first_item_rule = cd["first_item_rule"];
        Rcpp::List first_item_par = cd.containsElementNamed("first_item_par") ?
          cd["first_item_par"] : List();
        if (first_item_rule == "fixed_theta") {
          est_history_step["est_before"] = first_item_par["theta"];
        } else {
          est_history_step["est_before"] = 0;
        }
    }
    eh = Rcpp::List::create(est_history_step);
  } else { // Next item
    aa = Rcpp::as<List>(additional_args);
    eh = Rcpp::as<List>(est_history);

    int item_no = eh.size();  // The stage of the test.
    Rcpp::List est_history_last_step = eh[item_no-1];
    est_history_step["est_before"] = est_history_last_step["est_after"];
    est_history_step["se_before"] = est_history_last_step["se_after"];
    eh.push_back(est_history_step);

  }
  return List::create(Named("est_history") = eh,
                      Named("additional_args") = aa);
}

//#############################################################################@
//########################### generate_cat_resp_cpp ############################
//#############################################################################@
// Generate response for CAT
//
// @description This function generates response for CAT given est_history and
//   additional_pars
//
//
//
// [[Rcpp::export]]
Rcpp::List generate_cat_resp_cpp(Rcpp::List true_ability,
                                 Rcpp::List cd,
                                 Rcpp::List est_history,
                                 Rcpp::List additional_args) {
  // Rcout << "generate_cat_resp_cpp - Stage 1" << std::endl;
  double ta = as<double>(true_ability[0]);
  Rcpp::List eh = clone(est_history);
  int item_no = eh.size();  // The stage of the test.
  // Rcout << "generate_cat_resp_cpp - Stage 2" << std::endl;
  Rcpp::List est_history_last_step = eh[item_no-1];
  // Rcout << "generate_cat_resp_cpp - Stage 3" << std::endl;
  est_history_last_step["resp"] = sim_resp_bare_cpp(
    ta, est_history_last_step("item"));
  // Rcout << "generate_cat_resp_cpp - Stage 4" << std::endl;
  eh[item_no-1] = est_history_last_step; // Update est_history
  // Rcout << "generate_cat_resp_cpp - Stage 5" << std::endl;
  return List::create(Named("est_history") = eh,
                      Named("additional_args") = additional_args);
}


//#############################################################################@
//########################### est_ability_cat_cpp ##############################
//#############################################################################@
// This function estimates the ability so far in CAT given the cat design (cd)
// and the estimate history (est_history).
// last_estimate: an indicator for whether the this is the last estimate before
//   terminating the CAT.
// The output will be a list with two elements: "est" and "se"
// [[Rcpp::export]]
Rcpp::List est_ability_cat_cpp(Rcpp::List true_ability,
                               Rcpp::List cd,
                               Rcpp::List est_history,
                               Rcpp::List additional_args,
                               bool last_estimate = false) {
  // Rcout << "Stage 1" << "\n";
  Rcpp::List eh = clone(est_history);
  Rcpp::List aa = clone(additional_args);
  int item_no = eh.size();  // The stage of the test.
  Rcpp::List est_history_last_step = eh[item_no-1];

  Rcpp::List cd_steps = cd["step"];
  Rcpp::List cdi = cd_steps[item_no-1];
  // Get ability estimation rule
  std::string ability_est_rule = cdi["ability_est_rule"];
  // Get ability estimation parameter
  Rcpp::List ability_est_par;
  if (cdi.containsElementNamed("ability_est_par"))
    ability_est_par = cdi["ability_est_par"];
  // If this is the last ability estimation, than use the method that is
  // designatef for final ability estimate.
  if (last_estimate) {
    ability_est_rule = as<std::string>(cd["final_ability_est_rule"]);
    if (cd.containsElementNamed("final_ability_est_par"))
      ability_est_par = as<Rcpp::List>(cd["final_ability_est_par"]);
  }
  // Rcout << "Stage 2" << "\n";
  // Get Administered item pool and response string
  NumericVector resp(item_no);
  Rcpp::List administered_ip_list(item_no);
  Rcpp::List temp_list;
  for (int i=0; i < item_no; i++) {
    temp_list = eh[i];
    resp[i] = temp_list["resp"];
    administered_ip_list[i] = temp_list["item"];
  }
  // Rcout << "Stage 3" << "\n";
  S4 administered_ip("Itempool");
  administered_ip.slot("item_list") = administered_ip_list;
  // Rcout << "Stage 4" << "\n";
  if (ability_est_rule == "eap") {
    // Rcout << "  Stage 4.1 - EAP" << "\n";
    // Get parameters for the eap
    std::string prior_dist = ability_est_par("prior_dist");
    NumericVector prior_par = ability_est_par("prior_par");
    double min_theta = ability_est_par("min_theta");
    double max_theta = ability_est_par("max_theta");
    NumericVector theta_range = NumericVector::create(min_theta, max_theta);
    int no_of_quadrature = ability_est_par("no_of_quadrature");
    // Rcout << "    Stage 4.1.2" << "\n";
    temp_list = est_ability_eap_single_examinee_cpp(
      resp, administered_ip, theta_range, no_of_quadrature, prior_dist,
      prior_par);
    est_history_last_step["est_after"] = temp_list["est"];
    est_history_last_step["se_after"] = temp_list["se"];
    eh[item_no-1] = est_history_last_step; // Update est_history
    return List::create(Named("est_history") = eh,
                        Named("additional_args") = aa);
  } else if (ability_est_rule == "owen") {
    // Rcout << "  Stage 4.2 - Owen" << "\n";
    // Get parameters for the owen
    double prior_mean = ability_est_par("prior_mean");
    double prior_var = ability_est_par("prior_var");
    temp_list = est_ability_owen_cpp(administered_ip, resp, prior_mean,
                                     prior_var);

    est_history_last_step["est_after"] = temp_list["est"];
    est_history_last_step["se_after"] = temp_list["se"];
    eh[item_no-1] = est_history_last_step; // Update est_history
    return List::create(Named("est_history") = eh,
                        Named("additional_args") = aa);
  } else if (ability_est_rule == "ml") {
    // Rcout << "  Stage 4.2 - ml" << "\n";
    double min_theta = ability_est_par("min_theta");
    double max_theta = ability_est_par("max_theta");
    double criterion = ability_est_par("criterion");
    double est = est_ability_4pm_nr_itempool_cpp(
      resp, administered_ip, NumericVector::create(min_theta, max_theta), 
      criterion);
    Rcpp::NumericVector info = info_itempool_bare_cpp(est, administered_ip,
                                                       true, false, resp);
    // double se = 1/sqrt(info[0]);
    // return Rcpp::List::create(Named("est") = est , Named("se") = se);
    //
    est_history_last_step["est_after"] = est;
    est_history_last_step["se_after"] = 1/sqrt(info[0]);
    eh[item_no-1] = est_history_last_step; // Update est_history
    return List::create(Named("est_history") = eh,
                        Named("additional_args") = aa);
  } else if (ability_est_rule == "sum_score") {
    // Rcout << "  Stage 4.2 - sum_score" << "\n";
    double est = 0;
    for (int i=0; i < item_no; i++) {
      est = est + as<double>(as<Rcpp::List>(eh[i])["resp"]);
    }
    est_history_last_step["est_after"] = est;
    est_history_last_step["se_after"] = NA_REAL;
    eh[item_no-1] = est_history_last_step; // Update est_history
    return List::create(Named("est_history") = eh,
                        Named("additional_args") = aa);
  }
  // Rcout << "Stage 5" << "\n";
  return temp_list;
}

//#############################################################################@
//########################### terminate_cat_cpp ################################
//#############################################################################@
//' Function determines whether to terminate CAT.
//'
//' @description This function returns either \code{true} or \code{false} where
//'   \code{true} indicates to terminate the test and \code{false} indicates to
//'   terminate the test.
//'
//'   If there is only one condition, test will end when the condition
//'   satisfied. If there are multiple conditions, all of them should be
//'   satisfied in order for test to terminate.
//'
//'
//' @param true_ability True ability of the examinee.
//' @param cd A \code{cat_design} object that holds the test specifications
//'   of the CAT.
//' @param est_history is a \code{List} that holds each step of the adaptive
//'   test. The first element is "1" which represents the beginning of the test.
//'   The elements are:
//'   \describe{
//'     \item{\code{"est_before"}}{The estimated ability before the item's
//'       administration.}
//'     \item{\code{"se_before"}}{The estimated standard error before the
//'       item's administration.}
//'     \item{\code{"item"}}{The item object that will be administered.}
//'     \item{\code{"testlet"}}{The testlet object that the administered item
//'       belongs to.}
//'     \item{\code{"resp"}}{The response value of the item that is
//'       administered}
//'     \item{\code{"est_after"}}{The estimated ability after the item's
//'       administration.}
//'     \item{\code{"se_after"}}{The estimated standard error after the
//'       item's administration}
//'   }
//'
//' @param additional_args Additional arguments
//' 
//' @noRd
//' 
// [[Rcpp::export]]
bool terminate_cat_cpp(Rcpp::List true_ability,
                       Rcpp::List cd,
                       Rcpp::List est_history,
                       Rcpp::List additional_args) {
  // This function determines whether to end CAT test or not.
  // This function will run after each step of CAT test, i.e. runs after item
  // is selected and administered.
  Rcpp::List eh = clone(est_history);
  Rcpp::List cd_copy = clone(cd);
  Rcpp::List temp_list;
  int item_no = eh.size();  // The stage of the test.
  // Rcout << "  terminate_test  --   1  --  item_no = " << item_no << std::endl;
  int max_test_length = cd_copy["max_test_length"];
  // If test reached it's maximum value, terminate the test, this will obviate
  // the need for checking whether there is an available item unadministered
  // item pool because max_test_length bounded with the size of the item pool,
  // and if item_no = max_test_length, it means there are no items left in
  // the item pool.
  if (item_no == max_test_length) return true;
  bool terminate_test = true;

  // Get termination rule
  StringVector termination_rules = Rcpp::as<StringVector>(cd_copy("termination_rule"));
  // Get termination parameters
  Rcpp::List termination_pars = cd_copy["termination_par"];

  // Check whether there is at least one unadministered item


  int no_of_rules = termination_rules.size();  // Number of termination rules
  for (int i = 0; i < no_of_rules; i++) {
    if (termination_rules[i] == "min_item") {
      // Rcout << "  terminate_test  --   min_item " << std::endl;

      temp_list = termination_pars["min_item"];
      int min_item = as<int>(temp_list["min_item"]);
      if (item_no < min_item) {
        return false;
      } else
        terminate_test = true;
    } else if (termination_rules[i] == "max_item") {
      // Rcout << "  terminate_test  --   max_item " << std::endl;

      temp_list = termination_pars["max_item"];
      int max_item = as<int>(temp_list["max_item"]);
      if (item_no >= max_item) {
        return true;
      } else
        terminate_test = false;
    } else if (termination_rules[i] == "min_se") {
      // Rcout << "  terminate_test  --   min_se " << std::endl;

      temp_list = termination_pars["min_se"];
      double min_se = as<double>(temp_list["min_se"]);
      Rcpp::List est_history_step = eh[item_no-1];
      double current_se = as<double>(est_history_step["se_after"]);
      // Rcout << "  terminate_test  --   min_se = " << min_se << "  --  current_se = " << current_se << std::endl;
      if (current_se <= min_se) {
        return true;
      } else {
        terminate_test = false;
      }
    } else if (termination_rules[i] == "sprt") {
      // Extract the parameters needed by the method.
      // double cut_score = as<int>(termination_pars["cut_score"]);
      // theta_0 is the highest value for a failing examinee, i.e. examinees
      // whose theta values smaller than this value should fail.
      // theta_1 is the lowest value for a passing examinee, i.e. examinees
      // whose theta values larger than this value should pass.

      // Rcout << "  terminate_test  --   SPRT " << std::endl;

      temp_list = termination_pars["sprt"];
      double theta_0 = as<double>(temp_list["theta_0"]);
      double theta_1 = as<double>(temp_list["theta_1"]);
      double alpha = as<double>(temp_list["alpha"]);
      double beta = as<double>(temp_list["beta"]);
      // Rcout << "  terminate_test  --   theta_0 = " << theta_0 << "  --  theta_1 " << theta_1 << "  --  alpha = " << alpha << std::endl;
      // Calculate A and B:
      double A = (1 - beta) / alpha;
      double B = beta/(1-alpha);
      // Calculate the log-likelihood at theta_0 and theta_1
      // Rcout << "  terminate_test  --   log(theta_1) = " << loglik_est_history(eh, theta_1, true) << "  --  log(theta_0) = " << loglik_est_history(eh, theta_0, true) << std::endl;
      double loglik_ratio = loglik_est_history(eh, theta_1, true) -
        loglik_est_history(eh, theta_0, true);
      // Rcout << "  terminate_test  --   log(B) = " << log(B) << "  --  loglik_ratio = " << loglik_ratio << "  --  log(A) = " << log(A) << std::endl;
      if (loglik_ratio > log(A) || loglik_ratio < log(B)) {
        terminate_test = true;
      } else terminate_test = false;
    }
  }
  return terminate_test;
}



//#############################################################################@
//########################### calculate_exposure_rates_cpp #####################
//#############################################################################@
// This function calculates the exposure rate of a given list of cat_output.
// [[Rcpp::export]]
Rcpp::NumericVector calculate_exposure_rates_cpp(Rcpp::StringVector item_ids,
                                                 Rcpp::List cat_output_list) {
  // Get the size of ip
  unsigned int ip_size = item_ids.size(); // item pool size, testlets and items
  int n_sim = cat_output_list.size(); // number of simulations

  // The numeric vector that holds exposure rates:
  NumericVector er(ip_size);
  er.attr("names") = item_ids;
  List est_history, est_history_step;
  std::string item_id, testlet_id, old_testlet_id;
  S4 item("Item");
  S4 testlet("Testlet");
  // Rcout << "calculate_overlap_rates_cpp == 3 " << std::endl;
  // iterate through cat_output_list:
  for (int i = 0; i < n_sim; i++ ) {
    // Rcout << "  calculate_overlap_rates_cpp == 3.1 - i = " << i << std::endl;
    est_history = cat_output_list[i];
    est_history = est_history["est_history"];
    // old_testlet_id holds the id of previously administered testlet. If
    // it is "" it means that no testlet administered before.
    old_testlet_id = "";
    // Iterate through est_history:
    for (int j = 0; j < est_history.size(); j++ ) {
      // Rcout << "  calculate_overlap_rates_cpp == 3.2 - j = " << j << std::endl;
      est_history_step = est_history[j];
      // Check if testlet exists, if yes, all items in that testlet are
      // assumed to be administered, so only testlet exposure rates
      // will be calculated.
      if (Rf_isNull(est_history_step("testlet"))) { // there is no testlet
        // Rcout << "    calculate_overlap_rates_cpp == 3.2.1 " << std::endl;
        item = as<S4>(est_history_step("item"));
        item_id = as<std::string>(item.slot("id"));
        // Rcout << "    calculate_overlap_rates_cpp == 3.2.1.1 " << std::endl;
        //idx = std::find(item_ids.begin(), ip)
        er[item_id] = as<double>(er[item_id]) + 1;
        // set old_testlet_id to "" to indicate that previous item was not
        // a testlet.
        old_testlet_id = "";
      } else { // There is a testlet
        // Rcout << "    calculate_overlap_rates_cpp == 3.2.2 " << std::endl;
        testlet = as<S4>(est_history_step("testlet"));
        testlet_id = as<std::string>(testlet.slot("id"));
        if (old_testlet_id != testlet_id) {
          er[testlet_id] = as<double>(er[testlet_id]) + 1;
        }
        old_testlet_id = testlet_id;
      }
    }
  }
  for (unsigned int i = 0; i < ip_size; i++ )
    er[i] = er[i] / n_sim;
  return er;
}



//########################### nChoosek #########################################
unsigned nChoosek( unsigned n, unsigned k )
{
  // I copied this function from: https://stackoverflow.com/a/9331125/2275286
  if (k > n) return 0;
  if (k * 2 > n) /*return*/ k = n-k;  //remove the commented section
  if (k == 0) return 1;

  int result = n;
  for(unsigned int i = 2; i <= k; ++i ) {
      result *= (n-i+1);
      result /= i;
  }
  return result;
}


//#############################################################################@
//########################### calculate_overlap_rates_cpp ######################
//#############################################################################@
// [[Rcpp::export]]
Rcpp::NumericVector calculate_overlap_rates_cpp(Rcpp::StringVector item_ids,
                                                Rcpp::List cat_output_list) {
  // This function calculates the average overlap rate (i.e. proportion of
  // items shared by pairs of examinees) of a given list of cat_output.
  // Rcout << "calculate_overlap_rates_cpp == 1 " << std::endl;

  // Get the size of ip
  unsigned int ip_size = item_ids.size();
  unsigned n_sim = cat_output_list.size(); // number of simulations

  Rcpp::StringVector item_ids_j;
  // The numeric vector that holds exposure rates:
  NumericVector overlap_rate(ip_size);
  overlap_rate.attr("names") = item_ids;
  List est_history, est_history_step;
  std::string item_id, testlet_id, old_testlet_id, temp_id;
  S4 item("Item");
  S4 testlet("Testlet");
  // Rcout << "calculate_overlap_rates_cpp == 2 " << std::endl;
  // This list will hold all of the administered id's of each test
  List item_id_list(n_sim);
  // Get the administered id's of all tests and store it in item_id_list
  for (unsigned int i = 0; i < n_sim; i++ ) {
    est_history = cat_output_list[i];
    est_history = est_history["est_history"];
    // old_testlet_id holds the id of previously administered testlet. If
    // it is "" it means that no testlet administered before.
    old_testlet_id = "";
    Rcpp::StringVector item_ids_i;
    // Iterate within an examinee's test
    for (int j = 0; j < est_history.size(); j++) {
      est_history_step = est_history[j];
      if (Rf_isNull(est_history_step("testlet"))) { // there is no testlet
        item = as<S4>(est_history_step("item"));
        item_ids_i.push_back(as<std::string>(item.slot("id")));
        old_testlet_id = "";
      } else { // There is a testlet
        testlet = as<S4>(est_history_step("testlet"));
        testlet_id = as<std::string>(testlet.slot("id"));
        if (old_testlet_id != testlet_id) {
          item_ids_i.push_back(testlet_id);
          // item_ids_i[j] = testlet_id;
        }
        old_testlet_id = testlet_id;
      }
    }
    item_id_list[i] = item_ids_i;
  }
  Rcpp::StringVector item_ids_i;

  // Rcout << "calculate_overlap_rates_cpp == 3 " << std::endl;
  // iterate through cat_output_list:
  for (unsigned int i = 0; i < (n_sim - 1); i++ ) {
    // Rcout << "  calculate_overlap_rates_cpp == 3.1 - i = " << i << std::endl;
    item_ids_i = as<Rcpp::StringVector>(item_id_list[i]);
    for (unsigned int j = i+1; j < n_sim; j++ ) {
      // Rcout << "    calculate_overlap_rates_cpp == 3.1.1 - j = " << j << std::endl;
      item_ids_j = as<Rcpp::StringVector>(item_id_list[j]);
      // iterate through the items if item_ids_i to see whether it overlap
      //0 with items in item_ids_j
      for (int k = 0; k < item_ids_i.size(); k++) {
        temp_id = as<std::string>(item_ids_i[k]);
        // Rcout << "    calculate_overlap_rates_cpp == 3.1.1.1 - k = " << k << "/" << item_ids_i.size() <<" ;  id = " << temp_id << std::endl;
        // Checks whether item_ids_i[k] in item_ids_j
        if (std::find(item_ids_j.begin(), item_ids_j.end(), temp_id) != item_ids_j.end())
          overlap_rate[temp_id] = as<double>(overlap_rate[temp_id]) + 1;
        // Rcout << "    calculate_overlap_rates_cpp == 3.1.1.2 " << std::endl;
      }
    }
  }
  unsigned no_of_pairs = nChoosek(n_sim, 2);
  for (unsigned int i = 0; i < ip_size; i++ )
    overlap_rate[i] = overlap_rate[i] / no_of_pairs;
  return overlap_rate;
}


