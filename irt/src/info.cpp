#include <Rcpp.h>
#include "itempool_class_methods.h"
#include "prob.h"
#include "resp_loglik.h"
using namespace Rcpp;


//##############################################################################
//##############################################################################
//########################### irt4PMInfo #######################################
//##############################################################################
//##############################################################################

//##############################################################################
//########################### info_4pm_bare_cpp ################################
//##############################################################################

// [[Rcpp::export]]
double info_4pm_bare_cpp(double theta, Rcpp::S4 item)
{
  // This function calculates the expected information of a single item
  // for one theta.
  Rcpp::List parList = as<List>(item.slot("parameters"));
  std::string model = as<std::string>(item.slot("model"));
  double a = 1, b = as<double>(parList["b"]), c = 0, d = 1, D = 1;
  if (model != "Rasch") {
    D = as<double>(parList["D"]);
    if ((model == "2PL") || (model == "3PL") || (model == "4PL")) {
      a = as<double>(parList["a"]);
      if ((model == "3PL") || (model == "4PL")) {
        c = as<double>(parList["c"]);
        if (model == "4PL")
          d = as<double>(parList["d"]);
      }
    }
  }
  return ((D * a)*(D * a) * (d - c)*(d - c)) /
       ((c + d * exp(D * a * (theta - b))) *
       (1 - c + (1-d) * exp(D*a * (theta - b))) *
       (1 + exp(-D * a * (theta - b))) * (1 + exp(-D * a * (theta - b))));
  // return (pow(D * a, 2) * pow(d - c, 2)) /
  //      ((c + d * exp(D * a * (theta - b))) *
  //      (1 - c + (1-d) * exp(D*a * (theta - b))) *
  //      pow(1 + exp(-D * a * (theta - b)), 2));
}

//##############################################################################
//########################### info_4pm_matrix_cpp ##############################
//##############################################################################

// [[Rcpp::export]]
NumericMatrix info_4pm_matrix_cpp(NumericVector theta, NumericMatrix ip,
                                  bool tif)
{
  int num_of_theta = theta.size();
  int num_of_items = ip.nrow();
  int model = ip.ncol()-1;
  NumericVector a(num_of_items), b(num_of_items), c(num_of_items), d(num_of_items), D(num_of_items);
  NumericMatrix output(num_of_theta, num_of_items);
  for(int j = 0; j < num_of_items; j++)
  {
    if (model == 1) {
      a[j] = 1;
    } else {
      a[j] = ip(j, 0);
    }
    if (model == 1) {
      b[j] = ip(j, 0);
    } else {
      b[j] = ip(j, 1);
    }
    if ((model == 1) || (model == 2)) {
      c[j] = 0;
    } else {
      c[j] = ip(j, 2);
    }
    if ((model == 1) || (model == 2) || (model == 3)) {
      d[j] = 1;
    } else {
      d[j] = ip(j, 3);
    }
    D[j] = ip(j, model);
  }
   for(int i = 0; i < num_of_theta; i++)
   {
      for(int j = 0; j < num_of_items; j++)
      {
        output(i,j) = ((D[j] * a[j]) * (D[j] * a[j]) * (d[j] - c[j]) * (d[j] - c[j])) /
       ((c[j] + d[j] * exp(D[j] * a[j] * (theta[i] - b[j]))) *
       (1 - c[j] + (1-d[j]) * exp(D[j] * a[j] * (theta[i]-b[j]))) *
       (1 + exp(-D[j] * a[j] * (theta[i] - b[j]))) * 
       (1 + exp(-D[j] * a[j] * (theta[i] - b[j]))));
       // output(i,j) = (pow(D[j] * a[j], 2) * pow(d[j] - c[j], 2)) /
       // ((c[j] + d[j] * exp(D[j] * a[j] * (theta[i] - b[j]))) *
       // (1 - c[j] + (1-d[j]) * exp(D[j] * a[j] * (theta[i]-b[j]))) *
       // pow(1 + exp(-D[j] * a[j] * (theta[i] - b[j])), 2));
      }
   }
  if (tif == true)
  {
    NumericMatrix tifMatrix(num_of_theta, 1);
    for(int i = 0; i < num_of_theta; i++)
    {
     tifMatrix(i, 0) = 0;
      for(int j = 0; j < num_of_items; j++)
      {
        tifMatrix(i, 0) = tifMatrix(i, 0) + output(i,j);
      }
    }
   return tifMatrix;
  }
  return output;
}

//##############################################################################
//########################### info_grm_bare_cpp ################################
//##############################################################################

// [[Rcpp::export]]
double info_grm_bare_cpp(double theta, Rcpp::S4 item)
{
  // This function calculates the expected information for one item and one
  // theta for Graded Response Model.
  // Calculations use formula on Baker and Kim (2004) p. 226, Eq. 8.23.
  Rcpp::List parList = as<List>(item.slot("parameters"));
  // Item difficulty
  Rcpp::NumericVector b = parList["b"];
  // Item discrimination
  double a = as<double>(parList["a"]);
  double D = as<double>(parList["D"]);
  // Set the  number of choices
  int no_choices = b.size() + 1;
  double prob_cdf1, prob_cdf2;
  double info = 0;
  prob_cdf1 = 1;
  for(int i = 0; i < no_choices - 1; i++)
  {
    prob_cdf2 = 1 / (1 + exp(-D * a * (theta - b[i])));
    info = info + D*D * a*a * (prob_cdf1 * 
      (1-prob_cdf1) - prob_cdf2 * (1-prob_cdf2)) * 
      (prob_cdf1 * (1-prob_cdf1) - prob_cdf2 * (1-prob_cdf2)) / 
      (prob_cdf1-prob_cdf2);
    // info = info + pow(D, 2) * pow(a, 2) * pow(prob_cdf1 * (1-prob_cdf1) -
    //   prob_cdf2 * (1-prob_cdf2), 2) / (prob_cdf1-prob_cdf2);
    // Rprintf("%d: %f\n", i, info);
    prob_cdf1 = prob_cdf2;
  }
  info = info + D*D * a*a * (prob_cdf1 * (1-prob_cdf1) - 0 * (1-0)) * 
    (prob_cdf1 * (1-prob_cdf1) - 0 * (1-0)) / (prob_cdf1-0);
  // info = info + D*D * a*a * pow(prob_cdf1 * (1-prob_cdf1) - 0
  //                                             * (1-0), 2) / (prob_cdf1-0);
  return info;
}

//##############################################################################
//########################### info_gpcm_bare_cpp ###############################
//##############################################################################

// [[Rcpp::export]]
double info_gpcm_bare_cpp(double theta, Rcpp::S4 item)
{
  // This function calculates the expected information for one item and one
  // theta for Generalized Partial Credit Model.
  // Calculations use formula Donoghue (1994), p.299, Eq.4.
  Rcpp::List parList = as<List>(item.slot("parameters"));
  std::string model = as<std::string>(item.slot("model"));
  
  // Item discrimination, if PCM, set them to 1, else if GPCM get them
  double a = 1;
  double D = 1;
  // Item difficulty
  Rcpp::NumericVector b;
  //unsigned int no_choices;
  if (model == "GPCM2") {
    b = as<Rcpp::NumericVector>(parList["d"]);
  } else { // "GPCM" or "PCM" 
    b = as<Rcpp::NumericVector>(parList["b"]);
  }
  unsigned int no_choices = b.size() + 1;
  // if (model == "GPCM2") {
  //   //Rcpp::NumericVector d = as<Rcpp::NumericVector>(parList["d"]);
  //   //no_choices = d.size() + 1;
  //   b = as<Rcpp::NumericVector>(parList["d"]);
  //   no_choices = b.size() + 1;
  //   double b_loc = as<double>(parList["b"]);
  //   // b = clone(d);
  //   for (unsigned int i = 0; i < no_choices; i++)
  //     // b[i] = b_loc - d[i];
  //     b[i] = b_loc - b[i];
  // } else { // "GPCM" or "PCM"
  //   b = as<Rcpp::NumericVector>(parList["b"]);
  //   no_choices = b.size() + 1;
  // }
  
  if (model == "GPCM" || model == "GPCM2") {
    a = as<double>(parList["a"]);
    D = as<double>(parList["D"]);
  } 
  
  // // Item difficulty
  // Rcpp::NumericVector b = parList["b"];
  // // Item discrimination
  // double a = as<double>(parList["a"]);
  // double D = as<double>(parList["D"]);
  // // Set the  number of choices
  // int no_choices = b.size() + 1;
  
  Rcpp::NumericVector P = prob_gpcm_bare_cpp(theta, item, 0);
  double lambda1 = 0;
  double lambda2 = 0;
  for(unsigned int i = 0; i < no_choices; i++)
  {
    lambda1 = lambda1 + i * i * P[i];
    lambda2 = lambda2 + i * P[i];
  }
  return D * D * a * a * (lambda1 - lambda2 * lambda2);
}


//##############################################################################
//########################### info_item_bare_cpp ###############################
//##############################################################################

// [[Rcpp::export]]
double info_item_bare_cpp(double theta, Rcpp::S4 item, bool observed, 
                          double resp)
{
  // This function calculates the information of a single item for single theta.
  if (Rcpp::NumericVector::is_na(resp)) return NA_REAL;
  std::string model(as<std::string>(as<Rcpp::S4>(item).slot("model")));
  if (model == "GRM") {
    // TODO: for "GRM" both observed and expected information is the same.
    if (observed) {
      return info_grm_bare_cpp(theta, item);
    } else
      return info_grm_bare_cpp(theta, item);
  } else if (model == "GPCM" || model == "PCM" || model == "GPCM2") {
    if (observed) {
      return -1 * resp_loglik_bare_item_cpp(resp, theta, item, 2);
      //return -1 * resp_loglik_sd_gpcm_bare_cpp(resp, theta, item);
    } else
      return info_gpcm_bare_cpp(theta, item);
  } else {
    // This is default option, i.e. 1PM-4PM, if there are other models
    // add them using 'else if' above.
    if (observed) {
      return -1 * resp_loglik_bare_item_cpp(resp, theta, item, 2);
      // return -1 * resp_loglik_sd_bare_cpp(resp, theta, item);
    } else
      return info_4pm_bare_cpp(theta, item);
  }
  return 0;
}


//##############################################################################
//########################### info_testlet_bare_cpp ############################
//##############################################################################

// [[Rcpp::export]]
double info_testlet_bare_cpp(
  double theta, Rcpp::S4 testlet, bool observed,
  Rcpp::Nullable<Rcpp::NumericVector> resp = R_NilValue) {
  // This function calculates the information of a single testlet for single
  // theta. It returns the total information of all items in the testlet.
  Rcpp::List item_list = as<List>(testlet.slot("item_list"));
  int num_of_items = item_list.size();
  S4 item;
  double output = 0;
  bool return_na = !Rf_isNull(resp);
  // resp is also used to calculate whether to include an item in information
  // calculation. It is not just for 'observed' information.
  if (return_na) {
    // NumericVector resp_i = as<NumericVector>(resp);
    if (as<NumericVector>(resp).size() != num_of_items)
      throw std::range_error("Inadmissible 'resp' value. The length of the "
             "'resp' and number of items in the testlet should be the same.");
  }
  for (int i = 0; i < num_of_items; i++) {
    item = as<Rcpp::S4>(item_list(i));
    if (Rf_isNull(resp)  || !NumericVector::is_na(as<NumericVector>(resp)[i])) {
      //Rcout << "resp_i[i] is not null. i = " << i << " and  resp[i] = " << resp_i[i] << std::endl;
      output = output + info_item_bare_cpp(theta, item, false, 0);
      return_na = false;
    }
  }
  if (return_na) return NA_REAL;
  return output;
}


//##############################################################################
//########################### info_item_cpp ####################################
//##############################################################################

// [[Rcpp::export]]
Rcpp::NumericVector info_item_cpp(
    Rcpp::NumericVector theta, Rcpp::S4 item, bool observed,
    Rcpp::Nullable<Rcpp::NumericVector> resp = R_NilValue)
{
  // This function calculates the information of a single item for multiple
  // thetas.
  int num_of_theta = theta.size();
  Rcpp::NumericVector output(num_of_theta);
  // Define a new resp_ variable to deal with nullable nature of resp.
  // Rcpp::NumericVector resp_(num_of_theta);
  // Make sure the size of resp is equal to the size of theta.
  if (observed && resp.isNotNull()) {
    // If not null cast resp to the underlying type
    // See: https://stackoverflow.com/a/43391979/2275286
    Rcpp::NumericVector resp_(resp.get());
    if (num_of_theta != resp_.size()) {
      throw std::invalid_argument("The size of the 'resp' vector should be equal to the size of theta.");
    }
    for(int i = 0; i < num_of_theta; i++) {
      output[i] = info_item_bare_cpp(theta[i], item, true, resp_[i]);
    }
    return output;
  }

  for(int i = 0; i < num_of_theta; i++) {
    output[i] = info_item_bare_cpp(theta[i], item, false, 0);
  }
  return output;
}


//##############################################################################
//########################### info_itempool_bare_cpp ##########################
//##############################################################################

// [[Rcpp::export]]
Rcpp::NumericVector info_itempool_bare_cpp(
    double theta, Rcpp::S4 ip, bool tif, bool observed,
    Rcpp::Nullable<Rcpp::NumericVector> resp = R_NilValue)
{
  // This function calculates the information of multiple items for a single
  // thetas.
  Rcpp::List item_list = as<List>(ip.slot("item_list"));
  IntegerVector ip_size = get_itempool_size(ip);
  int num_of_items = ip_size["elements"];
  unsigned int testlet_size;
  S4 item;
  NumericVector output(num_of_items);
  // an index that tracts the column number to read for the resp vector

  if (Rf_isNull(resp)) {
    for (int i = 0; i < num_of_items; i++) {
      item = as<Rcpp::S4>(item_list(i));
      if (item.inherits("Item")) {
        output[i] = info_item_bare_cpp(theta, item, false, 0);
      } else if (item.inherits("Testlet")) {
        output[i] = info_testlet_bare_cpp(theta, item, false, R_NilValue);
      }
    }
  } else {  // Use resp to calculate info, if resp is NA return NA otherwise calculate info
    int resp_index = 0;
    for (int i = 0; i < num_of_items; i++) {
      item = as<Rcpp::S4>(item_list(i));
      if (item.inherits("Item")) {
        output[i] = info_item_bare_cpp(theta, item, false, as<NumericVector>(resp)[resp_index]);
        resp_index += 1;
      } else if (item.inherits("Testlet")) {
        testlet_size = as<List>(as<S4>(item.slot("item_list")).slot("item_list")).size();
        NumericVector resp_ = Rcpp::as<NumericVector>(resp)[Rcpp::Range(resp_index, resp_index + testlet_size - 1)];
        output[i] = info_testlet_bare_cpp(theta, item, false, resp_);
        resp_index += testlet_size;
      }
    }
  }
  if (tif) {
    NumericVector tif_output(1);
    for (int i = 0; i < num_of_items; i++) {
      if (!R_IsNA(output[i]))
        tif_output[0] = tif_output[0] + output[i];
    }
    return tif_output;
  }
  return output;
}


//##############################################################################
//########################### info_itempool_cpp ###############################
//##############################################################################

// [[Rcpp::export]]
Rcpp::NumericMatrix info_itempool_cpp(
  Rcpp::NumericVector theta, Rcpp::S4 ip, bool tif, bool observed,
  Rcpp::Nullable<Rcpp::NumericMatrix> resp = R_NilValue)
{
  // This function calculates the information of multiple items for multiple
  // thetas.
  int num_of_cols = as<List>(ip.slot("item_list")).size();
  int num_of_theta = theta.size();
  if (tif) num_of_cols = 1;
  NumericMatrix output(num_of_theta, num_of_cols);
  if (Rf_isNull(resp)) {
    for(int i = 0; i < num_of_theta; i++)
      output(i, Rcpp::_) = info_itempool_bare_cpp(theta[i], ip, tif, false,
                                                   R_NilValue);
  } else {
    NumericMatrix resp_ = Rcpp::as<NumericMatrix>(resp);
    NumericVector resp_row;
    for(int i = 0; i < num_of_theta; i++) {
      resp_row = resp_(i, Rcpp::_);
      output(i, Rcpp::_) = info_itempool_bare_cpp(theta[i], ip, tif, false,
                                                   resp_row);
    }
  }
  return output;
}
