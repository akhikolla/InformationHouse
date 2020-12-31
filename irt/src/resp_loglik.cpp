#include <Rcpp.h>
#include "itempool_class_methods.h"
#include "prob.h"
#include "resp_lik.h"
using namespace Rcpp;


//#############################################################################@
//#############################################################################@
//########################### resp_loglik ######################################
//#############################################################################@
//#############################################################################@

//##############################################################################
//########################### resp_loglik_bare_item_cpp ########################
//##############################################################################
// [[Rcpp::export]]
double resp_loglik_bare_item_cpp(double resp, double theta, Rcpp::S4 item,
                                 int derivative = 0)
{
  // Calculate response log-likelihood for one item and one theta's (and
  // responses)

  // Deal with missing responses, return NA directly
  if (NumericVector::is_na(resp))
    return NA_REAL;


  if (derivative == 0) { // Response log-likelihood
      double resp_lik = resp_lik_bare_item_cpp(resp, theta, item);
      // Here is an assumption that resp_lik is not NA. It can be if the model
      // is not implemented yet in "resp_lik_bare_cpp" function.
      return log(resp_lik);
      // // Get the Psychometric Model name
      // std::string model = as<std::string>(item.slot("model"));
      // if (model == "GPCM" || model == "GRM") {
      //   Rcpp::NumericVector P(resp);
      //   if (model == "GPCM") {
      //     P = prob_gpcm_bare_cpp(theta, item);
      //   } else if (model == "GRM") {
      //     P = prob_grm_bare_cpp(theta, item);
      //   }
      //   return log(P[resp]);
      // } else if (model == "1PL" || model == "2PL" || model == "3PL" ||
      //   model == "4PL") {
      //   double P = prob_4pm_bare_cpp(theta, item);
      //   return resp * log(P) + (1 - resp) * log(1 - P);
      // }
      // // This should return a more sensible value if model is none of the above.
      // return 0;
  } else if (derivative == 1) { // First derivative of response log-likelihood
      std::string model = as<std::string>(item.slot("model"));
      if (model == "GPCM" || model == "PCM" || model == "GPCM2") {
        // This function calculates the first derivative for one item and one
        // theta for Generalized Partial Credit Model.
        // Calculations use formula Penfield and Bergeron (2005), p. 220, Eq. (2)
        Rcpp::List parList = as<List>(item.slot("parameters"));
        std::string model = as<std::string>(item.slot("model"));
        
        // Item discrimination, if PCM, set them to 1, else if GPCM get them
        double a = 1;
        double D = 1;
        // Item difficulty
        Rcpp::NumericVector b;
        unsigned int no_choices;
        if (model == "GPCM2") {
          Rcpp::NumericVector d = as<Rcpp::NumericVector>(parList["d"]);
          no_choices = d.size() + 1;
          double b_loc = as<double>(parList["b"]);
          b = clone(d);
          for (unsigned int i = 0; i < no_choices; i++)
            b[i] = b_loc - d[i];
        } else { // "GPCM" or "PCM"
          b = as<Rcpp::NumericVector>(parList["b"]);
          no_choices = b.size() + 1;
        }        
        if (model == "GPCM" || model == "GPCM2") {
          a = as<double>(parList["a"]);
          D = as<double>(parList["D"]);
        }        
        
        Rcpp::NumericVector P = prob_gpcm_bare_cpp(theta, item, 0);
        double lambda1 = 0;
        for(unsigned int i = 0; i < no_choices; i++) {
          lambda1 = lambda1 + i * P[i];
        }
        double result = 0;
        for(unsigned int i = 0; i < no_choices; i++) {
          // It is assumed that scoring function U_j and j are the same.
          result = result + i * D * a * (i - lambda1);
        }
        return result;
      } else if (model == "Rasch" || model == "1PL" || model == "2PL" ||
                 model == "3PL" || model == "4PL") {
        double P = prob_4pm_bare_cpp(theta, item, 0);
        double dP =prob_4pm_bare_cpp(theta, item, 1);
        return  dP * (resp - P) / (P * (1 - P));
      }
      return NA_REAL;
  } else if (derivative == 2) { // Second derivative of response log-likelihood
      // Get the Psychometric Model name
      std::string model = as<std::string>(item.slot("model"));
      if (model == "GPCM" || model == "GPCM2" || model == "PCM") {
        // This function calculates the second derivative for one item and one
        // theta for Generalized Partial Credit Model.
        // Calculations use formula Penfield and Bergeron (2005), p. 220, Eq. (3)
        // Alternative formula can be found in Donoghue (1994), p.309, second
        // equation from the last
        
        Rcpp::List parList = as<List>(item.slot("parameters"));
        std::string model = as<std::string>(item.slot("model"));
        
        // Item discrimination, if PCM, set them to 1, else if GPCM get them
        double a = 1;
        double D = 1;
        // Item difficulty
        Rcpp::NumericVector b;
        unsigned int no_choices;
        if (model == "GPCM2") {
          Rcpp::NumericVector d = as<Rcpp::NumericVector>(parList["d"]);
          no_choices = d.size() + 1;
          double b_loc = as<double>(parList["b"]);
          b = clone(d);
          for (unsigned int i = 0; i < no_choices; i++)
            b[i] = b_loc - d[i];
        } else { // "GPCM" or "PCM"
          b = as<Rcpp::NumericVector>(parList["b"]);
          no_choices = b.size() + 1;
        }        
        if (model == "GPCM" || model == "GPCM2") {
          a = as<double>(parList["a"]);
          D = as<double>(parList["D"]);
        } 
        
        Rcpp::NumericVector P = prob_gpcm_bare_cpp(theta, item, 0);
        double lambda1 = 0;
        double lambda2 = 0;
        for(unsigned int i = 0; i < no_choices; i++) {
          lambda1 = lambda1 + i * P[i];
          lambda2 = lambda2 + i*i * P[i];
        }
        return D*D * a*a * (lambda1*lambda1 - lambda2);
      } else if (model == "Rasch" || model == "1PL" || model == "2PL" ||
                 model == "3PL" || model == "4PL") {
        // This function calculates the second derivative of the log likelihood of
        // a response string for one theta.
        double P = prob_4pm_bare_cpp(theta, item, 0);
        double dP = prob_4pm_bare_cpp(theta, item, 1);
        double d2P = prob_4pm_bare_cpp(theta, item, 2);
        return  (1 / (P * (1 - P))) * (d2P * (resp - P) - dP * dP *
          (1+(resp - P) * (1- 2*P) / (P * (1-P)) ) );
        // return  pow(P * (1 - P), -1) * (d2P * (resp - P) - pow(dP, 2) *
        //   (1+(resp - P) * (1- 2*P) / (P * (1-P)) ) );
      }
      return NA_REAL;
  } else
      stop("'derivative' value can take only values 0, 1 or 2.");
}



//##############################################################################
//########################### resp_loglik_item_cpp #############################
//##############################################################################
//' Calculate the response log-likelihood of a response string.
//' @param resp Response vector.
//' @param theta Theta value.
//' @param item An \code{Item-class} object.
//' @param derivative An integer indicating which derivative to calculate:
//'    0 = No derivative
//'    1 = First derivative
//'    2 = Second derivative
//' 
//' @noRd
//' 
// [[Rcpp::export]]
Rcpp::NumericVector resp_loglik_item_cpp(Rcpp::NumericVector resp,
                                         Rcpp::NumericVector theta,
                                         Rcpp::S4 item,
                                         int derivative = 0)
{
  // Calculate response log-likelihood for one item and multiple theta's (and
  // responses)
  unsigned int num_of_theta = theta.size();
  Rcpp::NumericVector output(num_of_theta);
  for(unsigned int i = 0; i < num_of_theta; i++)
    output[i] = resp_loglik_bare_item_cpp(resp[i], theta[i], item, derivative);
  return output;
}


//##############################################################################
//########################### resp_loglik_bare_testlet_cpp #####################
//##############################################################################
//' Calculate response log-likelihood for a testlet and a single theta (and a
//' response string)
//' @param resp Response vector.
//' @param theta Theta value.
//' @param testlet A \code{Testlet-class} object.
//' @param derivative An integer indicating which derivative to calculate:
//'    0 = No derivative
//'    1 = First derivative
//'    2 = Second derivative
//' 
//' @noRd
//' 
// [[Rcpp::export]]
double resp_loglik_bare_testlet_cpp(Rcpp::NumericVector resp, double theta,
                                    Rcpp::S4 testlet, int derivative = 0)
{
  // Calculate response log-likelihood for a testlet and single examinee
  Rcpp::List item_list = as<List>(testlet.slot("item_list"));
  unsigned int num_of_items = item_list.size();
  double output = 0;
  for(unsigned int i = 0; i < num_of_items; i++) {
    if (!NumericVector::is_na(resp[i]))
      output = output + resp_loglik_bare_item_cpp(
            resp[i], theta, as<Rcpp::S4>(item_list(i)), derivative);
  }
  return output;
}



//##############################################################################
//########################### resp_loglik_testlet_cpp ##########################
//##############################################################################
// [[Rcpp::export]]
Rcpp::NumericVector resp_loglik_testlet_cpp(Rcpp::NumericMatrix resp,
                                            Rcpp::NumericVector theta,
                                            Rcpp::S4 testlet,
                                            int derivative = 0)
{
  // Calculate response log-likelihood for an Itempool and multiple
  // theta's (and response strings)
  unsigned int num_of_theta = theta.size();
  NumericVector output(num_of_theta);
  for(unsigned int i = 0; i < num_of_theta; i++) {
    // Get the row belong to the examinee. It is assumed that each row represents
    // an examinee.
    NumericVector resp_vector = resp(i, _);
    output[i] = resp_loglik_bare_testlet_cpp(resp_vector, theta[i], testlet,
                                             derivative);
  }
  return output;
}



//##############################################################################
//########################### resp_loglik_bare_itempool_cpp ###################
//##############################################################################
// [[Rcpp::export]]
double resp_loglik_bare_itempool_cpp(Rcpp::NumericVector resp, double theta,
                                      Rcpp::S4 ip, int derivative = 0)
{
  // Calculate response log-likelihood for an Itempool and single examinee
  Rcpp::List item_list = as<List>(ip.slot("item_list"));
  unsigned int num_of_items = item_list.size();
  // The variable that will hold the number of items in a testlet:
  unsigned int testlet_size;
  double output = 0;
  Rcpp::S4 item; // This can be item or a testlet
  // an index that tracks the column number to read for the resp vector
  int resp_index = 0;
  for(unsigned int i = 0; i < num_of_items; i++) {
    item = as<Rcpp::S4>(item_list(i));
    if (item.inherits("Testlet")) {
      // Find the number of items within the testlet
      testlet_size = as<List>(as<S4>(item.slot("item_list")).slot("item_list")).size();
      output = output + resp_loglik_bare_testlet_cpp(
        resp[Rcpp::Range(resp_index, resp_index + testlet_size - 1)], theta,
        item, derivative);
      resp_index += testlet_size;
    } else if (item.inherits("Item")) {
      // Check the class of the item, if it is "Item"
      if (!NumericVector::is_na(resp[resp_index]))
        output = output + resp_loglik_bare_item_cpp(resp(resp_index), theta,
                                                    item, derivative);
      resp_index += 1;
    }
  }
  return output;
}


//##############################################################################
//########################### resp_loglik_itempool_cpp ########################
//##############################################################################
// [[Rcpp::export]]
Rcpp::NumericVector resp_loglik_itempool_cpp(Rcpp::NumericMatrix resp,
                                              Rcpp::NumericVector theta,
                                              Rcpp::S4 ip, int derivative = 0)
{
  // Calculate response log-likelihood for an Itempool and multiple
  // theta's (and response strings)
  unsigned int num_of_theta = theta.size();
  NumericVector output(num_of_theta);
  for(unsigned int i = 0; i < num_of_theta; i++) {
    // Get the row belong to the examinee. It is assumed that each row represents
    // an examinee.
    NumericVector resp_vector = resp(i, _);
    output[i] = resp_loglik_bare_itempool_cpp(resp_vector, theta[i], ip,
                                               derivative);
  }
  return output;
}

