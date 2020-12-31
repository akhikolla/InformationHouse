#include <Rcpp.h>
#include "itempool_class_methods.h"
using namespace Rcpp;


//##############################################################################
//##############################################################################
//########################### prob_cpp #########################################
//##############################################################################
//##############################################################################


//##############################################################################
//########################### prob_4pm_bare_cpp ################################
//##############################################################################

// [[Rcpp::export]]
double prob_4pm_bare_cpp(double theta, Rcpp::S4 item, int derivative = 0)
{
  // This function calculates the probability of correct response for one item
  // and one theta.
  Rcpp::S4 tempS4 = clone(item);
  Rcpp::List parList = as<List>(tempS4.slot("parameters"));
  std::string model = as<std::string>(tempS4.slot("model"));
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
  if (derivative == 0) {
    return c + (d-c) / (1+exp(-a * D * (theta - b)));
  } else if (derivative == 1) { // First Derivative
    return ((d-c) * D * a) / (exp(D * a * (theta - b)) + 2 +
                              exp(-D * a * (theta - b)));
  } else if (derivative == 2) { // Second Derivative
    return -(d-c) * D * D * a * a *
      (exp(D * a * (theta - b)) - exp(-D * a * (theta - b))) /
      ((exp(D * a * (theta - b)) + 2 + exp(-D * a * (theta - b))) * 
      (exp(D * a * (theta - b)) + 2 + exp(-D * a * (theta - b))));
      // pow((exp(D * a * (theta - b)) + 2 + exp(-D * a * (theta - b))), 2);
  } else
    stop("'derivative' value can take only values 0, 1 or 2.");
}

//##############################################################################
//########################### prob_4pm_item_cpp ################################
//##############################################################################

// [[Rcpp::export]]
Rcpp::NumericVector prob_4pm_item_cpp(NumericVector theta, Rcpp::S4 item,
                                      int derivative = 0)
{
  // This function calculates the probability of correct response for one item
  // and multiple thetas.;
  unsigned int num_of_theta = theta.size();
  Rcpp::NumericVector output(num_of_theta);
  for(unsigned int i = 0; i < num_of_theta; i++) {
    output[i] = prob_4pm_bare_cpp(theta[i], item, derivative);
  }
  return output;
}

//##############################################################################
//########################### prob_4pm_itempool_cpp ############################
//##############################################################################

// [[Rcpp::export]]
Rcpp::NumericMatrix prob_4pm_itempool_cpp(NumericVector theta, Rcpp::S4 ip,
                                          int derivative = 0)
{
  // This function calculates the probability of correct response for multiple
  // items and multiple thetas.;
  Rcpp::List item_list = as<List>(ip.slot("item_list"));
  unsigned int num_of_items = item_list.size();
  unsigned int num_of_theta = theta.size();
  NumericMatrix output(num_of_theta, num_of_items);
  for(unsigned int i = 0; i < num_of_theta; i++)
  {
    for(unsigned int j = 0; j < num_of_items; j++)
    {
      output(i,j) = prob_4pm_bare_cpp(theta(i), as<Rcpp::S4>(item_list(j)),
             derivative);
    }
  }
  return output;
}

//##############################################################################
//########################### prob_grm_bare_cpp ################################
//##############################################################################

// [[Rcpp::export]]
Rcpp::NumericVector prob_grm_bare_cpp(double theta, Rcpp::S4 item,
                                      int derivative = 0)
{
  // This function calculates the probability of correct response for one item
  // and one theta for Graded Response Model.
  //
  // This function calculates the probability of correct response for one item
  // and one theta for Graded Response Model.
  // Based on the first equation on page 219 of Baker and Kim (2004)
  Rcpp::List parList = as<List>(item.slot("parameters"));
  // Item difficulty
  Rcpp::NumericVector b = parList["b"];
  // Item discrimination
  double a = as<double>(parList["a"]);
  double D = as<double>(parList["D"]);
  // Set the  number of choices
  unsigned int no_choices = b.size() + 1;
  if (derivative == 0) {
    // Vector holding the probabilities of each response
    Rcpp::NumericVector probs(no_choices);
    // prob_cdf1 is P*(k) and prob_cdf2 is P*(k+1)
    double prob_cdf1, prob_cdf2;
    prob_cdf1 = 1;
    for(unsigned int i = 0; i < no_choices - 1; i++)
    {
      prob_cdf2 = 1 / (1 + exp(-D * a * (theta - b[i])));
      // Rprintf("1: %f\n", prob_cdf1);
      // Rprintf("2: %f\n", prob_cdf2);
      probs[i] = prob_cdf1 - prob_cdf2;
      prob_cdf1 = prob_cdf2;
    }
    probs[no_choices - 1] = prob_cdf1;
    return probs;
  } else if (derivative == 1) {
    // Vector holding the first derivatives of each response
    Rcpp::NumericVector fd(no_choices);
    // prob_cdf1 is P*(k) and prob_cdf2 is P*(k+1)
    double prob_cdf1, prob_cdf2;
    prob_cdf1 = 1;
    for(unsigned int i = 0; i < no_choices - 1; i++)
    {
      prob_cdf2 = 1 / (1 + exp(-D * a * (theta - b[i])));
      // Rprintf("1: %f\n", prob_cdf1);
      // Rprintf("2: %f\n", prob_cdf2);
      fd[i] = a * D * (prob_cdf1 * (1 - prob_cdf1) - prob_cdf2 * (1 - prob_cdf2));
      prob_cdf1 = prob_cdf2;
    }
    //  Since prob_cdf2 = 0, terms on the right disappeared.
    fd[no_choices - 1] = a * D * prob_cdf1 * (1 - prob_cdf1);
    return fd;
  } else
    stop("'derivative' value for 'GRM' can take only values 0 or 1.");
}


//##############################################################################
//########################### prob_gpcm_bare_cpp ###############################
//##############################################################################

// [[Rcpp::export]]
Rcpp::NumericVector prob_gpcm_bare_cpp(double theta, Rcpp::S4 item,
                                       int derivative = 0)
{
  // This function calculates the probability of correct response for one item
  // and one theta for Partial Credit Model and Generalized Partial Credit
  // Model.
  //
  // This function calculates the first derivative of probability of a response
  // for one item and one theta for Generalized Partial Credit Model.
  // Function is based on Donoghue (1994), p.309, Eq.2.

  // This function calculates the second derivative of probability of a response
  // for one item and one theta for Generalized Partial Credit Model.
  // Function is based on Donoghue (1994), p.309, Eq.3.

  Rcpp::List parList = as<List>(item.slot("parameters"));
  std::string model = as<std::string>(item.slot("model"));
  
  // Item discrimination, if PCM, set them to 1, else if GPCM get them
  double a = 1;
  double D = 1;
  
  Rcpp::NumericVector b; 

  if (model == "GPCM2") {      
    b = as<double>(parList["b"]) - as<Rcpp::NumericVector>(parList["d"]);
  } else {
    b = as<Rcpp::NumericVector>(parList["b"]);
  }      
  
  unsigned int no_choices = b.size() + 1; // Number of categories
  
  if (model == "GPCM" || model == "GPCM2") {
    a = as<double>(parList["a"]);
    D = as<double>(parList["D"]);
  } 
  
  if (derivative == 0) {
    // Vector holding the probabilities of each response
    Rcpp::NumericVector probs(no_choices);
    // Vector holding the numerator of the probability function
    Rcpp::NumericVector numerator(no_choices);
    double denominator; // The denominator of the probability function
    numerator[0] = 1; // The exponent of the first category (i.e. exp(0) = 1)             
    for (unsigned int i = 1; i < no_choices; i++) {
      for (unsigned int j = 0; j < i; j++) {
        numerator[i] = numerator[i] + D * a * (theta - b[j]);  // cumulative sum
      }
      numerator[i] = exp(numerator[i]);
    }    
    // Calculate the denominator
    denominator = sum(numerator);
    // Calculate the probability of each category
    probs = numerator/denominator;
    return probs;    
  } else if (derivative == 1) {
    // Vector holding probabilities of each response catergory
    Rcpp::NumericVector P = prob_gpcm_bare_cpp(theta, item);
    // Vector holding the first derivative of probabilities of each response
    Rcpp::NumericVector fd(no_choices);
    double lambda1 = 0;
    for (unsigned int i = 0; i < no_choices; i++)
    {
      lambda1 = lambda1 + i * P[i];
    }
    for (unsigned int i = 0; i < no_choices; i++)
    {
      fd[i] = D * a * P[i] * (i - lambda1);
    }
    return fd;
  } else if (derivative == 2) {
    // Vector holding probabilities of each response catergory
    Rcpp::NumericVector P = prob_gpcm_bare_cpp(theta, item);
    // Vector holding the second derivative of probabilities of each response
    Rcpp::NumericVector sd(no_choices);
    double lambda1 = 0;
    double lambda2 = 0;
    for (unsigned int i = 0; i < no_choices; i++)
    {
      lambda1 = lambda1 + i * P[i];
      lambda2 = lambda2 + i*i * P[i];
    }
    for (unsigned int i = 0; i < no_choices; i++)
    {
      sd[i] = D * D * a * a * P[i] * (i * i - 2 * i * lambda1 +
        2 * lambda1 * lambda1 - lambda2);
    }
    return sd;
  } else
    stop("'derivative' value can take only values 0, 1 or 2.");
}



//##############################################################################
//########################### prob_poly_bare_cpp ###############################
//##############################################################################

// [[Rcpp::export]]
Rcpp::NumericVector prob_poly_bare_cpp(double theta, Rcpp::S4 item,
                                       int derivative = 0,
                                       bool expected_value = false) {
  std::string model = as<std::string>(item.slot("model"));
  Rcpp::NumericVector result;
  if (model == "GPCM" || model == "PCM" || model == "GPCM2") {
    result = prob_gpcm_bare_cpp(theta, item, derivative);
  } else if (model == "GRM") {
    result = prob_grm_bare_cpp(theta, item, derivative);
  } else
    stop("This model has not been implemented in 'prob_poly_bare_cpp()' function yet.");
  if (expected_value) {
    Rcpp::NumericVector temp_result(1);
    for(int i = 0; i < result.size(); i++) {
      temp_result[0] = temp_result[0] + i * result[i]; 
    }
    return temp_result;
  } 
  return result;
}

//##############################################################################
//########################### prob_mirt_bare_cpp ###############################
//##############################################################################

// [[Rcpp::export]]
double prob_mirt_bare_cpp(Rcpp::NumericVector theta, Rcpp::S4 item,
                          int derivative = 0)
{
  Rcpp::List parList = as<List>(item.slot("parameters"));
  std::string model = as<std::string>(item.slot("model"));
  unsigned int num_of_theta = theta.size();
  std::vector<double> a(num_of_theta, 0);
  double d = as<double>(parList["d"]);
  double c = 0;
  double D = as<double>(parList["D"]);
  double upperAsymptote = 1;
  if (model != "M1PL")
  {
    std::vector<double> aTemp = parList["a"];
    if (aTemp.size() != num_of_theta) {
      stop("Number of theta's should be equal to the number of 'a' parameters");
    }
    a = aTemp;
  }
  if (model == "M3PL" || model == "M4PL") {
    c = as<double>(parList["c"]);
  }
  if (model == "M4PL") {
    upperAsymptote = as<double>(parList["upperAsymptote"]);
  }
  double exponent = d;
  for (unsigned int i = 0; i < num_of_theta; i++)
  {
    exponent = exponent + a[i] * theta[i];
  }
  return c + (upperAsymptote - c) / (1 + exp(-D * exponent));
}

//##############################################################################
//########################### prob_mirt_item_cpp ###############################
//##############################################################################

// [[Rcpp::export]]
Rcpp::NumericVector prob_mirt_item_cpp(Rcpp::NumericMatrix theta,
                                       Rcpp::S4 item,
                                       int derivative = 0)
{
   unsigned int num_of_theta = theta.nrow();
  unsigned int noDim = theta.ncol();
  Rcpp::NumericVector singleTheta(noDim);
  Rcpp::NumericVector output(num_of_theta);
  for(unsigned int i = 0; i < num_of_theta; i++) {
    for(unsigned int j = 0; j < noDim; j++) {
      singleTheta[j] = theta(i,j);
    }
    output[i] = prob_mirt_bare_cpp(singleTheta, item);
  }
  return output;
}

//##############################################################################
//########################### prob_mirt_itempool_cpp ###########################
//##############################################################################

// [[Rcpp::export]]
Rcpp::NumericMatrix prob_mirt_itempool_cpp(Rcpp::NumericMatrix theta,
                                           Rcpp::S4 ip,
                                           int derivative = 0)
{
  Rcpp::List item_list = as<List>(ip.slot("item_list"));
  unsigned int num_of_items = item_list.size();
  unsigned int num_of_theta = theta.nrow();
  unsigned int noDim = theta.ncol();
  Rcpp::NumericVector singleTheta(noDim);
  Rcpp::NumericMatrix output(num_of_theta, num_of_items);
  for(unsigned int i = 0; i < num_of_theta; i++)
  {
    for(unsigned int j = 0; j < num_of_items; j++)
    {
      for(unsigned int k = 0; k < noDim; k++) {
        singleTheta[k] = theta(i,k);
      }
      output(i,j) = prob_mirt_bare_cpp(singleTheta, as<Rcpp::S4>(item_list(j)));
    }
  }
  return output;
}


//##############################################################################
//########################### prob_bare_item_cpp ###############################
//##############################################################################
// [[Rcpp::export]]
Rcpp::NumericVector prob_bare_item_cpp(Rcpp::NumericVector theta, Rcpp::S4 item,  
                                       int derivative = 0, 
                                       bool expected_value = false) {
  std::string model = as<std::string>(item.slot("model"));
  Rcpp::NumericVector result;
  if (model == "GPCM" || model == "PCM" || model == "GRM" || model == "GPCM2") {
    result = prob_poly_bare_cpp(theta[0], item, derivative, expected_value);
  } else if (model == "M1PL" || model == "M2PL" || model == "M3PL") {
    result = prob_mirt_bare_cpp(theta, item, derivative);    
  } else if (model == "Rasch" || model == "1PL" || model == "2PL" ||
             model == "3PL" || model == "4PL") {
    // Rcout << "UIRT model " << theta[0] << " - " << std::endl;
    result = NumericVector::create(prob_4pm_bare_cpp(theta[0], item, 
                                                     derivative));
  } else
    stop("This model has not been implemented in 'prob()' function yet.");
  return result;
}


//##############################################################################
//########################### prob_bare_itempool_cpp ###########################
//##############################################################################
//' This function calculates the probability of each response option for one
//' theta value. It returns a matrix, each column represents a response
//' option and each row represents an item.
//'
//' @noRd
// [[Rcpp::export]]
Rcpp::NumericMatrix prob_bare_itempool_cpp(Rcpp::NumericVector theta, 
                                           Rcpp::S4 ip, int derivative = 0, 
                                           bool expected_value = false)
{
  Rcpp::List item_list = flatten_itempool_cpp(ip);
  Rcpp::NumericVector max_scores = get_max_possible_score_itempool_cpp(ip);
  int max_score = Rcpp::max(max_scores);
  if (expected_value) 
    max_score = 0;
  int num_of_sa_items = item_list.size();
  Rcpp::NumericMatrix output(num_of_sa_items, max_score + 1);
  Rcpp::S4 item("Item");
  Rcpp::NumericVector temp_prob;
  double num_of_probs;
  if (!expected_value)
    std::fill( output.begin(), output.end(), NumericVector::get_na() );
  for(int i = 0; i < num_of_sa_items; i++) {
    item = as<Rcpp::S4>(item_list[i]);
    temp_prob = prob_bare_item_cpp(theta, item, derivative, expected_value);
    if (expected_value) { 
      output(i, 0) = temp_prob[0];
    } else {
      num_of_probs = temp_prob.size();
      if (num_of_probs == 1) {
        output(i, 0) = 1 - temp_prob[0];
        output(i, 1) = temp_prob[0];
      } else if (num_of_probs > 1) {
        for (int j = 0; j < num_of_probs; j++) {
          output(i, j) = temp_prob[j];
        }
      }
    }
  }
  return output;
}
