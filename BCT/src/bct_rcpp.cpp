#include <Rcpp.h>
#include <string>
#include "utils.h"

using namespace Rcpp;
//' @importFrom Rcpp evalCpp
//' @exportPattern "^[[:alpha:]]+"


//' @title k-Bayesian Context Trees (kBCT) algorithm
//' @description Finds the top k a posteriori most likely tree models.
//'
//' @param input_data the sequence to be analysed. 
//' The sequence needs to be a "character" object. See the examples section on how to transform any dataset to a "character" object.
//' @param depth maximum memory length.
//' @param beta hyper-parameter of the model prior. 
//' Takes values between 0 and 1. If not initialised in the call function, the default value is \ifelse{html}{\out{1-2<sup>-m+1</sup>}}{\eqn{1 - 2^{-m+1}}}, 
//' where \ifelse{html}{\out{m}}{\eqn{m}} is the size of the alphabet; for more information see: \href{https://arxiv.org/pdf/2007.14900.pdf}{Kontoyiannis et al. (2020)}.
//' @param k number of the a posteriori most likely tree models to be identified.
//' 
//' @return a list object which includes:
//' \item{Contexts}{top k a posteriori most likely models. Each model given as a list object containing the contexts of its leaves.}
//' \item{Results}{a dataframe with the following columns: prior probability, log(prior probability), posterior probability, log(posterior probability), posterior odds, number of leaves, maximum depth, BIC score, AIC score and maximum log-likelihood. }
//'
//' 
//' @export
//' @seealso \code{\link{BCT}}
//' 
//' @examples
//' # Finding the first 5 a posteriori most likely models with maximum depth <= 5 
//' # for the SP500 dataset (with default value beta):
//' 
//' kBCT(SP500, 5, 2)  
//' 
//' # For custom beta (e.g. 0.8):
//' 
//' kBCT(SP500, 5, 2, 0.8)  
//' 
//' # The type of the input dataset is "character"
//' # If the dataset is contained within a vector:
//' 
//' q <- c(1, 0, 1, 0, 0, 1, 0, 0, 1, 1, 1, 0, 1, 0)
//' 
//' # Convert a vector to a "character" object:
//' s <- paste(q, collapse = "")
//'
//' kBCT(s, 2, 2)
//' 
//' # Reading a file using the readChar function 
//' 
//' # Example 1: The dataset is stored in a .txt file
//' 
//' # fileName <- '~/example_data.txt' # fileName stores the path to the dataset
//' 
//' # s<-readChar(fileName, file.info(fileName)$size)
//' 
//' # Make sure that s does not contain any "\n" at the end of the string
//' # To remove last entry:
//' # s<-gsub('.$', '', s)
//' 
//' # To remove any unwanted characters (e.g. "\n"):
//' # s<-gsub('\n', '', s)
//' 
//' # Example 2: The dataset is stored in a .csv file
//' 
//' # fileName <- '~/example_data.csv' # fileName stores the path to the dataset
//' 
//' # s<-readChar(fileName, file.info(fileName)$size)
//' 
//' # Depending on the running environment, 
//' # s might contain unwanted characters such as: "\n" or "\r\n".
//' # Remove any unwanted characters (e.g. "\r\n"):
//' # s<-gsub('\r\n', '', s)
//' 
//' # Always make sure that s does not contain any unwanted characters
//' 
// [[Rcpp::export]]
List kBCT(CharacterVector input_data, IntegerVector depth, IntegerVector k, Nullable<NumericVector> beta = R_NilValue){
  int D;  // maximum depth
  unsigned int k_max; // number of desired models
  double b;  // hyper-parameter used in prior
  k_max = k[0];
  D = depth[0];
  std::string s = Rcpp::as<std::string>(input_data);
  
  
  if(beta.isNotNull()){
    NumericVector beta_(beta);
    b = beta_[0];
    set_global_parameters(s, D, k_max, b); 
  }
  else
    set_global_parameters(s, D, k_max);
  
  vector <Tree_properties> tp_vec = build_kbct(); // get the properties of the maximum a posteriori models
  List lst;
  
  NumericVector prior;
  NumericVector log_prior;
  NumericVector posterior;
  NumericVector log_posterior;
  NumericVector odd_posterior;
  NumericVector n_leaves;
  NumericVector max_depth;
  NumericVector bic;
  NumericVector aic;
  NumericVector mle;
  
  for(unsigned int i = 0; i<k_max; i++){
    StringVector lst_context;
    
    for(unsigned int j = 0; j<tp_vec[i].context.size(); j++){
      lst_context.push_back(tp_vec[i].context[j]);
    }
    
    // append the results to dataframe
    prior.push_back(tp_vec[i].prior);
    log_prior.push_back(tp_vec[i].log_prior * log(2));
    posterior.push_back(tp_vec[i].posterior);
    log_posterior.push_back(tp_vec[i].log_posterior * log(2));
    odd_posterior.push_back(tp_vec[i].odd_posterior);
    n_leaves.push_back(tp_vec[i].n_leaves);
    max_depth.push_back(tp_vec[i].max_depth);
    bic.push_back(tp_vec[i].bic);
    aic.push_back(tp_vec[i].aic);
    mle.push_back(tp_vec[i].mle);
    
    lst.push_back(lst_context);
  }
  
  DataFrame df = DataFrame::create( _["prior"] = prior,
                                    _["log_prior"] = log_prior,
                                    _["posterior"] = posterior,
                                    _["log_posterior"] = log_posterior,
                                    _["posterior_odds"] = odd_posterior,
                                    _["number_leaves"] = n_leaves,
                                    _["max_depth"] = max_depth,
                                    _["BIC"] = bic,
                                    _["AIC"] = aic,
                                    _["max_log_lik"] = mle);
  
  List results = List::create(Named("Contexts") = lst, _["Results"] = df);
  
  return results;
}

//==================================================================================

//' @title Calculating the log-loss incurred in prediction
//' @description Compute the log-loss incurred in BCT prediction with memory length D. Given an initial context
//' \ifelse{html}{\out{(x<sub>-D+1</sub>, ..., x<sub>0</sub>)}}{\eqn{(x_{-D+1},...,x_0)}} and training data \ifelse{html}{\out{(x<sub>1</sub>, ..., x<sub>n</sub>)}}{(\eqn{x_{1},...,x_n)}}, the log-loss is computed in sequentially predicting
//' the test data \ifelse{html}{\out{(x<sub>n+1</sub>, ..., x<sub>n+T</sub>)}}{\eqn{(x_{n+1},...,x_{n+T})}}. The function outputs the cummulative, normalized (per-sample) log-loss, at each prediction step; for more information see \href{https://arxiv.org/pdf/2007.14900.pdf}{Kontoyiannis et al.(2020)}.
//' @param input_data the sequence to be analysed. 
//' The sequence needs to be a "character" object. See the examples section of BCT/kBCT functions on how to transform any dataset to a "character" object.
//' @param depth maximum memory length.
//' @param train_size number of samples used in the training set. The training set size should be at least equal to the depth.
//' @param beta hyper-parameter of the model prior. 
//' Takes values between 0 and 1. If not initialised in the call function, the default value is \ifelse{html}{\out{1-2<sup>-m+1</sup>}}{\eqn{1 - 2^{-m+1}}}, 
//' where \ifelse{html}{\out{m}}{\eqn{m}} is the size of the alphabet; for more information see \href{https://arxiv.org/pdf/2007.14900.pdf}{Kontoyiannis et al. (2020)}.
//' @return returns a vector containing the averaged log-loss incurred in the sequential prediction at each time-step.
//' 
//' @seealso \code{\link{prediction}}, \code{\link{zero_one_loss}}
//' 
//' @export
//' @examples
//' # Compute the log-loss in the prediction of the last 10 elements 
//' # of a dataset. 
//' log_loss(pewee, 5, nchar(pewee) - 10)
//' 
//' # For custom beta (e.g. 0.7):
//' log_loss(pewee, 5, nchar(pewee) - 10, 0.7)
// [[Rcpp::export]]
NumericVector log_loss(CharacterVector input_data, IntegerVector depth, IntegerVector train_size, Nullable<NumericVector> beta = R_NilValue){
  int ttrain_size; 
  ttrain_size = train_size[0];
  D = depth[0];
  double b;
  
  std::string s = Rcpp::as<std::string>(input_data);
  if(beta.isNotNull()){
    NumericVector beta_(beta);
    b = beta_[0];
    set_global_parameters(s, D, 0, b);
  }
  else
    set_global_parameters(s, D, 0);
  
  return(compute_log_loss(xn, ttrain_size));
}
//==================================================================================

//' @title Bayesian Context Trees (BCT) algorithm 
//' @description Finds the maximum a posteriori probability (MAP) tree model. 
//'
//' @param input_data the sequence to be analysed. 
//' The sequence needs to be a "character" object. See the examples section on how to transform any dataset to a "character" object.
//' @param depth maximum memory length.
//' @param beta hyper-parameter of the model prior. 
//' Takes values between 0 and 1. If not initialised in the call function, the default value is \ifelse{html}{\out{1-2<sup>-m+1</sup>}}{\eqn{1 - 2^{-m+1}}}, 
//' where \ifelse{html}{\out{m}}{\eqn{m}} is the size of the alphabet; for more information see \href{https://arxiv.org/pdf/2007.14900.pdf}{Kontoyiannis et al. (2020)}.
//' 
//' @return returns a list object which includes:
//' \item{Contexts}{MAP model given as a list object containing the contexts of its leaves.}
//' \item{Results}{a dataframe with the following columns: prior probability, log(prior probability), posterior probability, log(posterior probability), number of leaves, maximum depth, BIC score, AIC score and maximum log-likelihood. }
//'
//'
//' @export
//' 
//' @seealso \code{\link{kBCT}}
//' 
//' @examples
//' # Finding the MAP model with maximum depth <= 10 
//' # for the SP500 dataset (with default value beta):
//' 
//' BCT(SP500, 10)  
//' 
//' # For custom beta (e.g. 0.7):
//' 
//' BCT(SP500, 10, 0.7)  
//' 
//' # The type of the input dataset is "character"
//' # If the dataset is contained within a vector:
//' 
//' q <- c(1, 0, 1, 0, 0, 1, 0, 0, 1, 1, 1, 0, 1, 0)
//' 
//' # Convert a vector to a "character" object:
//' s <- paste(q, collapse = "")
//'
//' BCT(s, 2)
//' 
//' # Reading a file using the readChar function 
//' 
//' # Example 1: The dataset is stored in a .txt file
//' 
//' # fileName <- '~/example_data.txt' # fileName stores the path to the dataset
//' 
//' # s<-readChar(fileName, file.info(fileName)$size)
//' 
//' # Make sure that s does not contain any "\n" at the end of the string
//' # To remove last entry:
//' # s<-gsub('.$', '', s)
//' 
//' # To remove any unwanted characters (e.g. "\n"):
//' # s<-gsub('\n', '', s)
//' 
//' # Example 2: The dataset is stored in a .csv file
//' 
//' # fileName <- '~/example_data.csv' # fileName stores the path to the dataset
//' 
//' # s<-readChar(fileName, file.info(fileName)$size)
//' 
//' # Depending on the running environment, 
//' # s might contain unwanted characters such as: "\n" or "\r\n".
//' # Remove any unwanted characters (e.g. "\r\n"):
//' # s<-gsub('\r\n', '', s)
//' 
//' # Always make sure that s does not contain any unwanted characters
// [[Rcpp::export]]
List BCT(CharacterVector input_data, IntegerVector depth, Nullable<NumericVector> beta = R_NilValue){
  int D;
  double b;
  D = depth[0];
  std::string s = Rcpp::as<std::string>(input_data);
  
  
  if(beta.isNotNull()){
    NumericVector beta_(beta);
    b = beta_[0];
    set_global_parameters(s, D, 0, b); // set k_max to 0 as the bct does not require k_max
  }
  else
    set_global_parameters(s, D, 0); // set k_max to 0 as the bct does not require k_max
  
  Tree_properties tp_vec = build_bct();
  List lst;
  
  NumericVector prior;
  NumericVector log_prior;
  NumericVector posterior;
  NumericVector log_posterior;
  NumericVector n_leaves;
  NumericVector max_depth;
  NumericVector bic;
  NumericVector aic;
  NumericVector mle;
  
 
  StringVector lst_context;
    
  for(unsigned int j = 0; j<tp_vec.context.size(); j++){
    lst_context.push_back(tp_vec.context[j]);
  }
    
  prior.push_back(tp_vec.prior);
  log_prior.push_back(tp_vec.log_prior * log(2));
  posterior.push_back(tp_vec.posterior);
  log_posterior.push_back(tp_vec.log_posterior * log(2));
  n_leaves.push_back(tp_vec.n_leaves);
  max_depth.push_back(tp_vec.max_depth);
  bic.push_back(tp_vec.bic);
  aic.push_back(tp_vec.aic);
  mle.push_back(tp_vec.mle);
  
  
  DataFrame df = DataFrame::create( _["prior"] = prior,
                                    _["log_prior"] = log_prior,
                                    _["posterior"] = posterior,
                                    _["log_posterior"] = log_posterior,
                                    _["number_leaves"] = n_leaves,
                                    _["max_depth"] = max_depth,
                                    _["BIC"] = bic,
                                    _["AIC"] = aic,
                                    _["max_log_lik"] = mle);

  List results = List::create(Named("Contexts") = lst_context, _["Results"] = df);
  return results;
} 

//==================================================================================
//' @title Calculating the 0-1 loss incurred in prediction
//' @description Compute the 0-1 loss, i.e., the proportion of incorrectly predicted values,
//' incurred in BCT prediction with memory length D. Given an initial context
//' \ifelse{html}{\out{(x<sub>-D+1</sub>, ..., x<sub>0</sub>)}}{\eqn{(x_{-D+1},...,x_0)}}  and training data \ifelse{html}{\out{(x<sub>1</sub>, ..., x<sub>n</sub>)}}{\eqn{(x_{1},...,x_n)}}, the 0-1 loss is computed in sequentially predicting
//'  the test data \ifelse{html}{\out{(x<sub>n+1</sub>, ..., x<sub>n+T</sub>)}}{\eqn{(x_{n+1},...,x_{n+T})}}. The function outputs the cummulative, normalized (per-sample) 0-1 loss, at each prediction step; for more information see \href{https://arxiv.org/pdf/2007.14900.pdf}{Kontoyiannis et al. (2020)}.
//' @param input_data the sequence to be analysed. 
//' The sequence needs to be a "character" object. See the examples section of kBCT/BCT functions on how to transform any dataset to a "character" object.
//' @param train_size number of samples used in the training set. The training set size should be at least equal to the depth.
//' @param depth maximum memory length.
//' @param beta hyper-parameter of the model prior. 
//' Takes values between 0 and 1. If not initialised in the call function, the default value is \ifelse{html}{\out{1-2<sup>-m+1</sup>}}{\eqn{1 - 2^{-m+1}}}, 
//' where \ifelse{html}{\out{m}}{\eqn{m}} is the size of the alphabet; for more information see \href{https://arxiv.org/pdf/2007.14900.pdf}{Kontoyiannis et al. (2020)}
//' 
//' 
//' @return returns a vector containing the averaged number of errors at each timestep. 
//' @export
//' @seealso \code{\link{log_loss}}, \code{\link{prediction}}
//' 
//' @examples
//' # Use the pewee dataset and look at the last 8 elements:
//'   substring(pewee, nchar(pewee)-7, nchar(pewee)) 
//' # [1] "10001001"
//' 
//' # Predict last 8 elements using the prediction function
//' pred <- prediction(pewee, 10, nchar(pewee)-8)[["Prediction"]] 
//' # Taking only the "Prediction" vector:
//' 
//' pred
//' # [1] "1" "0" "0" "1" "1" "0" "0" "1"
//' 
//' # To transform the result of the prediction function into a "character" object:
//' paste(pred, collapse = "")
//' # [1] "10011001"
//' 
//' # As observed, there is only 1 error (the sixth predicted element is 1 instead of a 0). 
//' # Thus, up to the 4th place, the averaged error is 0 
//' # and the sixth averaged error is expected to be 1/4. 
//' # Indeed, the zero_one_loss function yields the expected answer: 
//' 
//' zero_one_loss(pewee, 10, nchar(pewee)-8) 
//' # [1] 0.0000000 0.0000000 0.0000000 0.2500000 0.2000000 0.1666667 0.1428571 0.1250000
//' 
// [[Rcpp::export]]
NumericVector zero_one_loss(CharacterVector input_data, IntegerVector depth, IntegerVector train_size, Nullable<NumericVector> beta = R_NilValue){
  int ttrain_size;
  ttrain_size = train_size[0];
  D = depth[0];
  double b;
  
  std::string s = Rcpp::as<std::string>(input_data);
  if(beta.isNotNull()){
    NumericVector beta_(beta);
    b = beta_[0];
    set_global_parameters(s, D, 0, b);
  }
  else
    set_global_parameters(s, D, 0);
  
  List onl_pred = online_predict(ttrain_size); //uses the online_predict function to predict the most likely characters at each time step
  
  CharacterVector pred  = onl_pred["Prediction"];
  
  NumericVector out(xn.size()- ttrain_size); 
  int sum = 0;
  for(int i = 0; i< pred.size(); i++){
    string pr = std::string(pred[i]);
    string pr2(1, decoder[xn[ttrain_size+i]]);
    
    if(pr != pr2) // compares the prediction to the actual value; if different, an error occured
      sum ++;
    
    out[i] = sum/(i+1.0); // averaged number of errors
  }
  return(out);
}
//==================================================================================

//' @title Prediction
//' @description Computes the posterior predictive distribution at each time step, and predicts the next symbol as its most likely value. Given an initial context
//' \ifelse{html}{\out{(x<sub>-D+1</sub>, ..., x<sub>0</sub>)}}{\eqn{(x_{-D+1},...,x_0)}} and training data \ifelse{html}{\out{(x<sub>1</sub>, ..., x<sub>n</sub>)}}{\eqn{(x_{1},...,x_n)}}, the posterior predictive distribution is computed sequentially for the test data \ifelse{html}{\out{(x<sub>n+1</sub>, ..., x<sub>n+T</sub>)}}{\eqn{(x_{n+1},...,x_{n+T})}}. The function outputs the predicted distribution at each time step, along with the most likely symbol; for more information see \href{https://arxiv.org/pdf/2007.14900.pdf}{Kontoyiannis et al.(2020)}.
//' @param input_data the sequence to be analysed. 
//' The sequence needs to be a "character" object. See the examples section on how to transform any dataset to a "character" object.
//' @param depth maximum memory length.
//' @param train_size number of samples used for training.
//' @param beta hyper-parameter of the model prior. 
//' Takes values between 0 and 1. If not initialised in the call function, the default value is \ifelse{html}{\out{1-2<sup>-m+1</sup>}}{\eqn{1 - 2^{-m+1}}}, 
//' where \ifelse{html}{\out{m}}{\eqn{m}} is the size of the alphabet; for more information see: \href{https://arxiv.org/pdf/2007.14900.pdf}{Kontoyiannis et al. (2020)}.
//' 
//' @return returns a "list" containing the posterior predictive distribution at each time step. The last entry in the list, named "Prediction", contains the most likely character 
//' at each time step according to the posterior predictive distribution. 
//' 
//' @seealso \code{\link{log_loss}}, \code{\link{zero_one_loss}}
//' @export
//' @examples
//' # Predicting the 2 last characters of a dataset using a model with a maximum depth of 5
//' # The training size is the total number of characters within the dataset minus 2: nchar(pewee) - 2
//' 
//' q <- prediction(pewee, 5, nchar(pewee) - 2)
//' 
//' q
//' # [[1]]
//' # [1] 0.56300039 0.05899728 0.37800233
//' 
//' # [[2]]
//' # [1] 0.08150306 0.76293065 0.15556628
//' 
//' # $Prediction
//' # [1] "0" "1"
//' 
//' # To access the "Prediction" from result list q:
//' q[["Prediction"]]
//' 
//' # For custom beta (e.g. 0.8):
//' prediction(pewee, 5, nchar(pewee) - 10, 0.8)
// [[Rcpp::export]]
Rcpp::List prediction(CharacterVector input_data, IntegerVector depth, IntegerVector train_size, Nullable<NumericVector> beta = R_NilValue){
  int ttrain_size;
  ttrain_size = train_size[0];
  D = depth[0];
  double b;
  
  std::string s = Rcpp::as<std::string>(input_data);
  if(beta.isNotNull()){
    NumericVector beta_(beta);
    b = beta_[0];
    set_global_parameters(s, D, 0, b); // sets k_max = 0 
  }
  else
    set_global_parameters(s, D, 0); // sets k_max = 0 
  
  return(online_predict(ttrain_size));
}
//=====================================================================================================
//' @title Maximum Likelihood
//' @description Computes the logarithm of the likelihood of the observations, maximised over all models and parameters.
//' @param input_data the sequence to be analysed. 
//' The sequence needs to be a "character" object. See the examples section of the BCT/kBCT functions on how to transform any dataset to a "character" object.
//' @param depth maximum memory length.
//' @return returns the natural logarithm of the maximum likelihood.
//'
//' @seealso \code{\link{BCT}}, \code{\link{kBCT}}
//'
//' @export
//' 
//' @examples
//' # Computing the maximum likelihood of the gene_s dataset 
//' # with a maximum depth of 5:
//' ML(gene_s, 5)
// [[Rcpp::export]]
long double ML(CharacterVector input_data, IntegerVector depth){
  int D;
  D = depth[0];
  std::string s = Rcpp::as<std::string>(input_data);
  set_global_parameters(s, D, 0, 0);
  long double result =  mle_tree(); 
  return result;
} 
//==================================================================================

//' @title  Compute empirical frequencies of all contexts
//' @description Computes the count vectors of all contexts up to a certain length (D) for a given dataset. The first D characters are used to construct the initial context and the counting is performed on the remaining characters.
//' These counts are needed for intermediate computations in BCT and kBCT, 
//' and can also be viewed as maximum likelihood estimates of associated parameters; see \href{https://arxiv.org/pdf/2007.14900.pdf}{Kontoyiannis et al. (2020)}.
//' @param input_data the sequence to be analysed. 
//' The sequence needs to be a "character" object. See the examples section of the BCT/kBCT functions on how to transform any dataset to a "character" object.
//' @param depth maximum memory length.
//' @return a list containing the counts of all contexts of length \ifelse{html}{\out{&le;}}{\eqn{\le}} depth. 
//' If a context with a smaller length than the maximum depth is not contained in the output, its associated count vector is 0. 'Root' indicates the empty context.
//' @export
//' @seealso \code{\link{BCT}}, \code{\link{generate_data}}
//' @examples
//' # For the pewee dataset:
//' compute_counts(pewee, 3)
// [[Rcpp::export]]
List compute_counts(CharacterVector input_data, IntegerVector depth){
  int D;
  D = depth[0];
  std::string s = Rcpp::as<std::string>(input_data);
  
  set_global_parameters(s, D, 0);
  
  map<string, vector<int>> lst =  dictionary_counts();
  
  List out;
  for(std::map<string, vector<int>>::iterator it = lst.begin(); it != lst.end(); ++it) {
    out.push_back(it->second, it->first);
  }
  return(out);
} 

//==================================================================================

//' @title Context Tree Weighting (CTW) algorithm
//' @description Computes the prior predictive likelihood of the data.
//' @param input_data the sequence to be analysed. 
//' The sequence needs to be a "character" object. See the examples section of the BCT/kBCT functions on how to transform any dataset to a "character" object.
//' @param depth maximum memory length.
//' @param beta hyper-parameter of the model prior. 
//' Takes values between 0 and 1. If not initialised in the call function, the default value is \ifelse{html}{\out{1-2<sup>-m+1</sup>}}{\eqn{1 - 2^{-m+1}}}, 
//' where \ifelse{html}{\out{m}}{\eqn{m}} is the size of the alphabet; for more information see \href{https://arxiv.org/pdf/2007.14900.pdf}{Kontoyiannis et al. (2020)}.
//' 
//' @return returns the natural logarithm of the prior predictive likelihood of the data. 
//'
//' @seealso \code{\link{BCT}}, \code{\link{kBCT}}
//'
//' @export
//' 
//' @examples
//' # For the gene_s dataset with a maximum depth of 10 (with dafault value of beta):
//' CTW(gene_s, 10)
//' 
//' # For custom beta (e.g. 0.8):
//' CTW(gene_s, 10, 0.8)
// [[Rcpp::export]]
long double CTW(CharacterVector input_data, IntegerVector depth, Nullable<NumericVector> beta = R_NilValue){
  int D;
 // int k_max;
  double b;
  D = depth[0];
  std::string s = Rcpp::as<std::string>(input_data);
  
  
  if(beta.isNotNull()){
    NumericVector beta_(beta);
    b = beta_[0];
    set_global_parameters(s, D, 0, b);
  }
  else
    set_global_parameters(s, D, 0);
  
 long double prob = build_ctw_rcpp();
 Rcout<<"log-Prior predictive likelihood:"<<endl;
 return(prob*log(2)); // in the C files the logarithm is in base 2
}

// ========================================================================================================

//' @title Parameters of the MAP model
//' @description Returns the parameters of each leaf contained in the MAP model.
//' @param input_data the sequence to be analysed. 
//' The sequence needs to be a "character" object. See the examples section of BCT/kBCT functions on how to transform any dataset to a "character" object.
//' @param depth maximum memory length.
//' @param beta hyper-parameter of the model prior. 
//' Takes values between 0 and 1. If not initialised in the call function, the default value is \ifelse{html}{\out{1-2<sup>-m+1</sup>}}{\eqn{1 - 2^{-m+1}}}, 
//' where \ifelse{html}{\out{m}}{\eqn{m}} is the size of the alphabet; for more information see \href{https://arxiv.org/pdf/2007.14900.pdf}{Kontoyiannis et al. (2020)}.
//' @return list of parameters for each of the context within the MAP model.
//'
//' @export
//' 
//' @seealso \code{\link{BCT}}, \code{\link{kBCT}}, \code{\link{generate_data}}
//' 
//' @examples
//' 
//' # Use the gene_s dataset:
//' q <- BCT(gene_s, 10)
//' expected_contexts <- q[['Contexts']]
//' expected_contexts
//' # [1] "3"  "1"  "0"  "23" "20" "21" "22"
//' 
//' # For default beta:
//' v <- MAP_parameters(gene_s, 10)
//' 
//' # For custom beta (e.g. 0.8):
//' MAP_parameters(gene_s, 10, 0.8)
//' 
//' # generate a sequence of data using the generate_data function
//' s <- generate_data(v, 20000)
//' 
//' # Use BCT:
//' r <- BCT(s, 10)
//' 
//' # Check the resulting contexts:
//' r[['Contexts']]
//' # [1] "3"  "0"  "1"  "20" "22" "23" "21"
//' 
//' # The resulting contexts are as expected 
//' 
// [[Rcpp::export]]
List MAP_parameters(CharacterVector input_data, IntegerVector depth, Nullable<NumericVector> beta = R_NilValue){
  int D;
  // int k_max;
  double b;
  D = depth[0];
  std::string s = Rcpp::as<std::string>(input_data);
  
  if(beta.isNotNull()){
    NumericVector beta_(beta);
    b = beta_[0];
    set_global_parameters(s, D, 0, b);
  }
  else
    set_global_parameters(s, D, 0);
  
  List theta = map_param(); // stores the result in a list
  
  return(theta);
}
