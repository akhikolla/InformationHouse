#include <Rcpp.h>
#include "CreateRandomizedResponseBloomfilters.h"
using namespace Rcpp;

/**
 * Wrapper function for CreateBitFlippingBF
 */

// [[Rcpp::export]]
DataFrame CreateBitFlippingBF(DataFrame data, SEXP password, float f) {
  if (0>f ||1<f){
    Rcpp::Rcerr << "f has to be between 0 and 1 (including)! Otherwise the original CLK will be returned.";
    return data;
  }
  vector<string> RR(data.nrows());
  CharacterVector ID = data[0];
  if(TYPEOF(data[1])==STRSXP)
    RR=Rcpp::as<std::vector<string> > (data[1]);
  else
    Rcpp::Rcerr << "CLKs must be of type string/character";


  string password_ = as<string>(password);
  for (int i = 0; i < data.nrows(); ++i){
    RR[i]= CreateBitFlippingBFHelper(RR[i], password_, f);
  }
  return DataFrame::create(Named("ID") = ID ,
                           Named("BalancedBloomfilter") = RR,
                           Named("stringsAsFactors") = false);
}



/**
 * Creats Permanent Randomized Response.
 * he randomized response technique is used on each
 *bit position B[i] of a Bloom filter B:

 *B[i]= 1 with probability 1/2 f

 *B[i]= 0 with probability 1/2 f

 *B[i]= B[i] with probability 1 - f

 *with 0 <= f <= 1.
 */
string CreateBitFlippingBFHelper(string data, string password, float f) {
  string B=data;
  int k = data.size();
  //create vector of length of the clk with random values between 0 and 1
  //set seed accourding to the user defined password
  seed_seq seed(password.begin(), password.end());
  //create a random engine
  default_random_engine engine(seed); // or other engine as std::mt19937
  //vector for k random numbers
  vector<float> rands(k);
  uniform_real_distribution<> distr(0,1);
  //fill vector with k random numbers
  generate(begin(rands), end(rands), [&]() {return distr(engine);});

  //create "Permanent randomized response" for each posotion in clk
  for (int position =0; position < k; position++){
    if (rands[position] <= (f- (f/2))) {
      B[position] ='1';
    } else if (rands[position] > (f- (f/2)) && rands[position] <= f) {
      B[position]= '0';
    } else {
      B[position] = data[position];
    }
  }
  return B;
} //end CreateBitFlippingBFHelper
