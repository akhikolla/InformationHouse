//============================================================================
// Name        : CreateBalancedBF.cpp
// Author      : Rukasz
// Version     :
// Copyright   : Your copyright notice
// Description : Ansi-style
//============================================================================

// Enable C++11 via this plugin
// [[Rcpp::plugins("cpp11")]]

#include <Rcpp.h>
#include "CreateBalancedBF.h"
#include "checkVectors.h"
using namespace Rcpp;

/**
 * Creates a Balanced Bloom filters of length 2 * l, which
 have the Hamming weight l.
 */
string CreateBalancedBloomfilterHelper(string clk, string password) {

  string negClk = clk;
  string bb;
  for (size_t i = 0; i < clk.size(); ++i)
  {
    if (clk[i]=='0')
      negClk[i]='1';
    else if (clk[i]=='1')
      negClk[i]='0';
    else{
      Rcpp::Rcerr << "Bloomfilter has to consist of zeros and ones only." << endl;
      return 0;
    }
  }
   bb =clk+negClk;
    //permutate binary vector
  //set seed accourding to the user defined password
  seed_seq seed(password.begin(), password.end());
  //create a random engine
  default_random_engine engine(seed); // or other engine as std::mt19937
  shuffle(bb.begin(), bb.end(), default_random_engine(seed));

  return bb;

} //end CreateBalancedBloomfilter


/**
 * Wrapper function for CreateBalancedBF from CLK or BF
 */
// [[Rcpp::export]]
DataFrame CreateBalancedBF(CharacterVector ID, CharacterVector data, SEXP password){
  if (data.size()!=ID.size()){
    Rcpp::Rcerr << " ID-Vector and Input-Data must have the same size. " << endl;
    return NULL;
  }
  string password_ =  as<string>(password);
  vector<string> BB = Rcpp::as<std::vector<string> > (data);

  for (int i =0; i <data.size();i++){
    //create Balanced Bloomfilter
    BB[i] = CreateBalancedBloomfilterHelper(BB[i], password_);
  }
  //shuffle(begin(ID)+ID.length(), end(ID), default_random_engine(seed));
  return DataFrame::create(Named("ID") = ID ,
                           Named("BB") = BB,
                           Named("stringsAsFactors") = false);

} // end CreateBalancedBF


/**
 * Wrapper function for CreateDoubleBalancedBF from CLK or BF
 */
// [[Rcpp::export]]
DataFrame CreateDoubleBalancedBF(CharacterVector ID, CharacterVector data, SEXP password){
  if (data.size()!=ID.size()){
    Rcpp::Rcerr << " ID-Vector and Input-Data must have the same size. " << endl;
    return NULL;
  }
  string password_ =  as<string>(password);
  vector<string> BB = Rcpp::as<std::vector<string> > (data);
  string neg;
  vector<string> BBneg(data.size());

  for (int i =0; i <data.size();i++){
    //create Balanced Bloomfilter
    BB[i] = CreateBalancedBloomfilterHelper(BB[i], password_ );
    //copy BB
   // BBneg[i] = BB[i];
    //negate BB
    string neg=BB[i];
      for (size_t j = 0; j< neg.size(); ++j)
      {
        if (BB[i][j]=='0')
          neg[j]='1';
        else if (BB[i][j]=='1')
          neg[j]='0';
        else{
          Rcpp::Rcerr << "Bloomfilter has to consist of zeros and ones only." << endl;
          return 0;
        }
      }
      BBneg[i]=neg;

      //permute BBneg[i]
      seed_seq seed(password_.begin(), password_.end());
      //create a random engine
      default_random_engine engine(seed); // or other engine as std::mt19937
      // permute positions of neg
      shuffle(BBneg[i].begin(), BBneg[i].end(), default_random_engine(seed));
      //BBneg[i]=CreateBalancedBloomfilterHelper(neg, passwordWithSalt);

    }


  //shuffle negated
  seed_seq seed(password_.begin(), password_.end());
  //shuffle(begin(BBneg), end(BBneg), default_random_engine(seed));
  //concatenate
  BB.insert(BB.end(),BBneg.begin(), BBneg.end());
  shuffle(begin(BB), end(BB), default_random_engine(seed));
  //adapt ID
  CharacterVector ID_(ID.length()*2);
  for (int i = 0; i< ID.length(); i++){
    ID_[i]= ID[i];
    ID_[i+ID.length()] = "f"+ID[i];
  }
  shuffle(begin(ID_), end(ID_), default_random_engine(seed));
  return DataFrame::create(Named("ID") = ID_ ,
                           Named("DBB") = BB,
                           Named("stringsAsFactors") = false);


} // CreateDoubleBalancedBF
