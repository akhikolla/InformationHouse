// Enable C++11 via this plugin
// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::plugins("unwindProtect")]]

#include <Rcpp.h>
#include "FunctionWrapper.h"

using namespace Rcpp;
string defaultFunction = "unorderedPairing";
/**
* Wrapper function for Pairing functions
* Unordered Pairing Function applied on a dataframe with two columns:
*
* elegantPairing(x,y)=x*y+floor( ( (abs(x-y)-1)**2)/4)
*
* elegantPairing(x,y) = elegantPairing(y,x)
*
* When x and y are non−negative integers, elegantPairing(x,y) outputs a single
* non−negative integer that is uniquely associated with that unorderd pair.
*/
// [[Rcpp::export]]
DataFrame ElegantPairingVec(CharacterVector ID, DataFrame data) {
  vector<string> dataString(data.nrows());
  vector<int>  dataV(data.nrows());
  vector<int> inputVector(data.size());
  vector<vector<int> > dataM ;
  IntegerVector Pairing(data.nrows());

  //Sizes of Vectors are checked and adapted
  if (data.size()!=2){
    Rcpp::Rcerr << "The input data.frame must have two columns. Please check."<< endl;
    return NULL;
  }

  //Matrix is filled with input
  for (int i = 0 ; i< data.size(); i++){
    //if an input is a charactor vector
    if(TYPEOF(data[i])==STRSXP) {
      Rcpp::Rcerr << "Vector["<<i<<"] of the input data.frame is a charactor vector. It will be converted to an integer vector. Please check."<< endl;
      dataString=Rcpp::as<std::vector<string> > (data[i]);
      for (int j =0; j < data.nrows(); j++){
        dataV[j] = atoi(dataString[j].c_str());
      }
    }
     if(TYPEOF(data[i])!=INTSXP && TYPEOF(data[i])!=STRSXP) {
     Rcpp::Rcerr << "The input data.frame must be of type integer or character. Please check."<< endl;
       Rcpp::Rcerr << "Type: " << TYPEOF(data[i]) << endl;
       return NULL;
     }
     if(TYPEOF(data[i])==INTSXP) //input are integer vectors
    {

      dataV=Rcpp::as<std::vector<int> > (data[i]);
    }
    dataM.push_back(dataV);
  }

  for (int i = 0 ; i<data.nrows(); i++){
    for (int j = 0 ; j<data.size(); j++){
      inputVector[j]=dataM[j][i];
    }
    Pairing[i]=unorderedPairing(inputVector[0], inputVector[1]);
  }
  return DataFrame::create(Named("ID") = ID ,
                           Named("Pairing") = Pairing,
                           Named("stringsAsFactors") = false);

}

// [[Rcpp::export]]
SEXP ElegantPairingInt(int int1 ,int int2) {
  SEXP result = 0;
  result = Rcpp::wrap(unorderedPairing(int1, int2));
  return result;
}
