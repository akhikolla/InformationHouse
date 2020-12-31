#include "prepareData.h"

/**
* gets data from different types
*/
vector<string> prepareData(SEXP data_, string dataName, bool em) {
  vector < string > res;
  vector<float> varsV1Float;
  vector<int> varsV1Int;
  bool isAsciiString = false;

  // Transform dataA_ to vector<string>
  if (TYPEOF(data_) == STRSXP) { // case data is of type characters
    res = Rcpp::as < std::vector<string> > (data_);
    for (unsigned j = 0; j < res.size() && !isAsciiString; j++) {
      isAsciiString = isAscii(res[j], dataName); //in case of non ascii characters errors may occure, so a flag is set
    }
  }
  //Input (floats) is tranformed to strings
  else if (TYPEOF(data_) == REALSXP) {
    Rcpp::Rcerr
    << "Warning: vars1 contains floats. Data will be transformed to characters."
    << endl;
    varsV1Float = Rcpp::as < std::vector<float> > (data_);
    if (varsV1Float.size() > 0) {
      for (unsigned j = 0; j < varsV1Float.size(); j++) {
        res[j] = to_string(varsV1Float[j]); //Copy the vector to the string
        if (!isAsciiString) { //in case of non ascii characters errors may occure, so a flag is set
          isAsciiString = isAscii(res[j], dataName);
        }
      }
    }
  }
  //Input (int) is tranformed to strings
  else if (TYPEOF(data_) == INTSXP) {
    Rcpp::Rcerr
    << "Warning: data contains integers or factors. Make sure the input data are not factors when you want to use characters. Data will be transformed to characters."
    << endl;
    varsV1Int = Rcpp::as < std::vector<int> > (data_);
    if (varsV1Int.size() > 0) {
      for (unsigned j = 0; j < varsV1Int.size(); j++) {
        res[j] = to_string(varsV1Int[j]); //Copy the vector to the string
        if (!isAsciiString) { //in case of non ascii characters errors may occure, so a flag is set
          isAsciiString = isAscii(res[j], dataName);
        }
      }
    }
  } else if (!em) {
    Rcpp::Rcerr << "Error: data_ must be of type characters, int or float."
                << endl;
    return res;
  }
  return res;
}


bool isAscii(string str, string dataName) {
  bool res = (str.find_first_not_of(
    " \"!#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[//]^_`abcdefghijklmnopqrstuvwxyz{|}~")
                != std::string::npos);
  if (res) {
    Rcpp::Rcerr << dataName
                << ": Some strings contain illegal characters which might cause errors, e.g. '"
                << str << "'.\nThe preprossessing routine StandardizeString("<<dataName<<") can be used.\n";
  }
  return res;
}
