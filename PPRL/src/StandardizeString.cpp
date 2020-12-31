#include <Rcpp.h>
#include <vector>
#include "Standardisation.h"
using namespace Rcpp;
using namespace std;


/**
 * Preprocessing routine
 */
//[[Rcpp::export]]
CharacterVector StandardizeString(CharacterVector strings) {
  CharacterVector res = NULL;
  if (TYPEOF(strings) == STRSXP) { // case input is of type character
    vector < string > vars = Rcpp::as < std::vector<string> > (strings);
    for (unsigned i = 0;  i< vars.size(); i++) {
    replaceNonAscii(vars[i]);
    }
    res = vars;
    return res;
  }
  else {
    Rcpp::Rcerr << "Please enter a character vector!" << endl;
  }

  return res;
}


