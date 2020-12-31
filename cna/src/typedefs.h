
#include <Rcpp.h>
using namespace Rcpp;

typedef ListOf<IntegerVector> intList;
typedef ListOf<intList> recIntList; // recursive intList = list of intList
typedef ListOf<NumericMatrix> numMatList;

typedef ListOf<CharacterVector> charList;
typedef ListOf<charList> recCharList;
