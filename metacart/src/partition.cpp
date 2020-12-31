#include <Rcpp.h>
using namespace Rcpp;


bool contain_(CharacterVector x1, CharacterVector x2){
  // Test if x1 contains the first element of x2 
  return std::find(x1.begin(), x1.end(), x2[0]) != x1.end();
}


//' Partition the test set based on a trained tree
//' 
//' @param x1 the tree component of the REmrt object
//' @param x2 the moderators in the test set
//' @param x3 indicates whether a moderator is numeric or not
//' @param x4 the index vector of the spliting moderators
//' @param x5 the list of split points
//' @param x6 the moderators in the training set
//' @keywords internal
// [[Rcpp::export(".partition")]]
IntegerMatrix partition(DataFrame x1, DataFrame x2, 
                    LogicalVector x3, IntegerVector x4,
                    List x5, DataFrame x6) {
  IntegerVector pleaf = x1["pleaf"];
  CharacterVector mod = x1["mod"];
  CharacterVector NewModNames = x2.names();
  IntegerVector pnode;
  IntegerMatrix res(x2.nrows(), x1.nrows());
  int i; 
  // int j;
  for (i = 0; i < x2.nrows(); i++) {// put all observations in the root node
    pnode.push_back(1);
    }
  res(_, 0) = pnode;
  int j;
  for (j = 1; j < x1.nrows(); j++) {
    if (x3[j] == true) {
      NumericVector sv = x2[x4[j]-1]; // spliting moderators in new data
      NumericVector tempSP = x5[j-1]; // split points
      for (i = 0; i < x2.nrows(); i++) {
        if (pnode[i] == pleaf[j]) {
          pnode[i] = 2*j;
          if (sv[i] > tempSP[0]) {
            pnode[i] = pnode[i]+1;
          }
        }

      }
    } else {
      CharacterVector sv = x2[x4[j]-1];
      CharacterVector tempSP = x5[j-1];
      CharacterVector tempModOld = x6[x4[j]-1]; //in training set
      for (i = 0; i < x2.nrows(); i++) {
        if (pnode[i] == pleaf[j]) {
          pnode[i] = 2*j;
          CharacterVector tempMod;
          tempMod.push_back(sv[i]);
          if (!contain_(tempModOld, tempMod)) {
            pnode[i] = NA_REAL; // if a new category is observed 
            // in the test set, assign NA 
          } else {
            if (!contain_(tempSP, tempMod)) {
              pnode[i] = pnode[i]+1;
            }
          }

        }

    }

    }
    res(_, j) = pnode;
  }
  return res;
    
}
  



