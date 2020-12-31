#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export(.basal_metabolic_pow)]]
double basal_metabolic_pow(double airframeMass,
                                  double muscleMass,
                                  int taxon,
                                  double alphaPasserines,
                                  double alphaNonPasserines,
                                  double deltaPasserines,
                                  double deltaNonPasserines) {

  //double body_mass = airframe_mass + muscle_mass;
  // double bmp = alpha_passerines * pow(body_mass,delta_passerines);

  //return x * 2;
  if (taxon == 1) {
    return alphaPasserines * pow((airframeMass + muscleMass),deltaPasserines);
  }else if (taxon == 2) {
    return alphaNonPasserines * pow((airframeMass + muscleMass),deltaNonPasserines);
  }else {
    stop("taxon must be of level 1/2");
  }
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

// /*** R
// timesTwo(42)
// */
