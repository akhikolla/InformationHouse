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

// [[Rcpp::export(.mechanical_power)]]
double power_curve(double bm, double ws, double wa, double tas, double g,
                   double airDensity, double ipf, double bdc, double ppc) {
  // reduce speed
  double speed = tas * 0.5;
  double bodyFrontArea = (0.00813 * pow(bm, 0.666));
  double diskArea = (M_PI * pow(ws, 2))/4;
  double vmr = (pow(ipf, 0.25) * pow(bm, 0.5) * pow(g, 0.5))/
    ( pow(airDensity, 0.5) * pow(bodyFrontArea * bdc, 0.25) * pow(diskArea, 0.25));
  double profilePower = (ppc / (pow(ws, 2) / wa)) *
    (1.05 * pow(ipf, 0.75) * pow(bm, 1.5) * pow(g, 1.5) * pow(bodyFrontArea, 0.25) * pow(bdc, 0.25)) /
      (pow(airDensity, 0.5) * pow(ws, 1.5));

  NumericVector totalMechPower;
  // while loop till vmr
  while (speed < vmr) {
    // calculate mechanical power
    double inducedPower = 2 * ipf * pow((bm * g), 2) / (speed * M_PI * pow(ws, 2) * airDensity);
    double parasitePower =  (airDensity * pow(speed, 3) * bodyFrontArea * bdc) / 2;
    totalMechPower.push_back(profilePower + inducedPower + parasitePower);
    speed = speed + 0.1;
  }

  double mechPower = Rcpp::min(totalMechPower);

  return mechPower;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

///*** R
//timesTwo(42)
// */
