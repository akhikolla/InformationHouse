#include <Rcpp.h>
using namespace Rcpp;

// Function calculates total mechanical power
// @param bm body mass
// @param ws wing span
// @param wa wing area
// @param vt true speed
// @return totalMechPower total mechanical power
// @export

// [[Rcpp::export(.total_Mech_Pow_cpp)]]
double total_Mech_Pow_cpp(double bm, double ws, double wa,
                                 double vt, double g, double airDensity,
                                 double ipf, double bdc, double ppc){
  double bodyFrontArea = (0.00813 * pow(bm, 0.666));

  double inducedPower = 2 * ipf * pow((bm * g), 2) / (vt * M_PI * pow(ws, 2) * airDensity);

  double parasitePower =  (airDensity * pow(vt, 3) * bodyFrontArea * bdc) / 2;

  double profilePower = (ppc / (pow(ws, 2) / wa)) *
    (1.05 * pow(ipf, 0.75) * pow(bm, 1.5) * pow(g, 1.5) * pow(bodyFrontArea, 0.25) * pow(bdc, 0.25)) /
      (pow(airDensity, 0.5) * pow(ws, 1.5));

  double totalMechPower = profilePower + inducedPower + parasitePower;

  return totalMechPower;
}
