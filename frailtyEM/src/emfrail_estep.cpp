#include <Rcpp.h>
#include <vector>
#include <cmath>
using namespace Rcpp;

void printvec(std::vector<int> v) {
  Rcout<<std::endl;
  Rcout<<"( ";
  for(std::vector<int>::iterator it = v.begin(); it != v.end(); it++)
    Rcout<<*it<<" ";
  Rcout<<")";
}


double logfactorial(const int &k) {
  double res = 0.0;
  if(k==0 || k==1) return res; else
    for(int i = 2; i<= k; i++)
      res += std::log(i);
  return res;
}

/*
 * this function takes a configuration vector and
 * then it does some krzy sht.
 * first, it calculates a product of the log factorial thing.
 */

double repeats2(std::vector<int> &vec, const int &nev) {
  double res = 0.0;
  int mult = 1;

  // if the first one is more than one, then add it in the numerator of this
  if(vec[0] > 1) res += logfactorial(vec[0]);

  // the end is either at the end of the vector or when we get a 0.
  for(unsigned int i = 1; i < vec.size(); i++) {

    if(vec[i] == 0) break;

    // if at position i there is a new one, then reset mult and correct for that;
    res += logfactorial(vec[i]); // for sort of numerator I guess
    //Rcout<<" * "<<vec[i];
    if(vec[i] != vec[i-1]) {
      res += logfactorial(mult); // this also makes a lot of useless calls when mult is 1
      //Rcout<<" * "<<mult<<"!";
      mult = 1;
    } else { // if they are equal, we just increase the multiplicity.
      mult++;
    }
  }

  if(mult > 1) res += logfactorial(mult); // there might be a leftover(or not)
  //Rcout<<" * "<<mult<<"!) = exp ";

  return res;
}


// this is the derivative of the Laplace exponent of the gamma distribution;
// same thing as alpha * psi basically and their derivatives

// the -1 signs here are all the ones except the - from the alpha. that is added in findsums
// then the -1 from the derivative - laplace connection is also later
double exponent_gamma(const double& alpha, const double& bbeta, const double& c, const double &c_lt, const int& nderiv) {
  if(nderiv == 0) {
    return alpha * (std::log(bbeta + c + c_lt) - std::log(bbeta + c_lt));
    //return alpha * (std::log(bbeta + c) - std::log(bbeta));

  } else
    return alpha * pow(bbeta + c + c_lt, -nderiv) * std::pow(std::exp(1.0), logfactorial(nderiv - 1)) *
      pow(-1.0, nderiv - 1);

    // return alpha * pow(bbeta + c, -nderiv) * std::pow(std::exp(1.0), logfactorial(nderiv - 1)) *
    //   pow(-1.0, nderiv - 1);
}

double exponent_stab(const double& alpha, const double& bbeta, const double& c, const double &c_lt, const int& nderiv) {
  if(nderiv == 0) {
    return alpha * (std::pow(c + c_lt, bbeta) - std::pow(c_lt, bbeta));
    //return alpha * std::pow(c, bbeta);
  } else
   //
    // problem here: lgamma can't handle negative values, beta is in (0,1).
    // for nderiv = 1, we want the product to be bbeta
    // nderiv =2 it will be bbeta * (bbeta - 1) = bbeta * (1 - bbeta) * (-1)  // already in the negative side
    // nderiv = 3 it will be bbeta * (bbeta-1) *(bbeta - 2) = bbeta * (1 - bbeta) * (2 - bbeta) * (-1^2) ///// this is essentially
    // bbeta * (1-bbeta) * (2 - bbeta) * ()
    // return alpha * pow(c, bbeta - nderiv) * bbeta * std::exp(lgamma(nderiv - bbeta) - lgamma(1-bbeta)) * pow(-1.0, nderiv + 1) ;
    return alpha * pow(c + c_lt, bbeta - nderiv) * bbeta * std::exp(lgamma(nderiv - bbeta) - lgamma(1-bbeta)) * pow(-1.0, nderiv + 1);
  //*
  //    pow(-1.0, nderiv - 1);

}

/*
 * pvfm is a constant, it is not estimated.
 */
double exponent_pvf(const double& alpha, const double& bbeta, const double &pvfm,
                    const double& c, const double &c_lt, const int& nderiv) {

  double sign = 1.0;
  if(pvfm < 0) sign = -1.0;

  if(nderiv == 0) {
    return alpha * (std::pow(bbeta / (bbeta + c_lt), pvfm) - std::pow(bbeta / (bbeta + c + c_lt), pvfm)) * sign;
    //return alpha * (1 - std::pow( bbeta / (bbeta + c), pvfm)) * sign;
  } else  {
    // Rcout<<pow(bbeta, pvfm)<<" * "<<std::pow(bbeta + c, -pvfm - nderiv)<<" * "<<std::exp(lgamma(pvfm + nderiv) - lgamma(pvfm));
    // Rcout<<" = "<<alpha * std::pow(bbeta, pvfm) * std::pow(bbeta + c, -pvfm - nderiv) * std::exp(lgamma(pvfm + nderiv) - lgamma(pvfm))<<std::endl;
    // Rcout<<"sign = "<<pow(-1.0, nderiv) * sign;
    //return alpha * std::pow(bbeta, pvfm) * std::pow(bbeta + c, -pvfm - nderiv) * std::exp(lgamma(pvfm + nderiv) - lgamma(pvfm)) *  pow(-1.0, nderiv + 1);// * sign;
    return alpha * std::pow(bbeta, pvfm) * std::pow(bbeta + c + c_lt, -pvfm - nderiv) *
       std::exp(lgamma(pvfm + nderiv) - lgamma(pvfm)) *  pow(-1.0, nderiv + 1);
  }

}





/*
 * this returns all combinations of numbers that sum up to n.
 * rest should start at n - keeps track of how much there is left to "spend"
 * last should start at 1 - keeps track of what new number to try
 * pos should start at 0 and it determines the position at which I should search a solution for the recursive problem.
 * combs is the combination vector of length n.
 * ncalls counts how many calls to the function are actually happening.
 */
void findsums(int rest, int last, int pos, std::vector<int> combs, //int &ncalls,
              const double& alpha, const double& bbeta, const double& c, const double &c_lt, double& res,
              const double& pvfm, const int& dist) {
  //printvec(combs);
  if(rest == 0) {

    // printvec(combs);

    //calculate coefficient

    //Rcout<<std::endl;

    // coefficient



    // Rcout<<std::endl<<"unadjusted coefficient: "<<combs.size()<<"! / ";
    // Rcout<<"(repeats : ";
    // double adjust = repeats2(combs, combs.size());
    // Rcout<<")";
    //
    // Rcout<<"factorial numerator is "<<std::exp(logfactorial(combs.size()))<<std::endl;
    // Rcout<<"log of denominator is "<<repeats2(combs, combs.size())<<std::endl;
    // Rcout<<"factorial denominator is "<<std::exp(repeats2(combs, combs.size()));


    double element = std::exp(logfactorial(combs.size()) - repeats2(combs, combs.size()));

    //Rcout<<std::endl<<"res = "<<element;

    // the value itself (the product)
    for(std::vector<int>::iterator it = combs.begin(); it!=combs.end(); it++) {
      if(*it == 0) break; else {
        //Rcout<<" * "<<-1.0 * exponent_gamma(alpha, bbeta, c, *it);
        // this is the -1 that comes from - alpha all the time basically
        if(dist == 0) {
          element = element * -1.0 * exponent_gamma(alpha, bbeta, c, c_lt, *it);
        } else {
          if(dist == 1) {
            element = element * -1.0 * exponent_stab(alpha, bbeta, c, c_lt, *it);
          } else
            element = element * -1.0 * exponent_pvf(alpha, bbeta, pvfm, c, c_lt, *it);
        }

      }
    }
    // Rcout<<" = "<<element;
    // Rcout<<"repeats is..."<<repeats(combs, combs.size());
    // Rcout<<"coef is..." <<std::exp(coef);
    //
    // Rcout<<std::endl<<"contributing with..."<<element;

    res += element;


  } else if(last <= rest)
    for(int i = last; i <= rest; i++) {


      //ncalls++ ;
      combs[pos] = i;
      findsums(rest - i, i, pos + 1, combs, alpha, bbeta, c, c_lt, res, pvfm, dist);
    }

}


/*
 * This is just to make tings more clear in terms of arguments and what goes where.
 */

double wrap_integral(int n, const double& alpha, const double &bbeta, const double& c,
                     const double &c_lt,
                     const double& pvfm, const int& dist) {
  std::vector<int> v(n);

  double res = 0.0;

  findsums(n, 1, 0, v, alpha, bbeta, c, c_lt, res, pvfm, dist);

  //Rcout<<std::endl<<"wrap_integral("<<n<<", alpha="<<alpha<<", c="<<c<<") = "<<res<<std::endl;
  return res;
}


/*
 * Input: c vector with cumulative hazard for each indviidual / cluster
 * delta vector with number of events for each individual / cluster
 * alpha bbeta parameters of the frailty distribution
 * toDo: make a parameter to decide which distribution to use
 */



//' Perform the E step calculations
//'
//' This is an inner wrapper for the C++ functions which perform the E step and is not intended to be used directly.
//' This function does not check the input.
//' For a data set with \code{K} clusters,
//' @param c Vector of length \code{K} of cumulative hazards, i.e. total accumulated hazards within a cluster
//' @param c_lt Vector of length \code{K} of cumulative hazard from 0 to the left truncation time
//' @param delta Vector of integers of length \code{K} of the number of events for each cluster
//' @param alpha,bbeta Parameters of the frailty distribution
//' @param pvfm Parameter for the PVF distribution, only matters in that case
//' @param dist One of 0 (for gamma), 1 (for stable) or 2 (for PVF)
//'
//' @return A \code{K x 3} matrix where the first column and the second column are the numerators
//' and the denominators of the frailty fraction (without the Laplace transform) and the
//' last column is the log(denominator) + log-Laplace transform, i.e. the log-likelihood contribution
//' @importFrom Rcpp evalCpp
//' @keywords internal
// [[Rcpp::export]]
NumericMatrix Estep(NumericVector c, NumericVector c_lt, IntegerVector delta, double alpha, double bbeta,
                    const double& pvfm, const int& dist) {

  /*
   * column 1: numerators
   * column 2: denominators
   * column 3: Laplace transform
   * z is a matrix by row
   *
   * input: dist==0 is gamma
   * dist == 1 is stable
   * dist == 2 is PVF
   */

  NumericVector z(c.size() * 3);
  NumericMatrix Z(c.size(), 3);

  for(int i = 0; i < delta.size(); i++) {
    //Rcout<<i<<" ";
    int n = delta[i];
    //Rcout<<"n events = "<<n;
    if(n==0) {

      //Rcout<<std::endl<<"numerator: "<<std::endl;
      Z(i, 0) = pow(-1.0, 1) * wrap_integral(1, alpha, bbeta, c[i], c_lt[i], pvfm, dist);
      Z(i, 1) = 1;
       // log laplace transform (- alpha phi(c))
      //
      // z[3 * i] = pow(-1.0, 1) * wrap_integral(1, alpha, bbeta, c[i]);
      // z[3 * i + 1] = 1;
      // z[3* i + 2] = -1.0 * exponent_gamma(alpha, bbeta, c[i], 0); // log laplace transform
        // - alpha phi(c) basically ;; minus not included in exponent_gamma.
    } else {
      //Rcout<<std::endl<<"numerator: "<<std::endl;
      // Rcout<<" at i = "<<i<<std::endl;

      // here is the -1 that comes from the sum structure.
      // i.e. if there is 1 then we need a -1.
      // as a consequence, wrap_integral should generally be negative for odd first argument.

      Z(i, 0) = pow(-1.0, n + 1) * wrap_integral(n + 1, alpha, bbeta, c[i], c_lt[i], pvfm, dist);
      //Rcout<<std::endl<<"denominator: "<<std::endl;


      Z(i, 1) = pow(-1.0, n) * wrap_integral(n , alpha, bbeta, c[i], c_lt[i], pvfm, dist);


      // z[3 * i] =  pow(-1.0, n + 1) * wrap_integral(n + 1, alpha, bbeta, c[i]);
      // z[3 * i + 1] = pow(-1.0, n) * wrap_integral(n , alpha, bbeta, c[i]);
      // z[3 * i + 2] = -1.0 * exponent_gamma(alpha, bbeta, c[i], 0);
    }

      if(dist == 0) {
        Z(i,2) = -1.0 * exponent_gamma(alpha, bbeta, c[i], c_lt[i], 0) + std::log(Z(i,1));
      } else {
        if(dist == 1) {
          Z(i,2) = -1.0 * exponent_stab(alpha, bbeta, c[i], c_lt[i], 0) + std::log(Z(i,1));
        } else
          Z(i,2) =  -1.0 * exponent_pvf(alpha, bbeta, pvfm, c[i], c_lt[i], 0) + std::log(Z(i,1));
      }
    //Z(i, 2) is the contribution to the log-likelihood now

  }
  return Z;

}

/*** R
# Estep(c = rexp(1), c_lt = 0, delta = 2, alpha = 1, bbeta = 1, pvfm = -0.5, dist = 2)

*/

