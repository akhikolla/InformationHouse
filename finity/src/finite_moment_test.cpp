// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(BH)]]
// [[Rcpp::plugins("cpp11")]]

#define BOOST_DISABLE_ASSERTS

/* Following two lines to suppress deprecation warnings about 
 *  integer_log2.hpp.
 */
#define BOOST_PENDING_INTEGER_LOG2_HPP
#include <boost/integer/integer_log2.hpp>

#include <RcppArmadillo.h>
#include <boost/random.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/chi_squared.hpp>
#include <math.h>
#include <cmath>
using namespace Rcpp;

/*
// Fail structure for catching invalid values in overloaded function std_norm_moment
struct fail {};
*/

// Overloaded std_norm_moment function (to handle int and double inputs)
// Template
template <typename T>
double std_norm_moment(const T& t) {
    return -1;
}

// int instance
template <>
double std_norm_moment(const int& k) {
    /* Returns moment of order k of standard normal
     * Arguments: 
     *      k (int): moment order
     * Returns:
     *      (double): moment value
     */
    
    double moment = 1.0;
    for (int i=k-1; i>0; i-=2) {
        moment *= i;
    }
    if (k % 2 != 0) {   // if moment is odd
        moment *= pow(2.0 / M_PI, 0.5);
    }
    return(moment);
}

// double instance
template <>
double std_norm_moment(const double& k) {
    /* Returns moment of order k of standard normal
     * Arguments: 
     *      k (double): moment order
     * Returns:
     *      (double): moment value
     */
    
    double moment = pow(2, k / 2.0) * tgamma((k + 1.0) / 2.0) / pow(M_PI, 0.5);
    return(moment);
}

// Further function definitions.

double overloaded_std_norm_moment(double x) {
    /* Returns moment of order k of standard normal. Handles switch between int and double input.
     * Arguments: 
     *      k (any struct, but should be int or double): moment order
     * Returns:
     *      (double): moment value
     */
    if ((int)x == x) {
        return std_norm_moment<double>((int)x);
    } else {
        return std_norm_moment<double>(x);
    }
}

arma::uvec build_zeta_sorted(arma::vec xim, arma::vec u_series) {
    /* Computes vectors zeta (Trapani 2016 eq 6). Returns vector of sums over these vectors. 
     * This is more efficient than returning many long bool vectors.
     * Arguments: 
     *      xim (armadillo numeric vector): Modified artificial sample xi variates 
     *                                          exp(mu/2)*xi (Trapani 2016, Step 2)
     *      u_series (armadillo numeric vector): Sampling points
     * Returns:
     *      (armadillo unsigned int vector): Sums of Trapani's zeta indicator vectors
     */
    
    // Make sure the input vectors are sorted
    xim = arma::sort(xim);
    u_series = arma::sort(u_series);
    
    // Prepare computation
    unsigned int r = u_series.size();
    arma::uvec zeta(r);
    //zeta.fill(0);
    
    // Compute zeta
    unsigned int idx_xi = 0;
    for (unsigned int idx_u = 0; idx_u < r; idx_u++) {
        while ((idx_xi < r) && (xim(idx_xi) < u_series(idx_u))) {
            idx_xi++;
        }
        zeta(idx_u) = idx_xi;
    }
    return(zeta);
}

arma::uvec build_zeta_unsorted(arma::vec xim, arma::vec u_series) {     // This one is horribly inefficient. Do not use this.
    /* Computes vectors zeta (Trapani 2016 eq 6). Returns vector of sums over these vectors. 
     * This is more efficient than returning many long bool vectors.
     * 
     * This function is horrible inefficient and only here for comparison and 
     * verification purposes.
     * 
     * Arguments: 
     *      xim (armadillo numeric vector): Modified artificial sample xi variates 
     *                                          exp(mu/2)*xi (Trapani 2016, Step 2)
     *      u_series (armadillo numeric vector): Sampling points
     * Returns:
     *      (armadillo unsigned int vector): Sums of Trapani's zeta indicator vectors
     */
    // Compute zeta
    unsigned int r = u_series.size();
    arma::uvec zeta(r);
    for (unsigned int idx_u = 0; idx_u < r; idx_u++) {
        arma::uvec u_indicator = arma::find(xim < u_series(idx_u));
        //Rcpp::Rcout << "Indicator:" << u_indicator.t() << std::endl;
        zeta(idx_u) = u_indicator.size();
        //Rcpp::Rcout << "Sum:" << zeta(idx_u) << std::endl;
    }
    return(zeta);
}

//' Chi^2(1) Percentile
//'
//' @description Returns the Chi^2(1) percentile for the test statistic.
//' @param value Chi^2(1) value (type: double).
//' @return Chi^2(1) percentile (type: double).
//' @examples
//' get_chisq1_percentile(20.0)
//' @export
// [[Rcpp::export]]
double get_chisq1_percentile(double value) {
    /* Returns Chi^2(1) percentile for test.
     * Arguments: 
     *      value (double): Chi^2(1) value
     * Returns:
     *      (double): Chi^2(1) percentile
     */
    if (boost::math::isnan(value)) {
        return(NAN);
    } else {
        boost::math::chi_squared_distribution<> chi2(1);
        return(cdf(boost::math::complement(chi2, value)));
    }
}

//' Absolute Moment of Order k
//' 
//' @description Computes the absolute moment of order k of a sample of observations.
//' @param obs Observations (type: armadillo numeric vector).
//' @param k Moment order (type: double)
//' @return Moment value (type: double)
//' @examples
//' rvs <- stabledist::rstable(100000, 1.9, 0.5, 1, 0, pm = 0)
//' absolute_moment <- compute_absolute_moment(rvs, 2)
//' @export
// [[Rcpp::export]]
double compute_absolute_moment(arma::vec obs, double k) {
    /* Computes absolute moment of order k.
     * Arguments: 
     *      obs (armadillo numeric vector): Observations
     *      k (double): moment order
     * Returns:
     *      (double): moment value
     */

    // Remove nan's
    obs = obs.elem(arma::find_finite(obs));

    // Compute mu
    int N = obs.size();
    double mu = 0;
    for (unsigned int i = 0; i < N; i++) {
        mu += pow(fabs(obs(i)), k);
    }
    mu /= N;
    return(mu);
}

//' Finite Moment Test 
//' 
//' @description Computes Trapani's (2016) finite moment test for moment of order k of the distribution of a given the sample of observations obs. Knowledge of the identity of the distribution is not required. The null hypothesis is that the moment is infinite; the alternative is that it is finite. The function takes parameters of the test as optional arguments; some insights into the impact the choice of parameter values has are given in Trapani (2016).
//' @param obs Observations (type: armadillo numeric vector).
//' @param k Moment order (type: double)
//' @param r Artificial sample size (type: int). Default is N^0.8.
//' @param psi Pescaling moment (type: double). Must be <k. Default is 2.0.
//' @param u Sampling range width for sampling range [-u, u] (type: double) Default is 1.0.
//' @param force_random_variate_sample If True, draw random variates for xi and u_series. If False, use quantile function values from a regular percentile space grid. This represents the density function better. Defaiult is False.
//' @param ignore_errors Ignore errors caused by Inf and NaN results for too large absolute moments. If True, it will return test statistic=NA, pvalue=1. If False, it will stop with an error. Default is False. But normally this will indicate an infinite moment.
//' @param verbose If True, print detailed output for debugging. Default is False.
//' @param random_salting Salt number to be added to the random seed (type: int). This prevents identical random variate series if multiple instances are started and run in parallel. Default is 0.
//' @return Trapani's Theta test statistic (type: double).
//' @return Corresponding p-value (Chi^2(1) percentile) (type: double).
//' @examples
//' rvs <- stabledist::rstable(100000, 1.9, 0.5, 1, 0, pm = 0)
//' result <- finite_moment_test(rvs, 2)
//' @export
// [[Rcpp::export]]
arma::vec finite_moment_test(arma::vec obs, 
                             double k, 
                             unsigned int r=0, 
                             double psi=2, 
                             double u=1.0, 
                             bool force_random_variate_sample=0, 
                             bool ignore_errors=0, 
                             bool verbose=0, 
                             int random_salting=0) {
    /* Computes Trapani's (2016) finite moment test.
     * Arguments: 
     *      obs (armadillo numeric vector): Observations
     *      k (double): moment order
     *      r (int): artificial sample size. Default ovewridden with Trapani's recommended N^0.8 below.
     *      psi (double): rescaling moment. Must be <k. Default following Trapani 2016 2.0
     *      u (double): Sampling range width for sampling range [-u, u] (Trapani 2016, Spep 3).
     *                   Default following Trapani u=1.
     *      force_random_variate_sample (bool): Draw random variates for xi and u_series. Default is using quantile 
     *                                          function values from a regular percentile space grid. This represents
     *                                          the density function better.
     *      ignore_errors (bool): Ignore errors caused by Inf and NaN results for too large absolute moments. If True, 
     *                            it will return test statistic NA, pvalue 0. If False, it will stop with an error. 
     *                            Default is False. But normally this will indicate an infinite moment.
     *      verbose (bool): Print detailed output for debugging.
     *      random_salting (int): salt number to be added to the random seed. Prevents identical random variate series 
     *                             if multiple instances are started and run in parallel.
     * Returns:
     *      (armadillo numeric vector)(2):  (0) Trapani's \Theta test statistic,
     *                                      (1) Corresponding Chi^2(1) percentile
     */
    
    // Remove nan's
    obs = obs.elem(arma::find_finite(obs));
    int N = obs.size();
    if (verbose) Rcpp::Rcout << "obs is: " << obs.t() << std::endl;


    // Handle parameter values
    // Artificial sample size r
    if (r <= 0) {               // Handle default value (which is 0)
        // Use Trapani's suggestion of r = N^(4/5)
        r = (int)(pow(N, 0.8));
    }
    if (verbose) Rcpp::Rcout << "r is: " << r << std::endl;
    
    // Rescaling moment psi
    if (psi >= k) {
        psi = k - 1;
        if (psi <= 0) {
            psi = k / 2.;
        }
        if (verbose) Rcpp::Rcout << "Invalid psi value; chose instead psi=" << psi << std::endl;
    }
    
    // Sampling range width u
    if (u <= 0) {
        u = 1.0; 
        if (verbose) Rcpp::Rcout << "Invalid u value; chose instead u=" << u << std::endl;
    }
    
    
    // Compute absolute moment and rescale (Trapani 2016 Step 1)
    
    long double mu = compute_absolute_moment(obs, k);
    long double mu_psi = compute_absolute_moment(obs, psi);
    // rescaling to mu* (Trapani 2016 eq 16)
    if (verbose) Rcpp::Rcout << "mu is: " << mu << std::endl;
    mu = mu / pow(((long double)mu_psi), ((long double)(k / psi)));
    if (verbose) Rcpp::Rcout << "   ...rescaled to: " << mu << std::endl;
    long double rescaling_factor_2 = pow(overloaded_std_norm_moment(psi), k / psi) / overloaded_std_norm_moment(k);
    // rescaling to mu** (Trapani 2016 eq 17)
    mu = mu * rescaling_factor_2;
    if (verbose) Rcpp::Rcout << "   ...rescaled to: " << mu << std::endl;
    
    long double long_exp_mu_half = exp(mu/2);
    double exp_mu_half = long_exp_mu_half;
    if (boost::math::isinf(exp_mu_half)) {
        //Stop
        Rcpp::Rcout << "Error: Rescaled absolute moment is too large. exp(mu/2) cannot be represented as double, which we need to do in the armadillo vector." << std::endl;
        Rcpp::Rcout << "            However, at this kind of value, you can safely assume that your moment is infinite. The test would return p=0 for H0 (moment finite)." << std::endl;
        Rcpp::Rcout << "            Absolute moment mu was in long double: " << long_exp_mu_half << ". In double it was: " << exp_mu_half << std::endl;
        if (ignore_errors) {
            //arma::vec return_values = {NAN, 1.0};
            arma::vec return_values(2);
            return_values << NAN << 1.0;
            return(return_values);
        } else {
            Rcpp::stop("Infinite value detected!\n");
        }
    }
        
    
    // Initialize random number generator
    boost::random::mt19937 generator(time(0) + random_salting); 
    boost::normal_distribution<double> norm_dist(0, 1);
    boost::uniform_real<double> unif_dist(-u, u);
    boost::math::normal standard_normal_frozen(0, 1);
    

    // Obtain xi variates from normal (Trapani 2016 Step 2)
    arma::vec xi_bak(r);
    double xi;
    arma::vec xim(r);
    for (unsigned int i = 0; i < r; i++) {
        if (force_random_variate_sample) {      
            // Use sample of random variates as suggested in (Trapani 2016 Step 2).
            xi = norm_dist(generator);
        } else {
            // Do not use sample of random variates as suggested in (Trapani 2016 Step 2).
            // Instead using a quantile function values from a recular percentile space grid.
            // This represents the density function better.
            xi = quantile(standard_normal_frozen, (i + 1.0) / (r + 1.0));
        }
        //xim(i) = pow(exp(mu), 0.5) * xi;
        //xim(i) = exp(mu/2) * xi;
        xim(i) = exp_mu_half * xi;
        xi_bak(i) = xi;
    }
    if (verbose) {
        xi_bak = arma::sort(xi_bak);
        xim = arma::sort(xim);
        Rcpp::Rcout << "xi is: " << xi_bak.t() << std::endl;
        Rcpp::Rcout << "xi modified is: " << xim.t() << std::endl;
    }
    
    // Obtain sorted u variates from uniform (Trapani 2016 Step 3)
    arma::vec u_series(r);
    for (unsigned int i = 0; i < r; i++) {
        if (force_random_variate_sample) {      
            // Use sample of random variates from some bounded density function (here uniform) 
            // as suggested in (Trapani 2016 Step 3).
            u_series(i) = unif_dist(generator);
        } else {
            // Do not use sample of random variates as suggested in (Trapani 2016 Step 3).
            // Instead using a quantile function values from a recular percentile space grid.
            // This represents the density function better.
            u_series(i) = 2.0 * u * i / (r - 1.0) - u;
        }
    }
    // Sort
    u_series = arma::sort(u_series);
    if (verbose) Rcpp::Rcout << "u series is: " << u_series.t() << std::endl;


    // Compute zeta (Trapani 2016 eq 6/7)
    // For efficiency, we only deal with counts of values, not with long binary arrays as Trapani 2016 eq 6.
    arma::uvec zeta = build_zeta_sorted(xim, u_series);
    if (verbose) Rcpp::Rcout << "zeta is: " << zeta.t() << std::endl;

    // Compute theta (Trapani 2016 eq 7)
    arma::vec theta(r);
    for (unsigned int i = 0; i < r; i++) {
        theta(i) = 2 / pow(r, 0.5) * (zeta(i) - 0.5 * r);
    }
    if (verbose) Rcpp::Rcout << "theta is: " << theta.t() << std::endl;
    
    // Compute Theta (Trapani 2016 eq 8)
    // There seem to be two ways to do this.
    // Option 1: Sum and divide by r. With sufficiently large r, this should be sufficiently close to the integral in Trapani 2016 eq 8
    double Theta = 0.0;
    for (unsigned int i = 0; i < r; i++) {
        Theta += pow(theta(i), 2);
        //Rcpp::Rcout << "   Theta is: " << Theta << std::endl;
    }
    Theta /= r;
    if (verbose) Rcpp::Rcout << "Theta is: " << Theta << std::endl;

    // Trapezoid Theta
    // Option 2: Sum over trapezoids for any two adjacent values of u and corresponding theta
    double TTheta = 0.0;
    if (verbose) {
        for (unsigned int i = 1; i < r; i++) {
            // Internal u are each counted twice half 
            TTheta += ((pow(theta(i-1), 2) + pow(theta(i), 2)) / 2.) * (u_series(i) - u_series(i-1));
            //Rcpp::Rcout << "   TTheta is: " << TTheta << std::endl;
        }
        // Left and right margin are further counted once half
        TTheta += (pow(theta(0), 2) / 2.) * (u_series(0) - (-u));
        TTheta += (pow(theta(r-1), 2) / 2.) * (u - u_series(r-1));
        // And finally we rescale by width of the u sampling range 
        TTheta /= 2 * u;
        Rcpp::Rcout << "Trapezoid Theta is: " << TTheta << std::endl;
    }
    
    // Obtain chi^2 percentile
    double chisq1_percentile = get_chisq1_percentile(Theta);
    
    if (verbose) Rcpp::Rcout << "Theta Chi^2(1) percentile is: " << chisq1_percentile << std::endl;
    if (verbose) {
        double TT_chisq1_percentile = get_chisq1_percentile(TTheta);
        Rcpp::Rcout << "Trapezoid Theta Chi^2(1) percentile is: " << TT_chisq1_percentile << std::endl;
    }
    
    // Return
    //arma::vec return_values = {Theta, chisq1_percentile};
    arma::vec return_values(2);
    return_values << Theta << chisq1_percentile;
    return(return_values);
}

