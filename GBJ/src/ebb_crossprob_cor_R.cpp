// ebb_crossprob_cor.c
// Ryan Sun, Harvard University
// July 17, 2016

// This binary calculates an approximation to the p-value of GOF tests.

// Usage is ./ebb_crossprob_cor [NUM_STATS] [BOUNDS_FILE] [CORRELATIONS_FILE]
// where NUM_STATS is the size of your set, BOUNDS_FILE is the name of the file
// holding the boundaries (0<=b_1<=...<=b_NUM_STATS) which come from inversion
// of the GOF statistic, and CORRELATIONS_FILE is the name of the file holding
// all the NUM_STATS*(NUM_STATS-1)/2 pairwise correlations between the observations
// in your set.
// Both BOUNDS_FILE and CORRELATIONS_FILE should be a single column of numbers
// with no headers/labels/blank spaces/other characters.

// One difference from before is that we calculate the conditional moments upfront in
// a separate routine.
// We are also now using the EBB distribution of Prentice (1988) instead of the standard
// Beta-Binomial.

// The routines below are:
// (1) avg_cond_covar - Calculate the average conditional covariance between any two indicators
//  Y_i,Y_j where Y_i=P(|Z_i|>=t_k).
// (2) match_moments - Given an average conditional covariance and a conditional mean, calculate
// \lambda and \gamma for the EBB PMF.  Also check to ensure both are inside the allowable parameter space.
// (3) eval_EBB_PMF - Self-explanatory
// (4) calc_qka - Called to fill each entry of the d*d matrix leading to final
// p-value.  Sums over all m>=a in P(S(t_k)=a|S(t_k)=m).
// (5) calc_allq - Loop through each entry of the d*d matrix until we get the (d,1) element.

// Need this for the boost functions on windows machines
// [[Rcpp::depends(BH)]]
#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>
#include <numeric>
#include <boost/math/special_functions/erf.hpp>		// For normal PDF/CDF
#include <boost/math/constants/constants.hpp>			// For pi()
#include <Rcpp.h>

using namespace Rcpp;

inline double fast_square(const double &x);
inline double fast_cube(const double &x);
inline double dnorm(const double &x);
inline double surv(const double &x);

bool create_logftable(const int &d, std::vector<double> &log_ftable);
bool create_rtable(const int &d,
                   const std::vector<double> &r_vec,
                   std::vector< std::vector<double> >&rtable);
bool create_hermtable(const int &d,
                      const std::vector <double> &t_vec,
                      std::vector<std::vector<double> > &herm_table);
bool avg_cond_covar(const int &d,
                    const std::vector<double> &t_vec,
                    const std::vector<std::vector<double> > &r_table,
                    const std::vector<std::vector<double> > &herm_table,
                    std::vector<double> &covar_vec);
bool eval_EBB_PMF_allN(const int &max_n,
                       int &y,
                       const double &lambda,
                       const double &gamma,
                       std::vector<double> &PMF_vec);
double calc_qka(const int &d,
                const int &k,
                const int &a,
                const std::vector<double> &prev_row,
                const std::vector<double> &log_ftable,
                const bool &ind_flag,
                const double &lambda,
                const double &gamma);
double calc_allq(const int &d,
                 const std::vector <double> &t_vec,
                 const std::vector <double> &r_vec,
                 const bool &ind_flag);


////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////


// To square a number faster than pow().
inline double fast_square(const double &x)
{
    return (x*x);
}

// To cube a number faster than pow().
inline double fast_cube(const double &x)
{
    return (x*x*x);
}

// PDF of a standard normal RV.
inline double dnorm(const double &x)
{
    return (exp(-fast_square(x) / 2.0) /
            sqrt(2.0*boost::math::constants::pi<long double>()));
}

// Survival function of a standard normal, 1-F_x(x)
inline double surv(const double &x)
{
    return (1 - 0.5 * erfc(-x/sqrt(2.0)));
}



////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////


// For the independence case, good to have lookup table of log-factorials
// Assume d doesn't go past the limit of 'int,' which is 32,767.
bool create_logftable(const int &d, std::vector<double> &log_ftable)
{
    // 0! = 1! = 1
    log_ftable[0] = 0.0;
    log_ftable[1] = 0.0;

    for (int iii=2; iii<=d; ++iii)
    {
        log_ftable[iii] = lgamma(iii+1);
    }

    return 0;
}

// Create a d(d-1)/2 by 5 table of each \rho_ij to the 2/4/6/8/10 power.
// These are reused d times in the calculation of the variance of S(t).
bool create_rtable(const int &d,
                   const std::vector <double> &r_vec,
                   std::vector< std::vector<double> > &r_table)
{
    long double num_cor = d*(d-1)/2;
    std::vector <double> r_two(num_cor);
    std::vector <double> r_four(num_cor);
    std::vector <double> r_six(num_cor);
    std::vector <double> r_eight(num_cor);
    std::vector <double> r_ten(num_cor);

    // Long is approximately 2*10^9 which is enough to cover up to d~40000.
    for (long iii=0; iii<num_cor; ++iii)
    {
        r_two[iii] = fast_square(r_vec[iii]);
        r_four[iii] = pow(r_vec[iii], 4.0);
        r_six[iii] = pow(r_vec[iii], 6.0);
        r_eight[iii] = pow(r_vec[iii], 8.0);
        r_ten[iii] = pow(r_vec[iii], 10.0);
    }

    r_table[0] = r_two;
    r_table[1] = r_four;
    r_table[2] = r_six;
    r_table[3] = r_eight;
    r_table[4] = r_ten;

    return 0;
}


// Create a d*5 matrix holding the values of the 2/4/6/8/10th Hermite polynomials
// evaluated at t_k from k=1,...,d.
// Each value will be reused d(d-1)/2 times in the calculation of the variance of S(t).
// We also pre-divide each of the He_(2n-1) terms by (2n!) to save a little more time.
bool create_hermtable(const int &d,
                      const std::vector <double> &t_vec,
                      std::vector<std::vector<double> > &herm_table)
{
    /*
     // Could check here if same t, then don't need to do the calculations
     */

    // Hold each of the five Hermite polynomial terms.
    std::vector <double> herm_one(d);
    std::vector <double> herm_two(d);
    std::vector <double> herm_three(d);
    std::vector <double> herm_four(d);
    std::vector <double> herm_five(d);

    for (int iii=0; iii<d; ++iii)
    {
        herm_one[iii] = fast_square( t_vec[iii] ) / 2.0;
        herm_two[iii] = fast_square( fast_cube(t_vec[iii]) - 3.0*t_vec[iii] ) / 24.0;
        herm_three[iii] = fast_square( pow(t_vec[iii], 5) - 10.0*fast_cube(t_vec[iii]) +15.0*t_vec[iii] ) / 720.0;
        herm_four[iii] = fast_square( pow(t_vec[iii], 7) - 21.0*pow(t_vec[iii], 5)
                                     + 105.0*fast_cube(t_vec[iii]) - 105.0*t_vec[iii] ) / 40320.0;
        herm_five[iii] = fast_square( pow(t_vec[iii], 9) - 36.0*pow(t_vec[iii], 7)
                                     + 378.0*pow(t_vec[iii], 5) - 1260.0*fast_cube(t_vec[iii]) + 945.0*t_vec[iii] ) / 3628800.0;
    }

    // Populate first dimension.
    herm_table[0] = herm_one;
    herm_table[1] = herm_two;
    herm_table[2] = herm_three;
    herm_table[3] = herm_four;
    herm_table[4] = herm_five;

    return 0;
}



// Calculate the average conditional covariance term, P(|Z_i|,|Z_j| >= t_k)
// for k=1,...,d.
bool avg_cond_covar(const int &d,
                    const std::vector<double> &t_vec,
                    const std::vector< std::vector<double> > &r_table,
                    const std::vector< std::vector<double> > &herm_table,
                    std::vector< double > &covar_vec)
{

    /*
     // Could actually do this for all k=1:d and a=1:k.
     // Right now just do the same for all k,a
     // Right now indexed by t_1,....,t_d
     // But it wouldn't be so much harder, just keep the terms we have now
     // and do some additional multiplications/divisions.
     // O(d^2) instead of O(d) but the operations are very simple multiplications)
     */


    // Numerator and denominator vectors for the tricky infinite sum in the variance.
    double num_cor = d*(d-1)/2;
    double cond_mean_sq = 0.0;
    std::vector<double> numerator_vec(num_cor);
    std::vector<double> denominator_vec(num_cor);

    // First the unconditional covariance at k=1.
    double sum_term = 0.0;
    double surv_term = 4 * fast_square( surv(t_vec[0]) );
    double phi_term = 4 * fast_square( dnorm(t_vec[0]) );

    /*
     // Could pull out the surv_term and phi_term from the sum, don't have to add them every time
     */


    for (long cor_it=0; cor_it<num_cor; ++cor_it)
    {
        numerator_vec[cor_it] = surv_term + phi_term *
        (r_table[0][cor_it]*herm_table[0][0] +
         r_table[1][cor_it]*herm_table[1][0] +
         r_table[2][cor_it]*herm_table[2][0] +
         r_table[3][cor_it]*herm_table[3][0] +
         r_table[4][cor_it]*herm_table[4][0]);
        sum_term += numerator_vec[cor_it];
    }
    cond_mean_sq = fast_square( 2*surv(t_vec[0]) );
    covar_vec[0] = ( sum_term - num_cor*cond_mean_sq ) / num_cor;

    // Now calculate the conditional p_n(1) and p_n(2) at k=2,...,d.
    for (int kkk=2; kkk<=d; ++kkk)
    {
        // For a test like GBJ the smallest half t_k are the same, so we don't need to call
        // calc_qka or calculate conditional moments for first half of rows.
        // Make sure our tolerance here 10^-8 is the same as in calc_allq.
        if ( fabs(t_vec[kkk-1] - t_vec[kkk-2]) < pow(10.0, -8.0) )
        {
            covar_vec[kkk-1] = covar_vec[kkk-2];
            continue;
        }

        // The conditional variance at t_k.
        denominator_vec.swap(numerator_vec);
        sum_term = 0.0;
        surv_term = 4 * fast_square( surv(t_vec[kkk-1]) );
        phi_term = 4 * fast_square( dnorm(t_vec[kkk-1]) );
        for (long cor_it=0; cor_it<num_cor; ++cor_it)
        {
            numerator_vec[cor_it] = surv_term + phi_term *
            (r_table[0][cor_it]*herm_table[0][kkk-1] +
             r_table[1][cor_it]*herm_table[1][kkk-1] +
             r_table[2][cor_it]*herm_table[2][kkk-1] +
             r_table[3][cor_it]*herm_table[3][kkk-1] +
             r_table[4][cor_it]*herm_table[4][kkk-1]);
            sum_term += numerator_vec[cor_it] / denominator_vec[cor_it];
        }

        // Record
        cond_mean_sq = fast_square( surv(t_vec[kkk-1]) / surv(t_vec[kkk-2]) );
        covar_vec[kkk-1] = ( sum_term - num_cor*cond_mean_sq ) / num_cor;
    }

    return 0;
}


// Evaluate the EBB PMF for a range of n, so we don't have to repeat
// multiplications for Pr[S(t_k)=a|S(t_k-1)=m] for m=a:(d-k+1).  Here 'y' is a.
bool eval_EBB_PMF_allN(const int & max_n,
                       const int & y,                  // min_n = y
                       const double & lambda,
                       const double & gamma,
                       std::vector<double> & PMF_vec)
{
    // If (d-k+1) <=1 then we don't need to do these calculations
    // because we will be using straight binomial calculation for
    // P(S(t)=a|m=0 or 1).
    if (max_n < 2) { return 0;}

    double log_prob_mass = 0.0;

    // We only need to fill PMF_vec from a (or 'y' as it's called here)
    // to d-k+1.
    // The EBB pmf has a product with a terms, a product with d-a terms, and a product
    // with d terms.
    // We don't need to fill PMF_vec[0] or PMF_vec[1] ever because P(S(t)=a|m=0/1) always
    // uses straight binomial.
    // We don't need to fill PMF_vec for any indices less than 'y' because m has to be
    // greater than a.

    // The case of y=(d-k+1) has only two products, treat it differently
    // The case of y=0 is special because then we only have two products; treat it differently.
    if (y == max_n)
    {
        for (int jjj=0; jjj<max_n; ++jjj)
        {
            log_prob_mass += std::log(lambda+gamma*jjj) - std::log(1+gamma*jjj);
        }

        // Only need to fill the one entry of PMF_vec if a=(d-k+1)
        PMF_vec[y] = log_prob_mass;
        return 0;

    } else if (y == 0)          // Again only two products
    {
        for (int jjj=0; jjj<max_n; ++jjj)
        {
            log_prob_mass += std::log(1 - lambda + gamma*jjj) - std::log(1 + gamma*jjj);
            PMF_vec[jjj+1] = log_prob_mass;
        }
        return 0;

    } else {                    // Normal three products
        // Get two of the products here, the (v-1) and (part of) the (d-1).
        // Don't need to fill PMF_vec for this part since m>=a always.
        for (int jjj=0; jjj<y; ++jjj)
        {
            log_prob_mass += std::log(lambda + gamma*jjj) - std::log(1 + gamma*jjj);
        }
        PMF_vec[y] = log_prob_mass;
        // Get the rest of (d-1) and also the (d-v-1) product
        for (int jjj=y; jjj<max_n; ++jjj)
        {
            log_prob_mass += std::log(1 - lambda + gamma*(jjj-y)) - std::log(1 + gamma*jjj);
            PMF_vec[jjj+1] = log_prob_mass;
        }
    }

    return 0;
}


// Calculate q_k,a by summing over q_k,a|S(t)=m for m=a:(d-k)
double calc_qka(const int &d,
                const int &k,
                const int &a,
                const std::vector<double> &prev_row,
                const std::vector<double> &log_ftable,
                const bool &ind_flag,
                const double &lambda,
                const double &gamma)
{
    double min_gamma;
    double log_m_choose_a;
    std::vector<double> PMF_vec(d+1);
    std::vector<double> log_qka_calc_vec(d+k-1 - a + 1);

    // To hold the the EBB probability (w/o factorial part) for m=a:(d-k+1)
    if (!ind_flag)
    {
        bool all_EBB_status = eval_EBB_PMF_allN((d-k+1),
                                                a,
                                                lambda,
                                                gamma,
                                                PMF_vec);
    }

    // Sum over all possible values of S(t)=m
    double max_value = 0.0;
    for (int mmm=a; mmm<=(d-k+1); ++mmm)
    {
        log_m_choose_a = log_ftable[mmm] - log_ftable[a] - log_ftable[mmm-a];

        // Use binomial if m=0/1 or if independence flag or if gamma outside parameter space.
        min_gamma = std::max(-lambda/(mmm-1), -(1-lambda)/(mmm-1));
        if (mmm<=1 || ind_flag || gamma<min_gamma)
        {
            log_qka_calc_vec[mmm-a] = prev_row[mmm] + log_m_choose_a + a*std::log(lambda) + (mmm-a)*std::log(1.0-lambda);
        }
        else          // EBB PMF
        {
            log_qka_calc_vec[mmm-a] = prev_row[mmm] + log_m_choose_a + PMF_vec[mmm];
        }

        // Check for max value
        if (log_qka_calc_vec[mmm-a] > max_value) {
            max_value = log_qka_calc_vec[mmm-a];
        }
    }

    // Different for k=1 because no way to represent log(0)
    if (k == 1)
    {
        return log_qka_calc_vec[d-k+1-a];
    }

    // LogSumExp procedure
    double sum_term = 0.0;
    for (int mmm=a; mmm<=(d-k+1); ++mmm)
    {
        sum_term += std::exp(log_qka_calc_vec[mmm-a] - max_value);
    }
    return max_value + std::log(sum_term);
}




// The p-value calculation 'master' function, the interface between the math and main fn.
// Loop through all k, calculating q_k,a from a=0:(t-k) until we get to q_d,0/
// Return the p-value.
double calc_allq(const int &d,
                 const std::vector <double> &t_vec,
                 const std::vector <double> &r_vec,
                 const bool &ind_flag)
{
    // Fill the d(d-1)/2 by 5 table of correlations to the 2/4/6/8/10 power.
    // Fill the d*5 table of Hermite polynomial terms in the S(t) conditional variance calculation.
    long num_cor = d*(d-1)/2;
    std::vector<std::vector<double> > r_table(5, std::vector<double>(num_cor));
    std::vector<std::vector<double> > herm_table(5, std::vector<double>(d));
    bool filled_hermtable = create_hermtable(d,
                                             t_vec,
                                             herm_table);
    bool filled_rtable = create_rtable(d,
                                       r_vec,
                                       r_table);

    if (filled_hermtable || filled_rtable)
    {
        // We can't use cout or exit() in Rcpp - shouldn't come here ever anyway
        return(-1);
        //std::cout << "Problem creating lookup tables." << std::endl;
        //exit(1);
    }

    // Fill the d*1 vector of log-factorials
    std::vector<double> log_ftable(d+1);
    bool filled_logftable = create_logftable(d,
                                             log_ftable);

    // Calculate the conditional average pairwise correlation for k=1,...,d.
    std::vector<double> covar_vec(d);
    bool filled_condcovar = avg_cond_covar(d,
                                           t_vec,
                                           r_table,
                                           herm_table,
                                           covar_vec);
    if (filled_condcovar)
    {
        // We can't use cout or exit() in a Rcpp - shouldn't come here ever anyway
        return(-1);
        //std::cout << "Problem calculating moments." << std::endl;
        //exit(1);
    }

    // We don't actually need to hold the entire d*d matrix, just have to always
    // know the previous row.
    std::vector<double> prev_row(d+1);
    std::vector<double> current_row((d+1), 0.0);

    // Initial conditions.
    double prev_bound = 0.0;
    current_row[d] = 0.0;

    // Loop through k=1,...,d.
    double lambda = 1.0;
    double gamma = 1.0;
    double rho;
    double avg_cond_covar;
    double temp_qka;
    for (int kkk=1; kkk<=d; ++kkk)
    {
        prev_row = current_row;

        // If new bound same as previous bound (make sure tolerance is same as match_moments).
        if (fabs(t_vec[kkk-1] - prev_bound) < pow(10.0, -8.0))
        {
            current_row[d-kkk+1] = 0.0;             // We don't actually need to zero it out since we never get here.
            prev_bound = t_vec[kkk-1];
            continue;
        }

        // New bound, new probabilities for the row.
        std::fill(current_row.begin(), current_row.end(), 0.0);

        // Match moments once for each row
        lambda = surv(t_vec[kkk-1]) / surv(prev_bound);
        avg_cond_covar = covar_vec[kkk-1];
        rho = avg_cond_covar / (lambda*(1-lambda));
        gamma = rho / (1-rho);

        // For each k, we want q_k,a for a=0:(t-k).
        for (int aaa=0; aaa<=(d-kkk); ++aaa)
        {
            temp_qka = calc_qka(d,
                                kkk,
                                aaa,
                                prev_row,
                                log_ftable,
                                ind_flag,
                                lambda,
                                gamma);
            current_row[aaa] = temp_qka;

        }

        // Update so we can check if same as next bound, also for matching lambda.
        prev_bound = t_vec[kkk-1];
    }

    return (1-std::exp(current_row[0]));
}



////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
// Main function reads in our boundary points and correlations, the
// runs calc_allq and prints the p-value.
// If the `correlations` argument is less than -1, then this will trigger
// the independence flag.

// [[Rcpp::export]]

double ebb_crossprob_cor_R(int d, NumericVector bounds, NumericVector correlations) {

    // For the Rcpp version, just input the bounds as NumericVectors and transfer them
    // to std::vectors.
    long num_cor = d*(d-1)/2;
    std::vector<double> boundary_pts;
    std::vector<double> cors_vec;

    // Reserve memory to prevent fragmentation, easy since we know the exact size
    boundary_pts.reserve(d);

    // Put the boundary pts into array, one by one.
    // Should be sorted in increasing order.
    for (int iii=0; iii<d; ++iii)
    {
        boundary_pts.push_back(bounds[iii]);
    }

    // If the last argument is -999 then no need for correlation vector
    // Carry through the independence flag so that we don't ever need to matchMoments()
    bool indFlag = false;
    if (correlations[1] < -1)
    {
        indFlag = true;
        cors_vec.push_back(-1.0);
    }
    else
    {
        // Put the correlations into array, one by one.
        cors_vec.reserve(num_cor);
        for (long iii=0; iii<num_cor; ++iii)
        {
            cors_vec.push_back(correlations[iii]);
        }
    }

    // Calculate the p-value.
    double pvalue = calc_allq(d, boundary_pts, cors_vec, indFlag);

    // Return the p-value.
    return pvalue;
}




