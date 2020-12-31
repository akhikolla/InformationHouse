#include <Rcpp.h>        // Rcpp
#include <vector>        // vector, size_t
#include <numeric>       // iota, accumulate
#include <algorithm>     // sort, lower_bound, upper_bound
#include <cmath>         // abs(double), sqrt, |Do not use pow! It's SLOW!|
#include "geos_stuff.h"  // sort_indices, permutate_using_index



/* Calculate PCF
 *-----------------------------------------------------------------------------
 * calculate PCF from dists and ratios at all rs using the Epanechnikov kernel
 * with stoyan giving the smoothing and Ripley's edge correction
 *
 * requirements:
 *  dist and ratio are in the same order & of the same length
 *  rs is not empty
 *  area, lambda, stoyan must be plausible
 *
 * returns vector PCF being of the same length as rs
 */
std::vector<double>
do_pcf(std::vector<double> dist,
       std::vector<double> ratio,
       const std::vector<double> rs,
       const double area, const double lambda, const double stoyan)
{
    double omega, edge_cor, denominator;
    double delta = stoyan / sqrt(lambda);  // bandwith
    std::vector<double>::iterator db, de, pb, lwr, upr, iter;
    std::vector<double> diff4k, ratio4k, pcf(rs.size());

    // sort dist and ratio according to dist
    std::vector<size_t> index = sort_indices(dist, true);
    permutate_using_index(dist, index);
    permutate_using_index(ratio, index);
    db = dist.begin();
    de = dist.end();
    pb = ratio.begin();

    for(unsigned int i=0; i < rs.size(); i++){
        // distances within bandwith of kernel
        diff4k.clear();
        lwr = std::lower_bound(db, de, rs[i]-delta);
        upr = std::upper_bound(lwr, de, rs[i]+delta);
        for(iter = lwr; iter != upr; ++iter){
            diff4k.push_back(std::abs(*iter - rs[i]));
        }

        // copy the same range of vector ratio
        std::vector<double> ratio4k(pb+(lwr-db), pb+(upr-db));
        unsigned int n = diff4k.size();

        if(n > 0){
            // weighted distances (kernel function)
            omega = 0;
            edge_cor = 0;
            for(unsigned int k = 0; k < n; k++){
                omega += (3/(4*delta)) *
                    (1 - ((diff4k[k]*diff4k[k]) / ( delta*delta)));
                edge_cor += ratio4k[k];

            }
            denominator = (lambda * lambda) * 2 * M_PI * rs[i] * area * edge_cor/n;
            // w/o edge correction: (lambda**2) * 2 * pi * rs[i] * area

            pcf[i] = (omega / denominator);
        }
        else
        {
            pcf[i] = 0;
        }
    }

    return pcf;
}



/* Create Envelope
 *-----------------------------------------------------------------------------
 * get nrank smallest and largest value at every rs
 * from a column-first-ordered matrix stored as a vector
 *
 * requirements:
 *  nrs = rs.size()
 *  nrank must be >=1 && <= rs.size()/2
 *  gnull.size() must be a multiple of rs.size()
 *  gnull.size() must be at least rs.size()
 *  rs.size() must not be empty
 *  lwr, upr, mean must be of size rs.size()
 *
 * returns lwr, upr, mean via arguments passed in by reference
 */
void
do_env(const std::vector<double> gnull,
       const unsigned int nrs, const int nrank,
       std::vector<double> &lwr,
       std::vector<double> &upr,
       std::vector<double> &mean)
{
    unsigned int nsim = gnull.size() / nrs;
    if((nrank < 1) || (nrank >= nsim/2.0))
        throw std::range_error("n_rank must be >= 1 and < n_sim/2");
    double sum;
    std::vector<double> row(nsim);

    // iterate over rows
    for(unsigned int m = 0; m < nrs; m++){
        // iterate over cols
        for(unsigned int n = 0; n < nsim; n++){
            // assemble row
            row[n] = gnull[(m + n*nrs)];
        }

        // fetch nrank smallest/largest element
        if(nrank == 1){
            // this is always fastest
            auto mima = std::minmax_element (row.begin(), row.end());
            lwr[m] = *mima.first;
            upr[m] = *mima.second;
        }
        else
        {
            // if vector.size() < _ISORT_MAX (== 32)
            //   nth_element() does a complete sort and in this case TWICE
            // cannot access _ISORT_MAX so hard code 32
            if (nsim < 32){
                std::sort(row.begin(), row.end());
                lwr[m] = row[nrank-1];
                upr[m] = row[nsim-nrank];
            }
            else
            {
                nth_element(row.begin(), row.begin()+(nrank-1), row.end());
                lwr[m] = row[nrank-1];
                nth_element(row.begin(), row.begin()+(nsim-nrank), row.end());
                upr[m] = row[nsim-nrank];
            }
        }

        // calculate mean for bias correction
        sum = std::accumulate(row.begin(), row.end(), 0.0);
        mean[m] = sum/nsim;
    }

    return;
}



/* Fractionize, PCF, Envelope
 *-----------------------------------------------------------------------------
 * splits up the vectors sim, dist, ratio in homogeneous ranges of sim
 * calculate pcf for all ranges from dist & ratio at all rs
 * create envelope (lwr,upr) from all sims != 0 using nrank
 *
 * requirements:
 *  sim, dist, ratio are of same length and in the same order
 *  sdp contains at least 1 range with sim == 0
 *     AND >= 2 ranges with sim != 0
 *  rs not empty
 *  nrank must be >=1 && <= rs.size()/2
 *  area, lambda, stoyan must be plausible
 *  pcf, lwr, must be of size rs.size()
 *
 * returns nsim, pcf, lwr, upr via arguments passed in by reference
 */
void
fractionize(const std::vector<double> sim,
    const std::vector<double> dist,
    const std::vector<double> ratio,
    const std::vector<double> rs,
    const double area, const double lambda,
    const double stoyan, const int nrank,
    unsigned int &nsim,
    std::vector<double> &pcf,
    std::vector<double> &lwr,
    std::vector<double> &upr)
{
    unsigned int nrs = rs.size();
    double sim_ind;
    std::vector<double> gnull, t_g, t_dist, t_ratio, mean(nrs);

    // setup iterators
    std::vector<double>::const_iterator sb, se, db, pb, lo, hi;
    sb = sim.begin();
    se = sim.end();
    db = dist.begin();
    pb = ratio.begin();

    // go over vectors and copy out homogeneous ranges of sim
    hi = sb;
    nsim = 0;
    while(hi != se){
        lo = hi;
        sim_ind = *lo;
        hi = std::upper_bound(lo, se, sim_ind);
        // copy range
        t_dist.clear();
        t_dist.insert(t_dist.end(), db+(lo-sb), db+(hi-sb));
        t_ratio.clear();
        t_ratio.insert(t_ratio.end(), pb+(lo-sb), pb+(hi-sb));

        // calc pcf
        t_g = do_pcf(t_dist, t_ratio, rs, area, lambda, stoyan);
        if(sim_ind == 0)
        {
            pcf = t_g;
        }
        else
        {
            nsim++;
            gnull.insert(gnull.end(), t_g.begin(), t_g.end());
        }
    }

    do_env(gnull, nrs, nrank, lwr, upr, mean);

    // bias correction
    for(unsigned int i = 0; i < nrs; i++){
        pcf[i] = pcf[i] / mean[i];
        lwr[i] = lwr[i] / mean[i];
        upr[i] = upr[i] / mean[i];
    }
    return;
}



/* Wrap for R
 *-----------------------------------------------------------------------------
 * wrap fractionize fro usage by R
 *
 * requirements:
 *  sim, dist, ratio are of same length and in the same order
 *  sdp contains at least 1 range with sim == 0
 *     AND >= 2 ranges with sim != 0
 *  rs not empty
 *  nrank must be >=1 && <= rs.size()/2
 *  area, lambda, stoyan must be plausible
 *  pcf, lwr, must be of size rs.size()
 *
 * returns a dataframe
 */

// [[Rcpp::export]]
Rcpp::DataFrame
pcf_envelope(Rcpp::NumericVector sim,
             Rcpp::NumericVector dist,
             Rcpp::NumericVector ratio,
             Rcpp::NumericVector rs,
             double area,
             int nobj,
             double stoyan,
             int nrank)
{
    // convert Rcpp::NumericVector to std::vector<double> ---------------------
    // naive
    // double* psim  = sim.begin();
    // double* pdist = dist.begin();
    // double* pratio = ratio.begin();
    // double* prs   = rs.begin();

    // std::vector<double> vsim(sim.begin(), sim.end());
    // std::vector<double> vdist(dist.begin(), dist.end());
    // std::vector<double> vratio(ratio.begin(), ratio.end());
    // std::vector<double> vrs(rs.begin(), rs.end());

    std::vector<double> vsim   = Rcpp::as< std::vector<double> >(sim);
    std::vector<double> vdist  = Rcpp::as< std::vector<double> >(dist);
    std::vector<double> vratio = Rcpp::as< std::vector<double> >(ratio);
    std::vector<double> vrs    = Rcpp::as< std::vector<double> >(rs);

	// fractionize, PCF, envelope ---------------------------------------------
    unsigned int nsim, nrs = vrs.size();
    double lambda = nobj/area;
    std::vector<double> pcf(nrs), lwr(nrs), upr(nrs);

    fractionize(vsim, vdist, vratio, vrs,
                area, lambda, stoyan, nrank,
                nsim, pcf, lwr, upr);

    // create data.frame ------------------------------------------------------
    Rcpp::DataFrame df = Rcpp::DataFrame::create(
                            Rcpp::Named("r")   = rs,
                            Rcpp::Named("g")   = pcf,
                            Rcpp::Named("lwr") = lwr,
                            Rcpp::Named("upr") = upr
                         );

    df.attr("n_sim")  = nsim;
    df.attr("n_rank") = nrank;
    df.attr("alpha")  = (2.0 * nrank) / (nsim + 1.0);
    df.attr("correc") = "Ripley";
    df.attr("kernel") = "epanechnikov";
    df.attr("stoyan") = stoyan;
    df.attr("bw")     = stoyan/sqrt(lambda);
    df.attr("class")  = Rcpp::CharacterVector::create("fv_pcf", "data.frame");

	return df;
}
