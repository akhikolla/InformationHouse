#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

double hash_cz(std::vector<bool> &cz, const NumericVector &lprimes) {
  double pid = 0;
  unsigned int cz_size = cz.size();
  for(unsigned int i = 0; i < cz_size; i++) {
    if (cz[i] == true) {
      pid += lprimes[i];
    }
  }
  return pid;
}

std::list<std::vector<bool>> link_cz_nb(std::vector<bool> &cz,
                                        const IntegerVector &nb,
                                        std::unordered_set<double> &lz_hash_table,
                                        const NumericVector &lprimes) {
  unsigned int nb_size = nb.length();
  std::list<std::vector<bool>> lz;
  for(unsigned int i = 0; i < nb_size; i++) {
    cz[nb[i]] = true;
    if (lz_hash_table.insert(hash_cz(cz, lprimes)).second) {
      // insertion successful, add cz to list
      lz.emplace_back(cz);
    }
    cz[nb[i]] = false;
  }
  return lz;
}

IntegerVector colsums_sub_iv(IntegerMatrix &cw, IntegerVector sub) {
  unsigned int cw_ncol = cw.ncol();
  unsigned int sub_length = sub.length();
  IntegerVector msum(cw_ncol);
  for (unsigned int i = 0; i < cw_ncol; i++) {
    for (unsigned int j = 0; j < sub_length; j++) {
      msum[i] += cw(sub[j], i);
    }
  }
  return msum;
}

IntegerVector colsums_sub(IntegerMatrix &cw, std::vector<bool> sub) {
  unsigned int cw_ncol = cw.ncol();
  unsigned int sub_size = sub.size();
  IntegerVector msum(cw_ncol);
  for (unsigned int i = 0; i < cw_ncol; i++) {
    for (unsigned int j = 0; j < sub_size; j++) {
      if (sub[j]) {
        msum[i] += cw(j, i);
      }
    }
  }
  return msum;
}

IntegerVector add_biv(std::vector<bool> &lv,
                      IntegerVector iv) {
  unsigned int lv_size = lv.size();
  for (unsigned int i = 0; i < lv_size; i++) {
    if (lv[i] == true) {
      iv[i] += 1;
    }
  }
  return iv;
}

LogicalVector lmb(LogicalVector a, std::vector<bool> &b) {
  unsigned int a_length = a.length();
  for (unsigned int i = 0; i < a_length; i++) {
    if (b[i]) {
      a[i] = false;
    }
  }
  return a;
}

std::list<vector<bool>> csg2_cpp(std::vector<bool> &cz,
                                 IntegerVector &cnn,
                                 IntegerMatrix &cw,
                                 IntegerVector &s,
                                 std::unordered_set<double> &lz_hash_table,
                                 NumericVector &lprimes) {
  LogicalVector ab = (add_biv(cz, colsums_sub(cw, cz)) >= 1);
  return link_cz_nb(cz, s[lmb(ab, cz)], lz_hash_table, lprimes);
}

std::list<std::vector<bool>> lcsg2_cpp(std::list<std::vector<bool>> &lz,
                                       IntegerVector &cnn,
                                       IntegerMatrix &cw,
                                       IntegerVector &s,
                                       NumericVector &lprimes) {
  // initialize needed variables
  std::list<std::vector<bool>> new_lz;
  std::unordered_set<double> lz_hash_table;
  auto lz_end = lz.end();

  for (auto lit = lz.begin(); lit != lz_end; lit++) {
    // *lit has cz (current zone)
    // expand cz to lz and splice to end of new_lz
    new_lz.splice(new_lz.end(), csg2_cpp(*lit, cnn, cw, s, lz_hash_table, lprimes));
  }
  return new_lz;
}

IntegerMatrix sub_cnn(IntegerMatrix &x, IntegerVector &cnn) {

  // Determine the number of observations
  unsigned int nnn = cnn.size();

  // Create an output matrix
  IntegerMatrix out = no_init(nnn, nnn);

  // Loop through each column and copy the data.
  for(unsigned int i = 0; i < nnn; i++) {
    for (unsigned int j = 0; j < nnn; j++) {
      out(i, j) = x(cnn[i] - 1, cnn[j] - 1);
    }
  }
  return out;
}

// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <progress_bar.hpp>
// [[Rcpp::export]]
std::list<std::list<std::vector<bool>>> scsg2_cpp(List &nn,
                                                  IntegerMatrix &w,
                                                  IntegerVector &idx,
                                                  unsigned int &nlevel,
                                                  NumericVector &lprimes,
                                                  bool verbose=false) {
  // declare int
  unsigned int i, clevel, nidx = idx.length();

  // current nearest neighbors, all region ids
  IntegerVector cnn, s;
  IntegerMatrix cw;

  // current list of zones
  std::list<std::vector<bool>> clz;
  // list of zones for idx i
  std::list<std::list<std::vector<bool>>> zidx;
  // list of zones
  std::list<std::list<std::vector<bool>>> z;

  // display progress
  Progress p(nidx, verbose);

  // loop over desired indices
  for (i = 0; i < nidx; i++) {
    p.increment();
    clevel = 1;
    // get current nn
    cnn = nn[idx[i] - 1];
    // get size of cnn vector
    unsigned int cnn_size = cnn.size();

    // create current w
    cw = sub_cnn(w, cnn);
    // create sequence for cnn
    s = seq(0, cnn_size - 1);

    // create current zone
    std::vector<bool> cz(cnn_size);
    cz[0] = true;

    // clear than create new clz for current zone
    clz.clear();
    clz.push_back(cz);
    // add clz to end of zidx
    zidx.clear();
    zidx.emplace_back(clz);

    // while there are still levels to consider
    while (clevel < nlevel) {
      // increment counter
      clevel++;
      // determine new set of zones from
      // previous list of candidate zones
      // and add to back list of zones
      zidx.emplace_back(lcsg2_cpp(zidx.back(), cnn, cw, s, lprimes));
      //if there are no more connected neighbors for
      // the current vector of zones, end expansion
      if (zidx.back().empty()) {
        clevel = nlevel;
      }
    }
    // clear clz to prepare for storage
    clz.clear();
    // move zidx indices to stacked list of bool vectors
    auto zidx_end = zidx.end();
    for (auto zidx_it = zidx.begin(); zidx_it != zidx_end; zidx_it++) {
      clz.splice(clz.end(), *zidx_it);
    }
    // move complete list of zones for region i to z
    z.emplace_back(std::move(clz));
  }
  return z;
}

// [[Rcpp::export]]
NumericVector stat_poisson_cpp(NumericVector yin,
                               NumericVector yout,
                               NumericVector ein,
                               NumericVector eout,
                               double a,
                               NumericVector shape) {
  unsigned int yin_length = yin.length();
  double lrin, lrout;
  NumericVector tall(yin_length, 0);

  // determine if there will be any problematic statistics
  for (unsigned int i = 0; i < yin_length; i++) {
    // compute statistic for good locations
    // yin > 0 and yin/ein > yout/ein
    if (yin[i] > 0) {
      lrin = log(yin[i]) - log(ein[i]);
      lrout = log(yout[i]) - log(eout[i]);
      if (lrin > lrout) {
        tall[i] = yin[i] * lrin + yout[i] * lrout;
      }
    }
  }

  // update test statistic if a > 0 for all shape values > 1
  if (a > 0) {
    // elliptical regions
    for (unsigned int j = 0; j < yin_length; j++) {
      if (shape[j] > 1) {
        tall[j] = tall[j] * pow(((4 * shape[j]) / pow(shape[j] + 1, 2)), a);
      }
    }
  }
  return tall;
}

// [[Rcpp::export]]
NumericVector stat_poisson0_cpp(NumericVector &yin,
                                NumericVector yout,
                                NumericVector &ein,
                                NumericVector &eout) {
  unsigned int yin_length = yin.length();
  double lrin, lrout;
  NumericVector tall(yin_length, 0);

  // determine if there will be any problematic statistics
  for (unsigned int i = 0; i < yin_length; i++) {
    // compute statistic for good locations
    // yin > 0 and yin/ein > yout/ein
    if (yin[i] > 0) {
      lrin = log(yin[i]) - log(ein[i]);
      lrout = log(yout[i]) - log(eout[i]);
      if (lrin > lrout) {
        tall[i] = yin[i] * lrin + yout[i] * lrout;
      }
    }
  }
  return tall;
}

// [[Rcpp::export]]
NumericVector stat_binom_cpp(NumericVector yin,
                             NumericVector yout,
                             double ty,
                             NumericVector popin,
                             NumericVector popout,
                             double tpop) {

  unsigned int yin_length = yin.length();
  double lrin, lrout, py_in, py_out;
  NumericVector tall(yin_length, 0);

  // determine if there will be any problematic statistics
  for (unsigned int i = 0; i < yin_length; i++) {
    // compute statistic for good locations
    // yin > 0 and yin / popin > yout / popout
    if (yin[i] > 0) {
      lrin = log(yin[i]) - log(popin[i]);
      lrout = log(yout[i]) - log(popout[i]);
      if (lrin > lrout) {
        py_in = popin[i] - yin[i];
        py_out = popout[i] - yout[i];
        tall[i] = yin[i] * lrin +
          py_in * (log(py_in) - log(popin[i])) +
          yout[i] * lrout +
          py_out * (log(py_out) - log(popout[i])) -
          ty * log(ty) - (tpop - ty) * log(tpop - ty) +
          tpop * log(tpop);
        // tall = yin * (log(yin) - log(popin)) +
        //py_in * (log(py_in) - log(popin)) +
        // yout * (log(yout) - log(popout)) +
        // py_out * (log(py_out) - log(popout)) -
        // ty * log(ty) - (tpop - ty) * log(tpop - ty) +
        // tpop * log(tpop)
      }
    }
  }
  return tall;
}

// add elements of nv for which corresponding bv value is true
double sum_nv_bv(NumericVector &nv, std::vector<bool> &bv) {
  double tsum = 0;
  int i = 0;
  for (auto bv_it = bv.begin(); bv_it != bv.end(); bv_it++) {
    if (*bv_it) {
      tsum += nv[i];
    }
    i++;
  }
  return tsum;
}

// Take list of bool vector, add corresponding elements of yidx for each vector
NumericVector zidx_sum(std::list<std::vector<bool>> &zidx,
                       NumericVector &yidx) {
  NumericVector tsums(zidx.size());
  int i = 0;
  auto zidx_end = zidx.end();
  for (auto zidx_it = zidx.begin(); zidx_it != zidx_end; zidx_it++) {
    tsums[i] = sum_nv_bv(yidx, *zidx_it);
    i++;
  }
  return tsums;
}

IntegerVector rmultinom_1(unsigned int &size, NumericVector &probs, unsigned int &N) {
  IntegerVector outcome(N);
  rmultinom(size, probs.begin(), N, outcome.begin());
  return outcome;
}

IntegerMatrix rmultinom_rcpp(unsigned int &n, unsigned int &size, NumericVector &probs) {
  unsigned int N = probs.length();
  IntegerMatrix sim(N, n);
  for (unsigned int i = 0; i < n; i++) {
    sim(_,i) = rmultinom_1(size, probs, N);
  }
  return sim;
}

std::vector<NumericVector> rmultinom_alt(unsigned int &n,
                                         unsigned int &size,
                                         NumericVector &probs) {
  unsigned int N = probs.length();
  std::vector<NumericVector> sim(n);
  for (unsigned int i = 0; i < n; i++) {
    sim[i] = rmultinom_1(size, probs, N);
  }
  return sim;
}

IntegerVector rpois_rcpp(unsigned int &n, NumericVector &lambda) {
  unsigned int lambda_i = 0;
  IntegerVector sim(n);
  for (unsigned int i = 0; i < n; i++) {
    sim[i] = R::rpois(lambda[lambda_i]);
    // update lambda_i to match next realized value with correct mean
    lambda_i++;
    // restart lambda_i at 0 if end of lambda reached
    if (lambda_i == lambda.length()) {
      lambda_i = 0;
    }
  }
  return sim;
}

std::vector<NumericVector> rpois_alt(unsigned int &nsim,
                                     unsigned int &n,
                                     NumericVector &lambda) {
  std::vector<NumericVector> sim(n);
  for (unsigned int i = 0; i < n; i++) {
    sim[i] = rpois_rcpp(n, lambda);
  }
  return sim;
}

NumericVector col_maxes(NumericMatrix &x) {
  unsigned int x_ncol = x.ncol();
  NumericVector maxes(x_ncol);
  for (unsigned int i = 0; i < x_ncol; i++) {
    maxes[i] = max(x(_, i));
  }
  return maxes;
}

// Determine simulated test statistics for flexible scan method
// with poisson statistic and multinomial simulated data
NumericVector flex_sim_prmulti(std::list<std::vector<bool>> &zidx,
                               std::vector<NumericVector> &ysim,
                               IntegerVector &cnn,
                               NumericVector &ein,
                               NumericVector &eout,
                               NumericVector &ty_vec
) {
  unsigned int nsim = ysim.size();
  unsigned int zidx_size = zidx.size();
  NumericVector yidx, yin;
  NumericMatrix tsim_mat(zidx_size, nsim);
  NumericVector tsim;

  for (unsigned int i = 0; i < nsim; i++) {
    yidx = (ysim[i])[cnn - 1];
    yin = zidx_sum(zidx, yidx);
    tsim_mat(_,i) = stat_poisson0_cpp(yin, ty_vec - yin, ein, eout);
  }
  return col_maxes(tsim_mat);
}

// [[Rcpp::export]]
NumericVector mc_pvalue_cpp(NumericVector &tobs, NumericVector &tsim) {
  unsigned int nsim = tsim.length();
  unsigned int ntobs = tobs.length();
  NumericVector pvalues(ntobs);
  for (int i = 0; i < ntobs; i++) {
    pvalues[i] = (static_cast<double>(sum(tsim >= tobs[i])) + 1)/(nsim + 1);
  }
  return pvalues;
}

std::vector<NumericVector> mc_pvalues_nested(std::vector<NumericVector> &tobs_nested,
                                             NumericVector &tsim) {
  std::vector<NumericVector> pvalues_nested;
  unsigned tobs_size = tobs_nested.size();
  for (unsigned int i = 0; i < tobs_size; i++) {
    pvalues_nested.emplace_back(mc_pvalue_cpp(tobs_nested[i], tsim));
  }
  return pvalues_nested;
}

std::vector<NumericVector> pvalues_nested_1(std::vector<NumericVector> &tobs_nested) {
  std::vector<NumericVector> pvalues_nested;
  double one = 1;
  unsigned tobs_size = tobs_nested.size();
  for (unsigned int i = 0; i < tobs_size; i++) {
    NumericVector t1(tobs_nested[i].length(), one);
    pvalues_nested.emplace_back(t1);
  }
  return pvalues_nested;
}


void sig_zidx(std::list<std::vector<bool>> &zidx, LogicalVector &sig) {
  unsigned int sig_length = sig.length();
  auto zidx_it = zidx.begin();
  for(unsigned int i = 0; i < sig_length; i++) {
    if (!sig[i]) {
      zidx.erase(zidx_it++);
    } else {
      zidx_it++;
    }
  }
}

/*
 Take logical vector cz (same length as cnn) and cnn.
 Select elements of cnn corresponding to true
 elements of cz. This corresponds to the current zone, cz
 */
std::vector<int> cz2zone(std::vector<bool> &cz, IntegerVector &cnn) {
  unsigned int cz_size = cz.size();
  std::vector<int> zone;
  for(int i = 0; i < cz_size; i++) {
    if (cz[i]) {
      zone.push_back(cnn[i]);
    }
  }
  return zone;
}

/*
 Take list of logical vectors (list of cz) and cnn
 Turn each element of zidx into a zone (in terms of regions ids)
 */
std::vector<std::vector<int>> logical2zones_nested(std::list<std::vector<bool>> &zidx, IntegerVector &cnn) {
  auto zidx_end = zidx.end();
  std::vector<std::vector<int>> zones_nested;
  for(auto it = zidx.begin(); it != zidx_end; it++) {
    zones_nested.push_back(cz2zone(*it, cnn));
  }
  return zones_nested;
}

List sig_z_tobs_pvalue(std::list<std::list<std::vector<bool>>> &z,
                       std::vector<NumericVector> &tobs_nested,
                       std::vector<NumericVector> &pvalues_nested,
                       double &alpha,
                       List &nn) {
  auto z_end = z.end();
  // std::vector<NumericVector> tobs = tobs_nested;
  // std::vector<NumericVector> pvalues = pvalues_nested;
  unsigned int step = 0;
  LogicalVector sig;
  std::vector<std::vector<int>> zones_nested, zones;
  std::vector<double> tobs, pvalues;
  IntegerVector cnn;
  for (auto z_it = z.begin(); z_it != z_end; z_it++) {
    // determine significant pvalues
    sig = pvalues_nested[step] <= alpha;
    // retain significant pvalues
    pvalues_nested[step] = (pvalues_nested[step])[sig];
    // retain significant tobs
    tobs_nested[step] = (tobs_nested[step])[sig];
    // retain significant zidx
    sig_zidx(*z_it, sig);
    cnn = nn[step];
    zones_nested = logical2zones_nested(*z_it, cnn);
    zones.insert(zones.end(), zones_nested.begin(), zones_nested.end());
    tobs.insert(tobs.end(), tobs_nested[step].begin(), tobs_nested[step].end());
    pvalues.insert(pvalues.end(), pvalues_nested[step].begin(), pvalues_nested[step].end());

    // (*z_it).clear();
    // manually advance iterators for tobs and pvalues
    step++;
  }
  return Rcpp::List::create(Rcpp::Named("zones") = zones,
                            Rcpp::Named("tobs") = tobs,
                            Rcpp::Named("pvalues") = pvalues
  );
}


// [[Rcpp::export]]
List flex_test_cpp (List nn,
                    NumericVector cases,
                    NumericVector pop,
                    IntegerMatrix w,
                    unsigned int k,
                    NumericVector ex,
                    unsigned int type,
                    unsigned int nsim,
                    double alpha,
                    NumericVector lprimes,
                    bool verbose = true) {
  IntegerVector idx = seq(1, nn.length());
  if (verbose) {
    Rcout << "determining candidate zones" << std::endl;
  }
  // number of regions
  unsigned int N = cases.length();
  // total number of cases double and rounded down
  double ty = sum(cases);
  unsigned int uint_ty = floor(sum(cases));

  IntegerVector cnn;
  std::list<std::list<std::vector<bool>>> z = scsg2_cpp(nn, w, idx, k, lprimes, verbose);
  std::vector<NumericVector> tobs_nested, pvalues_nested;
  std::vector<NumericVector> ysim;

  // simulated cases
  if (nsim > 0) {
    NumericVector prob = clone(ex);
    for (unsigned int j = 0; j < N; j++) {
      prob[j] = prob[j]/ty;
    }
    ysim = rmultinom_alt(nsim, uint_ty, prob);
  }

  NumericVector yidx, eidx, yin, yout, ein, eout, tobs_vec;
  NumericMatrix tsim_mat(N, nsim);
  NumericVector tsim(nsim, 0);

  int step = 0;

  if (verbose) {
    Rcout << "computing test statistics" << std::endl;
  }
  // display progress bar
  Progress p_flex(z.size(), verbose);

  for(auto z_it = z.begin(); z_it != z.end(); z_it++) {
    p_flex.increment();
    cnn = nn[step];
    yidx = cases[cnn - 1];
    eidx = ex[cnn - 1];
    yin = zidx_sum(*z_it, yidx);
    NumericVector ty_vec(yin.length(), ty);
    // yout = ty_vec - yin;
    ein = zidx_sum(*z_it, eidx);
    eout = ty_vec - ein;
    tobs_vec = stat_poisson0_cpp(yin, ty_vec - yin, ein, eout);
    tobs_nested.emplace_back(tobs_vec);
    if (nsim > 0) {
      tsim_mat(step,_) = flex_sim_prmulti(*z_it, ysim, cnn, ein, eout, ty_vec);
    }
    step++;
  }
  // If nsim > 0, compute associated pvalues
  // otherwise, all pvalues are 1
  if (nsim > 0) {
    tsim = col_maxes(tsim_mat);
    pvalues_nested = mc_pvalues_nested(tobs_nested, tsim);
  } else {
    pvalues_nested = pvalues_nested_1(tobs_nested);
  }
  return sig_z_tobs_pvalue(z, tobs_nested, pvalues_nested, alpha, nn);
}
