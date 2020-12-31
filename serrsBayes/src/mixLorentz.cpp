// #  -------------------------------------------------------------------------
// #  This file contains C++ code to fit Gaussian or Lorentzian peaks to
// #  spectroscopic data. For more detail, see:
// #  Moores; Gracie; Carson; Faulds; Graham & Girolami (2016; v2 2018)
// #  "Bayesian modelling and quantification of Raman spectroscopy"
// #  https://arxiv.org/abs/1604.07299
// #
// #  Copyright (C) 2017,2018  University of Warwick
// # 
// #  This source code is free software: you can redistribute it and/or modify
// #  it under the terms of the GNU General Public License as published by
// #  the Free Software Foundation, either version 2 of the License, or
// #  (at your option) any later version.
// # 
// #  This software is distributed in the hope that it will be useful,
// #  but WITHOUT ANY WARRANTY; without even the implied warranty of
// #  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// #  GNU General Public License for more details.
// # 
// #  You should have received a copy of the GNU General Public License
// #  along with this program.  If not, see <http://www.gnu.org/licenses/>.
// #  -------------------------------------------------------------------------
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
using namespace Rcpp;
using namespace Eigen;
#include <sys/time.h>
// [[Rcpp::plugins(openmp)]]

//' Compute the spectral signature using Lorentzian peaks.
//' 
//' Calculates the value of the Lorentzian function at the given wavelengths,
//' given the parameters of the peaks. This function is thread-safe.
//' 
//' @param location Vector of location parameters of the peaks.
//' @param scale Vector of scale parameters of the peaks.
//' @param amplitude Vector of amplitudes of the peaks.
//' @param wavelengths Vector of wavenumbers at which to compute the function.
//' @return The value of the Lorentian function at the given wavelengths.
//' @examples
//'   Cal_V <- seq(300,400,by=5)
//'   loc <- c(320,350,375)
//'   sca <- c(10,5,18)
//'   amp <- c(1000,5000,2000)
//'   weightedLorentzian(loc,sca,amp,Cal_V)
// [[Rcpp::export]]
Eigen::VectorXd weightedLorentzian(Eigen::VectorXd location, Eigen::VectorXd scale, Eigen::VectorXd amplitude, Eigen::VectorXd wavelengths)
{
  VectorXd y = VectorXd::Zero(wavelengths.size());
  for (int i=0; i < wavelengths.size(); i++)
  {
    for (int j=0; j < location.size(); j++)
    {
      y[i] += amplitude[j] * pow(scale[j], 2)/(pow(wavelengths[i] - location[j], 2) + pow(scale[j], 2));
    }
  }
  return y;
}

// the log sum of a vector of logs
// http://jblevins.org/log/log-sum-exp
double sum_logs(NumericVector log_prob)
{
  double suml = 0.0;
  NumericVector::iterator it = std::max_element(log_prob.begin(), log_prob.end());
  double maxl = *it;
  for (int i=0; i < log_prob.length(); i++)
  {
    if (std::isfinite(log_prob(i)))
      suml += exp(log_prob(i) - maxl);
  }
  return log(suml) + maxl;
}

//' Compute the effective sample size (ESS) of the particles.
//' 
//' The ESS is a "rule of thumb" for assessing the degeneracy of
//' the importance distribution:
//' \deqn{ESS = \frac{(\sum_{q=1}^Q w_q)^2}{\sum_{q=1}^Q w_q^2}}
//' 
//' @param log_weights logarithms of the importance weights of each particle.
//' @return the effective sample size, a scalar between 0 and Q
//' @examples
//' x <- runif(100)
//' effectiveSampleSize(log(x))
//' @references
//' Liu, JS (2001) "Monte Carlo Strategies in Scientific Computing." Springer, NY, pp. 34--36.
// [[Rcpp::export]]
double effectiveSampleSize(NumericVector log_weights)
{
  double sum_wt = sum_logs(log_weights);
  double sum_sq = sum_logs(log_weights + log_weights);
  double res = exp(sum_wt + sum_wt - sum_sq);
  if (std::isfinite(res)) return res;
  else return 0;
}

// internal utility function that saves duplication of error-prone code
inline int getIdx3D(const int i, const int j, const int k, const int I, const int J)
{
  return i + j*I + k*I*J;
}

inline double diffTM(const struct timeval &x, const struct timeval &y)
{
  double xs = x.tv_sec;
  xs += x.tv_usec / 1.0E6;
  double xy = y.tv_sec;
  xy += y.tv_usec / 1.0E6;
  return xs - xy;
}


//' Update the importance weights of each particle.
//' 
//' @param spectra \code{n_y * nwl} Matrix of observed Raman spectra.
//' @param peaks \code{nwl * npart} Matrix containing the spectral signatures for each observation.
//' @param baselines \code{nwl * npart} Matrix containing the current values of the baselines.
//' @param i index of the current observation to use in calculating the likelihood
//' @param start index of the next wavelength to use in calculating the likelihood, permuted by \code{idx}
//' @param sigma Vector of \code{npart} standard deviations for each particle.
//' @param old_weights logarithms of the importance weights of each particle.
//' @param alpha the target learning rate for the reduction in effective sample size (ESS).
//' @param idx permutation of the indices of the wavelengths.
//' @return a List containing:
//' \describe{
//'   \item{\code{ess}}{The effective sample size, after reweighting.}
//'   \item{\code{weights}}{Vector of updated importance weights.}
//'   \item{\code{index}}{index of the last wavelength used.}
//'   \item{\code{evidence}}{SMC estimate of the logarithm of the model evidence.}
//' }
//' @references
//' Pitt, dos Santos Silva, Giordani & Kohn (2012)
//' "On some properties of Markov chain Monte Carlo simulation methods based on the particle filter"
//' J. Econometrics 171(2): 134--151,
//' DOI: \href{http://dx.doi.org/10.1016/j.jeconom.2012.06.004}{10.1016/j.jeconom.2012.06.004}
//' 
//' Zhou, Johansen & Aston (2015) "Towards Automatic Model Comparison: An Adaptive Sequential Monte Carlo Approach"
//' \href{http://arxiv.org/abs/1303.3123}{arXiv:1303.3123} [stat.ME]
// [[Rcpp::export]]
List reWeightParticles(NumericMatrix spectra, NumericMatrix peaks, NumericMatrix baselines, int i, int start,
                       NumericVector sigma, NumericVector old_weights, double alpha, IntegerVector idx)
{
  int npart = peaks.ncol();
  int nwl = peaks.nrow();
  int n_y = spectra.nrow();
  NumericVector logWt(old_weights); // copy ctor
  double ess, oldESS = effectiveSampleSize(old_weights);
  double oldSumWt = sum_logs(old_weights); // equals zero if weights were already normalised
  Rcpp::Rcout << "previous ESS " << oldESS << " (target: " << alpha*oldESS << " for observation ";
  Rcpp::Rcout << i << " of " << n_y << "; wavenumber " << start << " of " << nwl << ")\n";
  
  int j=start-1;
  for (; j < nwl; j++)
  {
    for (int q=0; q < npart; q++)
    {
      logWt[q] += R::dnorm(spectra(i-1,idx(j)-1), baselines(idx(j)-1,q) + peaks(idx(j)-1,q), sigma[q], 1);
    }
    ess = effectiveSampleSize(logWt);
    if (ess < alpha*oldESS)
    {
      Rcpp::Rcout << "Required " << j-start+2 << " iterations to reduce ESS from " << oldESS << " to " << ess << "\n";
      break;
    }
  }
  // normalize the weights
  double sum_wt = sum_logs(logWt);
  return Rcpp::List::create(
    Rcpp::Named("ess")      = ess,             // effective sample size
    Rcpp::Named("weights")  = logWt - sum_wt,  // importance weights
    Rcpp::Named("index")    = j+1,
    Rcpp::Named("evidence") = sum_wt - oldSumWt
  );
}

//' Compute an ancestry vector for residual resampling of the SMC particles.
//' 
//' @param log_wt logarithms of the importance weights of each particle.
//' @return Vector of indices to the particles that will be propagated forward to the next generation (i.e. the parents)
//' @references
//' Liu & Chen (1998) "Sequential Monte Carlo methods for dynamic systems," JASA 93(443): 1032-1044,
//' DOI: \href{http://dx.doi.org/10.1080/01621459.1998.10473765}{10.1080/01621459.1998.10473765}
//' 
//' Douc, Cappe & Moulines (2005) "Comparison of resampling schemes for particle filtering"
//' In Proc. 4th IEEE Int. Symp. ISPA, pp. 64-69,
//' DOI: \href{http://dx.doi.org/10.1109/ISPA.2005.195385}{10.1109/ISPA.2005.195385}
// [[Rcpp::export]]
Eigen::ArrayXi residualResampling(NumericVector log_wt)
{
  const int n = log_wt.size();
  ArrayXi idx(n);

  // first loop is deterministic: only accept particles with n*weight > 1
  int r=0;
  for (int i=0; i<n; i++)
  {
    if (std::isfinite(log_wt(i)))
    {
      int tW = (int)trunc(exp(log_wt(i) + log((double)n)));
      for (int j=0; j < tW; j++)
      {
        idx[r+j] = i;
      }
      r += tW;
      log_wt(i) = log(exp(log_wt(i) + log((double)n)) - (double)tW);
    }
  }
  // renormalize the weights
  log_wt = log_wt - log(((double)n-r));

  // second loop uses multinomial resampling for the remaining n-r particles
  const NumericVector randU = runif(n-r);
  for (int i=r; i<n; i++)
  {
    // select a particle at random, according to the weights
    double total = 0.0;
    for (int j=0; j < n && total <= randU[i-r]; j++)
    {
      if (std::isfinite(log_wt(j)))
      {
        total += exp(log_wt(j));
      }
      idx[i] = j;
    }
  }

  // permute the index vector to ensure Condition 9 of Murray, Lee & Jacob (2015)
  for (int i=0; i<n; i++)
  {
    if ((idx[i] != i) && (idx[idx[i]] != idx[i]))
    {
      int old = idx[i];
      idx[i] = idx[idx[i]];
      idx[idx[i]] = old;
    }
  }
  return idx + 1;
}

//' Resample in place to avoid expensive copying of data structures, using a permutation
//' of the ancestry vector.
//' 
//' @param log_weights logarithms of the importance weights of each particle
//' @param ampMx \code{npeaks x npart} Matrix of amplitudes for each particle.
//' @param scaleMx \code{npeaks x npart} Matrix of scale parameters for each particle.
//' @param peaks \code{nwl x npart} Matrix containing the expectation of the Lorentzian mixture.
//' @param baselines \code{nwl x n_y x npart} Array of smoothing splines.
//' @param n_y number of observations
//' @param nwl number of wavenumbers
//' @return Vector of indices to the parents of the resampled particles.
//' @references
//' Murray, L.M., Lee, A. & Jacob, P.E. (2015) "Parallel resampling in the particle filter" \href{http://arxiv.org/abs/1301.4019}{arXiv:1301.4019v3}
//' @seealso \code{\link{residualResampling}}
// [[Rcpp::export]]
Eigen::ArrayXi resampleParticles(NumericVector log_weights, NumericMatrix ampMx, NumericMatrix scaleMx,
                                NumericMatrix peaks, NumericVector baselines, int n_y, int nwl)
{
  struct timeval t1,t2;
  ArrayXi idx = residualResampling(log_weights);

#pragma omp parallel for
  for (int p=0; p < idx.size(); p++)
  {
    // do nothing unless the particle has no offspring
    if (idx[p] != p+1)
    {
      for (int j=0; j < ampMx.rows(); j++)
      {
        ampMx(j,p) = ampMx(j,idx[p]-1);
        scaleMx(j,p) = scaleMx(j,idx[p]-1);
        peaks(j,p) = peaks(j,idx[p]-1);
      }
      for (int i=0; i < n_y; i++)
      {
        for (int w=0; w < nwl; w++)
        {
          int newIdx = getIdx3D(w,i,p,nwl,n_y);
          int oldIdx = getIdx3D(w,i,idx[p]-1,nwl,n_y);
          baselines[newIdx] = baselines[oldIdx];
        }
      }
    }
  }
  return idx;
}

//' Compute the weighted arithmetic means of the particles.
//' 
//' This SMC estimate of the means can be used to centre independent Metropolis-Hastings proposals.
//' 
//' @param particles \code{npeaks * npart} Matrix of parameter values for each particle.
//' @param log_weights logarithms of the importance weights of each particle.
//' @return A vector of means, one for each row.
//' @seealso \code{\link[stats:weighted.mean]{weighted.mean}}
// [[Rcpp::export]]
NumericVector weightedMean(NumericMatrix particles, NumericVector log_weights)
{
  NumericVector suml(particles.nrow());
  NumericVector::iterator it = std::max_element(log_weights.begin(), log_weights.end());
  double maxl = *it;
  for (int i=0; i < log_weights.size(); i++)
  {
    if (std::isfinite(log_weights[i]))
    {
      for (int j=0; j < particles.nrow(); j++)
      {
        suml[j] += exp(log_weights(i) - maxl) * particles(j,i);
      }      
    }
  }
  return suml * exp(maxl);
}

//' Compute the weighted variance of the particles.
//' 
//' This SMC estimate of the variance can be used to scale the bandwidth of adaptive,
//' Gaussian random walk Metropolis-Hastings proposals.
//' 
//' @param particles \code{npeaks * npart} Matrix of parameter values for each particle.
//' @param log_weights logarithms of the importance weights of each particle.
//' @param mean Vector of weighted means of each particle.
//' @return A vector of variances, one for each row.
//' @seealso \code{\link[Hmisc:wtd.stats]{wtd.var}}
// [[Rcpp::export]]
NumericVector weightedVariance(NumericMatrix particles, NumericVector log_weights, NumericVector mean)
{
  NumericVector suml(particles.nrow());
  NumericVector::iterator it = std::max_element(log_weights.begin(), log_weights.end());
  double maxl = *it;
  for (int i=0; i < log_weights.size(); i++)
  {
    if (std::isfinite(log_weights[i]))
    {
      for (int j=0; j < particles.nrow(); j++)
      {
        suml[j] += exp(log_weights(i) - maxl) * pow(particles(j,i) - mean[j], 2.0);
      }      
    }
  }
  return suml * exp(maxl);
}

//' Compute the spectral signature using Gaussian peaks.
//' 
//' Calculates the value of the squared exponential radial basis function at the given wavelengths,
//' given the parameters of the peaks. This function is thread-safe.
//' 
//' @param location Vector of location parameters of the peaks (mean).
//' @param scale Vector of scale parameters of the peaks (standard deviation).
//' @param amplitude Vector of amplitudes of the peaks.
//' @param wavelengths Vector of wavenumbers at which to compute the function.
//' @return The value of the Gaussian function at the given wavelengths.
//' @examples
//'   Cal_V <- seq(300,400,by=5)
//'   loc <- c(320,350,375)
//'   sca <- c(10,5,18)
//'   amp <- c(1000,5000,2000)
//'   weightedGaussian(loc,sca,amp,Cal_V)
// [[Rcpp::export]]
Eigen::VectorXd weightedGaussian(Eigen::VectorXd location, Eigen::VectorXd scale, Eigen::VectorXd amplitude, Eigen::VectorXd wavelengths)
{
  VectorXd y = VectorXd::Zero(wavelengths.size());
  for (int i=0; i < wavelengths.size(); i++)
  {
    for (int j=0; j < location.size(); j++)
    {
      y[i] += amplitude[j] * exp(-pow(wavelengths[i] - location[j], 2)/(2*pow(scale[j],2)));
    }
  }
  return y;
}

//' Update all of the parameters using a single Metropolis-Hastings step.
//' 
//' Updates all of the parameters using a single Metropolis-Hastings step, such that the
//' baseline cancels out in the MH ratio, using the marginalisation identity of Chib (1995).
//' If \code{npart > 1}, then multiple MCMC chains will be executed independently in parallel using OpenMP.
//' This means that all functions used for the proposal distributions and to evaluate the MH ratio
//' need to be thread-safe. Specifically, no calls to \code{R::rnorm}, \code{R::dnorm}, nor their
//' Rcpp equivalents, can be made from within the parallel portion of the code.
//' 
//' @param spectra \code{n_y * nwl} Matrix of observed Raman spectra.
//' @param n number of observations to use in calculating the likelihood
//' @param conc Vector of \code{n} nanomolar (nM) dye concentrations
//' @param wavelengths Vector of \code{nwl} wavenumbers at which the spetra are observed.
//' @param peakWL Vector of locations for each peak (cm^-1)
//' @param betaMx \code{npeaks * npart} Matrix of regression coefficients to update.
//' @param scaleMx \code{npeaks * npart} Matrix of scale parameters to update.
//' @param sigma Vector of \code{npart} standard deviations to update.
//' @param expMx \code{nwl * npart} Matrix of expectations of the Lorentzian or Gaussian function.
//' @param baselines \code{nKnots * n_y * npart} Array of smoothing splines.
//' @param sd_mh Vector of \code{2 * npeaks} bandwidths for the random walk proposals.
//' @param priors List of hyperparameters for the prior distributions.
//' @return The number of RWMH proposals that were accepted.
//' @references
//' Chib (1995) "Marginal Likelihood from the Gibbs Output," JASA 90(432): 1313--1321,
//' DOI: \href{http://dx.doi.org/10.1080/01621459.1995.10476635}{10.1080/01621459.1995.10476635}
//' 
//' Rosenthal (2000) "Parallel computing and Monte Carlo algorithms" Far East J. Theor. Stat. 4(2): 207--236,
//' URL: \href{http://www.pphmj.com/abstract/1961.htm}{http://www.pphmj.com/abstract/1961.htm}
// [[Rcpp::export]]
long marginalMetropolisUpdate(Eigen::MatrixXd spectra, unsigned n, Eigen::VectorXd conc, Eigen::VectorXd wavelengths,
              Eigen::VectorXd peakWL, NumericMatrix betaMx, NumericMatrix scaleMx, NumericVector sigma,
              NumericMatrix expMx, NumericVector baselines, Eigen::VectorXd sd_mh, List priors)
{
  // priors
  double priorNoiseNu = priors["noise.nu"];
  double priorNoiseSS = priors["noise.SS"];
  double priorScaleMu = priors["scale.mu"];
  double priorScaleSD = priors["scale.sd"];
  double priorBetaMu = priors["beta.mu"];
  double priorBetaSD = priors["beta.sd"];
  bool lorentzians = true;
  if (priors.containsElementNamed("peaks")) // && !List::is_na(priors["peaks"])) ISNA segfaults!
  {
    std::string szPeakType = Rcpp::as<std::string>(priors["peaks"]);
    lorentzians = (szPeakType.compare(0, szPeakType.length(), "Lorentzian", 0, szPeakType.length()) == 0);
  }
  int nWL = wavelengths.size();
  int nPart = betaMx.cols();
  MatrixXd basisMx = priors["bl.basis"];
  const MappedSparseMatrix<double> prPrecMx(as<MappedSparseMatrix<double> >(priors["bl.precision"]));
  int nBasis = prPrecMx.cols();
  SparseMatrix<double> xTx = (basisMx.transpose() * basisMx).sparseView();
  SparseMatrix<double> giPrecMx = xTx + prPrecMx;
  SimplicialLLT< SparseMatrix<double> > giChol(giPrecMx);
  MatrixXd giCovMx = giChol.solve(MatrixXd::Identity(nBasis, nBasis));
  MatrixXd smoother = giCovMx * basisMx.transpose();
  double nobs = n*nWL; // cast from int to double for log-likelihood

  // RNG is not thread-safe
  const NumericVector betaVec = rnorm(betaMx.rows() * betaMx.cols(), 0, 1);
  const NumericVector scaleVec = rnorm(scaleMx.rows() * scaleMx.cols(), 0, 1);
  const NumericVector blVec = rnorm(nBasis * nPart * n, 0, 1);
  const NumericVector sdVec = rgamma(nPart, (priorNoiseNu + nobs)/2.0); // scale = 1
  const NumericVector unifVec = runif(nPart);

  long accept = 0;
#pragma omp parallel for default(shared) reduction(+:accept)
  for (int pt = 0; pt < nPart; pt++)
  {
    // independent RWMH proposals for the parameters of the peaks
    VectorXd beta_prop(betaMx.rows());
    VectorXd scale_prop(scaleMx.rows());
    VectorXd mu_prop(nWL);
    bool negProp = false;
    for (int j=0; j < betaMx.rows(); j++)
    {
      beta_prop[j] = sd_mh[j]*betaVec(j + pt*betaMx.rows()) + betaMx(j,pt);
      if (beta_prop[j] < 0)
      {
        negProp = true;
      }
      scale_prop[j] = sd_mh[betaMx.rows() + j]*scaleVec(j + pt*scaleMx.rows()) + scaleMx(j,pt);
      if (scale_prop[j] < 0) negProp = true;
    }

    if (!negProp)
    {
      // translate the parameters into a continuous function, evaluated at the wavenumbers
      if (lorentzians)
      {
        mu_prop = weightedLorentzian(peakWL, scale_prop, beta_prop, wavelengths);
      }
      else
      {
        mu_prop = weightedGaussian(peakWL, scale_prop, beta_prop, wavelengths);
      }
      
      // subtract the peaks from the observed spectra
      MatrixXd diffMx(nWL, n), oldDiffMx(nWL, n);
      for (int i=0; i<n; i++)
      {
        for (int w=0; w<nWL; w++)
        {
          diffMx(w,i) = spectra(i,w) - conc(i)*mu_prop(w);
          oldDiffMx(w,i) = spectra(i,w) - conc(i)*expMx(w,pt);
        }
      }
      
      // Gibbs proposal for the standard deviation of the noise
      MatrixXd meanMx = smoother * diffMx;
      MatrixXd oldMux = smoother * oldDiffMx;
      double sqDiff = (diffMx.transpose() * diffMx).diagonal().sum() - (meanMx.transpose() * giPrecMx * meanMx).diagonal().sum();
      double oldSqDiff = (oldDiffMx.transpose() * oldDiffMx).diagonal().sum() - (oldMux.transpose() * giPrecMx * oldMux).diagonal().sum();
      double newSS = (priorNoiseSS + sqDiff)/2.0;
      double oldSS = (priorNoiseSS + oldSqDiff)/2.0;
      double newTau = sdVec[pt] / newSS;
      double sigma_prop = 1/sqrt(newTau);
      
      // Gibbs proposal for the parameters of the smoothing spline
      MatrixXd stdNorm(nBasis, n), oldAlpha(nBasis, n);
      for (int i=0; i<n; i++)
      {
        for (int b=0; b<nBasis; b++)
        {
          int randIdx = getIdx3D(b, i, pt, nBasis, n);
          stdNorm(b,i) = blVec[randIdx];
          int blIdx = getIdx3D(b, i, pt, nBasis, spectra.rows());
          oldAlpha(b,i) = baselines[blIdx];
        }
      }
      SimplicialLLT< SparseMatrix<double> > prChol(prPrecMx * newTau);
      SimplicialLLT< SparseMatrix<double> > blChol(giPrecMx * newTau);
      MatrixXd y = blChol.permutationP() * blChol.matrixL().solve(stdNorm);
      MatrixXd newAlpha = y + meanMx;
      MatrixXd blEst = basisMx * newAlpha;
      //TODO: reject if any blEst < 0
      MatrixXd blOld = basisMx * oldAlpha;
      MatrixXd resid = diffMx - blEst;
      sqDiff = resid.array().square().sum();
      oldSqDiff = (oldDiffMx - blOld).array().square().sum();
      
      // M-H ratio
      double u, logR = 0.0;
      for (int j=0; j < betaMx.rows(); j++)
      {
        logR += - pow(beta_prop[j] - priorBetaMu, 2)/(2*pow(priorBetaSD, 2));
        logR -= - pow(betaMx(j,pt) - priorBetaMu, 2)/(2*pow(priorBetaSD, 2));
        // lognormal prior for scale
        logR += -log(sqrt(2*M_PI)*priorScaleSD*scale_prop[j]) - pow(log(scale_prop[j]) - priorScaleMu, 2)/(2*pow(priorScaleSD, 2));
        logR -= -log(sqrt(2*M_PI)*priorScaleSD*scaleMx(j,pt)) - pow(log(scaleMx(j,pt)) - priorScaleMu, 2)/(2*pow(priorScaleSD, 2));
      }
      
      // inverse-gamma for std. dev.
      logR += (priorNoiseNu/2 - 1)*log(newTau) - newTau*priorNoiseSS/2;
      logR -= - (priorNoiseNu - 2)*log(sigma[pt]) - priorNoiseSS/(2*pow(sigma[pt], 2.0));
      logR += (priorNoiseNu + nobs)*log(oldSS)/2  - (priorNoiseNu + nobs - 2)*log(sigma[pt]) - oldSS/pow(sigma[pt], 2.0);
      logR -= (priorNoiseNu + nobs)*log(newSS)/2 + ((priorNoiseNu + nobs)/2 - 1)*log(newTau) - newSS*newTau;
      
      // multivariate normal for baselines
      SimplicialLLT< SparseMatrix<double> > oldPRchol(prPrecMx * pow(sigma[pt], -2.0));
      SimplicialLLT< SparseMatrix<double> > oldBLchol(giPrecMx * pow(sigma[pt], -2.0));
      MatrixXd yNew = prChol.matrixL().transpose().toDense() * prChol.permutationP() * newAlpha;
      MatrixXd yOld = oldPRchol.matrixL().transpose().toDense() * oldPRchol.permutationP() * oldAlpha;
      MatrixXd zNew = blChol.matrixL().transpose().toDense() * blChol.permutationP() * (newAlpha - meanMx);
      MatrixXd zOld = oldBLchol.matrixL().transpose().toDense() * oldBLchol.permutationP() * (oldAlpha - oldMux);
      logR += prChol.matrixL().toDense().diagonal().array().log().sum() -0.5 * yNew.array().square().sum();
      logR -= oldPRchol.matrixL().toDense().diagonal().array().log().sum() -0.5 * yOld.array().square().sum();
      logR -= blChol.matrixL().toDense().diagonal().array().log().sum() -0.5 * zNew.array().square().sum();
      logR += oldBLchol.matrixL().toDense().diagonal().array().log().sum() -0.5 * zOld.array().square().sum();
      
      // calculate log-likelihood of the observed spectrum
      logR += -nobs*log(sigma_prop) - 0.5*sqDiff*newTau;
      logR -= -nobs*log(sigma[pt]) - 0.5*oldSqDiff/pow(sigma[pt], 2.0);
      
      // accept/reject
      u = unifVec[pt];
      if (log(u) < logR)
      {
        // accept the joint proposal
        for (int j=0; j < betaMx.rows(); j++) {
          betaMx(j, pt) = beta_prop[j];
          scaleMx(j, pt) = scale_prop[j];
        }
        for (int w=0; w<nWL; w++) {
          expMx(w, pt) = mu_prop[w];
        }
        sigma(pt) = sigma_prop;
        for (int i=0; i<n; i++)
        {
          for (int b=0; b<nBasis; b++)
          {
            int blIdx = getIdx3D(b, i, pt, nBasis, spectra.rows());
            baselines[blIdx] = newAlpha(b, i);
          }
        }
        accept += 1;
      }
    }
  }
  return accept;
}

