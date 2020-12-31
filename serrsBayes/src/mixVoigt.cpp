// #  -------------------------------------------------------------------------
// #  This file contains C++ code to fit Voigt peaks to spectroscopic data.
// #  For more detail, see:
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
// [[Rcpp::plugins(openmp)]]

Eigen::VectorXd dCauchy(Eigen::VectorXd Cal_V, double loc, double scale)
{
  VectorXd Sigi = VectorXd::Zero(Cal_V.size());
  for (int i=0; i < Cal_V.size(); i++)
  {
    Sigi[i] = 1/(PI*scale*(1 + pow((Cal_V[i] - loc)/scale, 2)));
  }
  return Sigi;
}

Eigen::VectorXd dNorm(Eigen::VectorXd Cal_V, double loc, double sd)
{
  VectorXd Sigi = VectorXd::Zero(Cal_V.size());
  for (int i=0; i < Cal_V.size(); i++)
  {
    Sigi[i] = 1/(sqrt(2*PI)*sd) * exp(-pow(Cal_V[i] - loc, 2)/(2*pow(sd,2)));
  }
  return Sigi;
}

//' Sum log-likelihoods of Gaussian.
//' 
//' This is an internal function that is only exposed on the public API for unit testing purposes.
//' 
//' The sum of the log-likelihoods (log of the product of the likelihoods)
//' for independent, identically-distributed, Gaussian random variables.
//' Note: this Rcpp function is thread-safe, unlike the equivalent alternatives. 
//' 
//' @param x Vector of i.i.d. Gaussian random varibles
//' @param mean Vector of means
//' @param sd Vector of standard deviations
//' @return log-likelihood of x
//' @seealso \code{sum(dnorm(x, mean, sd, log=TRUE))}
//' @examples
//'   x <- rnorm(100)
//'   mu <- rep(0,length(x))
//'   sd <- rep(1,length(x))
//'   sumDnorm(x,mu,sd)
// [[Rcpp::export]]
double sumDnorm(Eigen::VectorXd x, Eigen::VectorXd mean, Eigen::VectorXd sd)
{
  double logLik = 0;
  for (int pk=0; pk < x.size(); pk++)
  {
    logLik += -pow(x[pk] - mean[pk], 2.0)/(2*pow(sd[pk],2.0)) - log(sd[pk] * sqrt(2*PI));
  }
  return logLik;
}

//' Sum log-likelihoods of i.i.d. lognormal.
//' 
//' This is an internal function that is only exposed on the public API for unit testing purposes.
//' 
//' The sum of the log-likelihoods (log of the product of the likelihoods)
//' for independent, identically-distributed, lognormal random variables. 
//' Note: this Rcpp function is thread-safe, unlike the equivalent alternatives. 
//' 
//' @param x Vector of i.i.d. lognormal random varibles
//' @param meanlog mean of the distribution on the log scale
//' @param sdlog standard deviation on the log scale
//' @return log-likelihood of x
//' @seealso \code{sum(dlnorm(x, meanlog, sdlog, log=TRUE))}
//' @examples
//' x <- rlnorm(100)
//' sumDlogNorm(x,0,1)
// [[Rcpp::export]]
double sumDlogNorm(Eigen::VectorXd x, double meanlog, double sdlog)
{
  double var = pow(sdlog, 2.0);
  ArrayXd sqDiff = (x.array().log() - meanlog).pow(2.0);
  ArrayXd logConst = (x.array() * sdlog * sqrt(2*PI)).log();
  ArrayXd logLik = -sqDiff/(2*var) - logConst;
  return logLik.sum();
}

double calcVoigtFWHM(double f_G, double f_L)
{
  // combined scale is the average of the scales of the Gaussian/Lorentzian components
  double Temp_d = pow(f_G,5)
    + 2.69269*pow(f_G,4)*f_L
    + 2.42843*pow(f_G,3)*pow(f_L,2)
    + 4.47163*pow(f_G,2)*pow(f_L,3)
    + 0.07842*f_G*pow(f_L,4)
    + pow(f_L,5);
  return pow(Temp_d, 0.2);
}

//' Compute the spectral signature using Voigt peaks.
//' 
//' Calculates the value of the pseudo-Voigt broadening function at the given wavenumbers,
//' given the parameters of the peaks. This function is thread-safe.
//' 
//' @param location Vector of location parameters of the peaks (\eqn{cm^{-1}})
//' @param scale_G Vector of standard deviations \eqn{\sigma_j} of the Gaussian components.
//' @param scale_L Vector of scale parameters \eqn{\phi_j} of the Lorentzian components.
//' @param amplitude Vector of amplitudes of the peaks (a.u.)
//' @param wavenum Vector of wavenumbers at which to compute the function.
//' @return The value of the pseudo-Voigt function at the given wavenumbers.
//' @examples
//'   Cal_V <- seq(300,400,by=5)
//'   loc <- c(320,350,375)
//'   scG <- c(10,5,1)
//'   scL <- c(3,20,7)
//'   amp <- c(100,500,200)
//'   mixedVoigt(loc,scG,scL,amp,Cal_V)
//' @references
//' Thompson, Cox & Hastings (1987) "Rietveld refinement of Debye--Scherrer synchrotron X-ray data from \eqn{Al_2 O_3},"
//' J. Appl. Crystallogr. 20(2): 79--83, DOI: \href{https://doi.org/10.1107/S0021889887087090}{10.1107/S0021889887087090}
// [[Rcpp::export]]
Eigen::VectorXd mixedVoigt(Eigen::VectorXd location, Eigen::VectorXd scale_G, Eigen::VectorXd scale_L, Eigen::VectorXd amplitude, Eigen::VectorXd wavenum)
{
  VectorXd Sigi = VectorXd::Zero(wavenum.size());
  for (int j=0; j < location.size(); j++)
  {
    // combined scale is the average of the scales of the Gaussian/Lorentzian components
    double f_G = 2.0*scale_G(j)*sqrt(2.0*PI);
    double f_L = 2.0*scale_L(j);
    double Temp_f = calcVoigtFWHM(f_G, f_L);

    // (0,1) Voigt parameter gives the mixing proportions of the two components
    double Temp_e = 1.36603*(f_L/Temp_f) - 0.47719*pow(f_L/Temp_f, 2) + 0.11116*pow(f_L/Temp_f, 3);

    // weighted additive combination of Cauchy and Gaussian functions
    Sigi += amplitude[j] * (Temp_e*dCauchy(wavenum,location[j],Temp_f/2.0) + (1.0-Temp_e)*dNorm(wavenum,location[j],Temp_f/(2.0*sqrt(2.0*log(2.0)))))/(Temp_e*(1.0/(PI*(Temp_f/2.0))) + (1.0-Temp_e)*(1.0/sqrt(2.0*PI*pow(Temp_f/(2.0*sqrt(2.0*log(2.0))), 2.0))));
  }
  return Sigi;
}  

//' Compute the pseudo-Voigt mixing ratio for each peak.
//' 
//' Calculates the mixing parameter \eqn{\eta_j} from the scales of the Gaussian/Lorentzian
//' components.
//' 
//' First, calculate a polynomial average of the scale parameters according to
//' the approximation of Thompson et al. (1987):
//' \deqn{f_{G,L} = (\sigma_j^5 + 2.69\sigma_j^4\phi_j + 2.42\sigma_j^3\phi_j^2 + 4.47\sigma_j^2\phi_j^3 + 0.07\sigma_j\phi_j^4 + \phi_j^5)^{1/5} }
//' 
//' Then the Voigt mixing parameter \eqn{\eta_j} is defined as:
//' \deqn{\eta_j = 1.36\frac{\phi_j}{f_{G,L}} - 0.47(\frac{\phi_j}{f_{G,L}})^2 + 0.11(\frac{\phi_j}{f_{G,L}})^3}
//' 
//' @param scale_G Vector of standard deviations \eqn{\sigma_j} of the Gaussian components.
//' @param scale_L Vector of scale parameters \eqn{\phi_j} of the Lorentzian components.
//' @return The Voigt mixing weights for each peak, between 0 (Gaussian) and 1 (Lorentzian).
//' @references
//' Thompson, Cox & Hastings (1987) "Rietveld refinement of Debye--Scherrer synchrotron X-ray data from \eqn{Al_2 O_3},"
//' J. Appl. Crystallogr. 20(2): 79--83, DOI: \href{https://doi.org/10.1107/S0021889887087090}{10.1107/S0021889887087090}
// [[Rcpp::export]]
Eigen::VectorXd getVoigtParam(Eigen::VectorXd scale_G, Eigen::VectorXd scale_L)
{
  VectorXd voigt = VectorXd::Zero(scale_G.size());
  for (int j=0; j < scale_G.size(); j++)
  {
    // combined scale is the average of the scales of the Gaussian/Lorentzian components
    double f_G = 2.0*scale_G(j)*sqrt(2.0*PI);
    double f_L = 2.0*scale_L(j);
    double Temp_f = calcVoigtFWHM(f_G, f_L);

    // (0,1) Voigt parameter gives the mixing proportions of the two components
    voigt[j] = 1.36603*(f_L/Temp_f) - 0.47719*pow(f_L/Temp_f, 2) + 0.11116*pow(f_L/Temp_f, 3);
  }
  return voigt;
}

//' Initialise the vector of Metropolis-Hastings proposals.
//' 
//' This is an internal function that is only exposed on the public API for unit testing purposes.
//' 
//' @param nPK number of Raman peaks in the spectral signature
//' @param T_Prop_Theta Vector of logarithms of the MH proposals
//' @return Vector of proposals
// [[Rcpp::export]]
Eigen::VectorXd copyLogProposals(int nPK, Eigen::VectorXd T_Prop_Theta)
{
  VectorXd Prop_Theta(4*nPK);
  for (int par = 0; par < 4; par++)
  {
    if (par != 2)
    {
      Prop_Theta.segment(par*nPK,nPK) = T_Prop_Theta.segment(par*nPK,nPK).array().exp();
    }
    else
    {
      Prop_Theta.segment(par*nPK,nPK) = T_Prop_Theta.segment(par*nPK,nPK);
    }
  }
  return Prop_Theta;
}

//' Compute the log-likelihood.
//' 
//' This is an internal function that is only exposed on the public API for unit testing purposes.
//' It computes the log-likelihood of the spline and the noise, once the spectral signature has
//' been subtracted from the observed data. Thus, it can be used with either Lorentzian, Gaussian,
//' or pseudo-Voigt broadening functions.
//' 
//' @param obsi Vector of residuals after the spectral signature has been subtracted.
//' @param lambda smoothing parameter of the penalised B-spline.
//' @param prErrNu hyperparameter of the additive noise
//' @param prErrSS hyperparameter of the additive noise
//' @param basisMx Matrix of B-spline basis functions
//' @param eigVal eigenvalues of the Demmler-Reinsch factorisation
//' @param precMx precision matrix for the spline
//' @param xTx sparse matrix cross-product
//' @param aMx orthoganal matrix A from the Demmler-Reinsch factorisation
//' @param ruMx product of Ru from the Demmler-Reinsch factorisation
//' @return The logarithm of the likelihood.
// [[Rcpp::export]]
double computeLogLikelihood(Eigen::VectorXd obsi, double lambda, double prErrNu, double prErrSS,
            Eigen::MatrixXd basisMx, Eigen::VectorXd eigVal, Eigen::SparseMatrix<double> precMx,
            Eigen::SparseMatrix<double> xTx, Eigen::MatrixXd aMx, Eigen::MatrixXd ruMx)
{
  double nWL = obsi.size();
  double a0_Cal = prErrNu/2.0;
  double ai_Cal = a0_Cal + nWL/2.0;

  SparseMatrix<double> g0_Cal = precMx * nWL * lambda;
  SparseLU < SparseMatrix<double> > g0LU(g0_Cal);
  SparseMatrix<double> gi_Cal = xTx + g0_Cal;
  SparseLU < SparseMatrix<double> > giLU(gi_Cal);
  double g0_Det = g0LU.logAbsDeterminant();
  double gi_Det = giLU.logAbsDeterminant();
//  Rcpp::Rcout << g0_Det << "; " << gi_Det << "; ";
  VectorXd b = aMx.transpose() * obsi;
  if (!b.allFinite()) b = aMx.transpose() * obsi.transpose(); // row vector
//  Rcpp::Rcout << b.allFinite() << "; " << b.mean() << "; ";
  ArrayXd ePlus = eigVal.array() * lambda * nWL + 1.0;
  VectorXd bRatio = (b.array() / ePlus).matrix();
  VectorXd mi_New = ruMx * bRatio;
  double sqDiff = obsi.squaredNorm() - mi_New.transpose() * gi_Cal * mi_New;
//  Rcpp::Rcout << obsi.squaredNorm() << "; " << mi_New.transpose() * gi_Cal * mi_New << "; ";
  double bi_Cal = (prErrSS + sqDiff)/2.0;
//  Rcpp::Rcout << sqDiff << "; " << bi_Cal << "; ";

  // log-likelihood:
  double L_Ev = -((nWL/2.0)*log(2.0*PI)) + 0.5*g0_Det - 0.5*gi_Det;
  L_Ev +=  a0_Cal*log(prErrSS) - ai_Cal*log(bi_Cal) + lgamma(ai_Cal) - lgamma(a0_Cal);
  return L_Ev;
}

//' Update the parameters of the Voigt peaks using marginal Metropolis-Hastings.
//' 
//' Updates all of the parameters (location, amplitude, std. dev., and scale) using a single Metropolis-
//' Hastings step, such that the baseline cancels out in the MH ratio, using the marginalisation identity
//' of Chib (1995).
//' Note: if \code{npart > 1}, then multiple MCMC chains will be executed independently in parallel using
//' OpenMP. This means that all functions used for the proposal distributions and to evaluate the MH ratio
//' need to be thread-safe. Specifically, no calls to \code{R::rnorm}, \code{R::dnorm}, nor their
//' Rcpp equivalents, can be made from within the parallel portion of the code.
//' 
//' @param spectra \code{n_y * nwl} Matrix of observed Raman spectra.
//' @param n number of observations to use in calculating the likelihood.
//' @param kappa likelihood tempering parameter.
//' @param conc Vector of \code{n_y} nanomolar (nM) dye concentrations
//' @param wavenum Vector of \code{nwl} wavenumbers at which the spetra are observed.
//' @param thetaMx \code{(4+npeaks*4) x npart} Matrix of parameter values for each peak.
//' @param logThetaMx \code{(4+npeaks*4) x npart} Matrix of logarithms of the parameters.
//' @param mhChol lower-triangular Cholesky factorisation of the covariance matrix for the random walk proposals.
//' @param priors List of hyperparameters for the prior distributions.
//' @return The number of RWMH proposals that were accepted.
//' @references
//' Chib (1995) "Marginal Likelihood from the Gibbs Output," JASA 90(432): 1313--1321,
//' DOI: \href{http://dx.doi.org/10.1080/01621459.1995.10476635}{10.1080/01621459.1995.10476635}
//' 
//' Rosenthal (2000) "Parallel computing and Monte Carlo algorithms" Far East J. Theor. Stat. 4(2): 207--236,
//' URL: \href{http://www.pphmj.com/abstract/1961.htm}{http://www.pphmj.com/abstract/1961.htm}
// [[Rcpp::export]]
long mhUpdateVoigt(Eigen::MatrixXd spectra, unsigned n, double kappa, Eigen::VectorXd conc, Eigen::VectorXd wavenum,
                   NumericMatrix thetaMx, NumericMatrix logThetaMx, Eigen::MatrixXd mhChol, List priors)
{
  // priors
  double prErrNu = priors["noise.nu"];
  double prErrSS = priors["noise.SS"];
  double prScaGmu = priors["scaG.mu"]; // squared exponential (Gaussian) peaks
  double prScaGsd = priors["scaG.sd"];
  double prScaLmu = priors["scaL.mu"]; // Lorentzian (Cauchy) peaks
  double prScaLsd = priors["scaL.sd"];
  VectorXd prLocMu = priors["loc.mu"];
  VectorXd prLocSD = priors["loc.sd"];
  VectorXd prAmpMu, prAmpSD;
  if (priors.containsElementNamed("beta.mu"))
  {
    prAmpMu = priors["beta.mu"];
    prAmpSD = priors["beta.sd"];
  }
  double lambda = priors["bl.smooth"];
  int nPK = prLocMu.size();
  int nPart = thetaMx.rows();
  int nWL = wavenum.size();

  // matrices for the cubic B-spline
  MatrixXd basisMx = priors["bl.basis"];
  VectorXd eigVal = priors["bl.eigen"];
  const MappedSparseMatrix<double> precMx(as<MappedSparseMatrix<double> >(priors["bl.precision"]));
  const MappedSparseMatrix<double> xTx(as<MappedSparseMatrix<double> >(priors["bl.XtX"]));
  MatrixXd aMx = priors["bl.orthog"]; // orthogonal, Demmler-Reinsch basis
  MatrixXd ruMx = priors["bl.Ru"];

  // RNG is not thread-safe
  const NumericVector stdNorm = rnorm(nPK * nPart * 4, 0, 1);
  const NumericVector rUnif = runif(nPart, 0, 1);

  long accept = 0;
#pragma omp parallel for default(shared) reduction(+:accept)
  for (int pt = 0; pt < nPart; pt++)
  {
    VectorXd theta(nPK*4), logTheta(nPK*4), stdVec(nPK*4);
    for (int pk = 0; pk < nPK*4; pk++)
    {
      stdVec(pk) = stdNorm[pt*nPK*4 + pk];
      theta(pk) = thetaMx(pt,pk);
      logTheta(pk) = logThetaMx(pt,pk);
    }
    VectorXd T_Prop_Theta = mhChol * stdVec + logTheta;
    VectorXd Prop_Theta = copyLogProposals(nPK, T_Prop_Theta);
    // enforce the boundary condition for proposed peak locations
    for (int pk = 0; pk < nPK; pk++)
    {
      if (Prop_Theta(2*nPK+pk) < wavenum(0) || Prop_Theta(2*nPK+pk) > wavenum(nWL-1))
      {
        Prop_Theta(2*nPK+pk) = theta(2*nPK+pk);
      }
    }
    std::sort(Prop_Theta.data() + 2*nPK,
              Prop_Theta.data() + 3*nPK - 1); // for identifiability

    VectorXd sigi = conc(n-1) * mixedVoigt(Prop_Theta.segment(2*nPK,nPK), Prop_Theta.segment(0,nPK),
       Prop_Theta.segment(nPK,nPK), Prop_Theta.segment(3*nPK,nPK), wavenum);
    VectorXd obsi = spectra.row(n-1).transpose() - sigi;

    // smoothing spline:
    //double lambda = thetaMx(pt,4*nPK+2) / thetaMx(pt,4*nPK+3);
    // log-likelihood:
    double L_Ev = computeLogLikelihood(obsi, lambda, prErrNu, prErrSS, basisMx, eigVal,
                                       precMx, xTx, aMx, ruMx);
    double lLik = kappa*L_Ev + sumDlogNorm(Prop_Theta.segment(0,nPK), prScaGmu, prScaGsd);
    lLik += sumDlogNorm(Prop_Theta.segment(nPK,nPK), prScaLmu, prScaLsd);
    lLik += sumDnorm(Prop_Theta.segment(2*nPK,nPK), prLocMu, prLocSD);
    lLik += -kappa*thetaMx(pt,4*nPK+n) - sumDlogNorm(theta.segment(0,nPK), prScaGmu, prScaGsd);
    lLik -= sumDlogNorm(theta.segment(nPK,nPK), prScaLmu, prScaLsd);
    lLik -= sumDnorm(theta.segment(2*nPK,nPK), prLocMu, prLocSD);
    if (priors.containsElementNamed("beta.mu"))
    {
      lLik += sumDnorm(Prop_Theta.segment(3*nPK,nPK), prAmpMu, prAmpSD);
      lLik -= sumDnorm(theta.segment(3*nPK,nPK), prAmpMu, prAmpSD);
    }

    // account for previous observations when n > 1
    VectorXd oldLogLik(n);
    for (int i=0; i < n-1; i++) {
      sigi = conc(i) * mixedVoigt(Prop_Theta.segment(2*nPK,nPK), Prop_Theta.segment(0,nPK),
                  Prop_Theta.segment(nPK,nPK), Prop_Theta.segment(3*nPK,nPK), wavenum);
      obsi = spectra.row(i).transpose() - sigi;
      oldLogLik(i) = computeLogLikelihood(obsi, lambda, prErrNu, prErrSS, basisMx, eigVal,
                precMx, xTx, aMx, ruMx);
      lLik += oldLogLik(i);
      lLik -= thetaMx(pt,4*nPK+i+1);
    }
    oldLogLik(n-1) = L_Ev;

    if (std::isfinite(lLik) && log(rUnif[pt]) < lLik)
    {
      for (int pk=0; pk < nPK*4; pk++)
      {
         logThetaMx(pt,pk) = T_Prop_Theta(pk);
         thetaMx(pt,pk) = Prop_Theta(pk);
      }
      for (int i=0; i < n; i++) {
        logThetaMx(pt,4*nPK+i+1) = oldLogLik(i);
        thetaMx(pt,4*nPK+i+1) = oldLogLik(i);
      }
      accept += 1;
    }
  }
  return accept;
}
