\name{graphicalVAR}
\alias{graphicalVAR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimate the graphical VAR model.
}
\description{
Estimates the graphical VAR (Wild et al., 2010) model through LASSO estimation coupled with extended Bayesian information criterion for choosing the optimal tuning parameters. The estimation procedure is outlined by Rothman, Levina and Zhu (2010) and is further described by Abegaz and Wit (2013). The procedure here is based on the work done in the R package SparseTSCGM (Abegaz and Wit, 2014).
}
\usage{
graphicalVAR(data, nLambda = 50, verbose = TRUE, gamma = 0.5, scale
                 = TRUE, lambda_beta, lambda_kappa, maxit.in = 100,
                 maxit.out = 100, deleteMissings = TRUE,
                 penalize.diagonal = TRUE, lambda_min_kappa = 0.05,
                 lambda_min_beta = lambda_min_kappa, mimic =
                 c("current", "0.1.2", "0.1.4", "0.1.5", "0.2"), vars,
                 beepvar, dayvar, idvar, lags = 1, centerWithin = TRUE,
                 likelihood = c("unpenalized", "penalized"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
A matrix or data frame containing repeated measures (rows) on a set of variables (columns). Must not contain missing data.
}
  \item{nLambda}{
The number of both lambda parameters to test. Defaults to 50, which results in 2500 models to evaluate.
}
  \item{verbose}{
Logical, should a progress bar be printed to the console?
}
  \item{gamma}{
The EBIC hyper-parameter. Set to 0 to use regular BIC.
}
\item{scale}{
Logical, should responses be standardized before estimation?
}
  \item{lambda_beta}{
An optional vector of lambda_beta values to test. Set \code{lambda_beta = 0} argument and \code{lambda_kappa = 0} for unregularized estimation.
}
  \item{lambda_kappa}{
An optional vector of lambda_kappa values to test. Set \code{lambda_beta = 0} argument and \code{lambda_kappa = 0} for unregularized estimation.
}
  \item{maxit.in}{
Maximum number of iterations in the inner loop (computing beta)
}
  \item{maxit.out}{
  Maximum number of iterations in the outer loop
  }
  \item{deleteMissings}{
    Logical, should missing responses be deleted?
  }
  \item{penalize.diagonal}{
  Logical, should the diagonal of beta be penalized (i.e., penalize auto-regressions)?
  }
  \item{lambda_min_kappa}{
      Multiplier of maximal tuning parameter for kappa
  }
    \item{lambda_min_beta}{
      Multiplier of maximal tuning parameter for beta
  }
  \item{mimic}{
      Allows one to mimic earlier versions of graphicalVAR
  }
\item{vars}{
Vectors of variables to include in the analysis
}
\item{beepvar}{
String indicating assessment beep per day (if missing, is added). Adding this argument will cause non-consecutive beeps to be treated as missing!
}
  \item{dayvar}{
String indicating assessment day. Adding this argument makes sure that the first measurement of a day is not regressed on the last measurement of the previous day. IMPORTANT: only add this if the data has multiple observations per day.
}
  \item{idvar}{
String indicating the subject ID
}
\item{lags}{
Vector of lags to include
}
\item{centerWithin}{
Logical, should subject data be within-person centered before estimating fixed effects?
}
\item{likelihood}{
Should likelihood be computed based on penalized contemporaneous matrix or unpenalized contemporaneous matrix. Set to \code{"penalized"} to mimic version 2.5 and later of \code{sparseTSCGM}.
}
}
\details{
Let y_t denote the vector of centered responses of a subject on a set of items on time point t. The graphical VAR model, using only one lag, is defined as follows:

y_t = Beta y_{y-1} + epsilon_t

In which epsilon_t is a vector of error and is independent between time points but not within time points. Within time points, the error is normally distributed with mean vector 0 and precision matrix (inverse covariance matrix) Kappa. The Beta matrix encodes the between time point interactions and the Kappa matrix encodes the within time point interactions. We aim to find a sparse solution for both Beta and Kappa, and do so by applying the LASSO algorithm as detailed by Rothman, Levina and Zhu (2010). The LASSO algorithm uses two tuning parameters, lambda_beta controlling the sparsity in Beta and lambda_kappa controlling the sparsity in Kappa. We estimate the model under a (by default) 50 by 50 grid of tuning parameters and choose the tuning parameters that optimize the extended Bayesian Information Criterion (EBIC; Chen and Chen,2008).

After estimation, the Beta and Kappa matrices can be standardized as described by Wild et al. (2010). The Kappa matrix can be standardized to partial contemporaneous correlations (PCC) as follows:

PCC(y_{i,t}, y_{j,t}) = - kappa_{ij} / ( sqrt{kappa_{ii} kappa_{jj}})

Similarly, the beta matrix can be standardized to partial directed correlations (PDC):

PDC(y_{i,t-1}, y_{j,t}) = beta_{ji} / sqrt{sigma_{jj} kappa_{ii} + beta_{ji}^2}

In which sigma is the inverse of kappa. Note that this process transposes the beta matrix. This is done because in representing a directed network it is typical to let rows indicate the node of origin and columns the node of destination.

Set \code{lambda_beta = 0} argument and \code{lambda_kappa = 0} for unregularized estimation.
}
\value{
A \code{graphicalVAR} object, which is a list containing:

\item{PCC}{The partial contemporaneous correlation network}
\item{PDC}{The  partial directed correlation network}
\item{beta}{The estimated beta matrix}
\item{kappa}{The estimated kappa matrix}
\item{EBIC}{The optimal EBIC}
\item{path}{Results of all tested tuning parameters}
\item{labels}{A vector containing the node labels}
}
\references{
Chen, J., & Chen, Z. (2008). Extended Bayesian information criteria for model selection with large model spaces. Biometrika, 95(3), 759-771.

Fentaw Abegaz and Ernst Wit (2013). Sparse time series chain graphical
models for reconstructing genetic networks. Biostatistics. 14, 3:
586-599.

Fentaw Abegaz and Ernst Wit (2014). SparseTSCGM: Sparse time series chain graphical
models. R package version 2.1.1. http://CRAN.R-project.org/package=SparseTSCGM

Rothman, A.J., Levina, E., and Zhu, J. (2010). Sparse multivariate
regression with covariance estimation. Journal of Computational and
Graphical Statistics. 19: 947-962.

Wild, B., Eichler, M., Friederich, H. C., Hartmann, M., Zipfel, S., &
Herzog, W. (2010). A graphical vector autoregressive modelling
approach to the analysis of electronic diary data. BMC medical
research methodology, 10(1), 28.
}
\author{
Sacha Epskamp <mail@sachaepskamp.com>
}

\examples{
# Simulate model:
Mod <- randomGVARmodel(4,probKappaEdge = 0.8,probBetaEdge = 0.8)

# Simulate data:
Data <- graphicalVARsim(100,Mod$beta,Mod$kappa)

# Estimate model:
Res <- graphicalVAR(Data, gamma = 0, nLambda = 10)

\donttest{
# Plot results:
layout(t(1:2))
plot(Mod, "PCC", layout = "circle")
plot(Res, "PCC", layout = "circle")

plot(Mod, "PDC", layout = "circle")
plot(Res, "PDC", layout = "circle")
}
}
