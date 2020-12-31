\name{ypsummary}
\alias{ypsummary}
\alias{ypsummary.default}
\alias{summary.ypsummary}
\alias{print.ypsummary}
\title{The main function of the package provides five summary measures of the treatment effect for clinical trials.}
\usage{
\method{ypsummary}{default}(time, event, group, tau, alpha = 0.05, tie = TRUE, 
    bound = 50, repnum = 2000, ...)
}
\arguments{
\item{...}{for S4 method only.}

\item{time}{A numeric vector of observations pooled from the two groups}

\item{event}{A numeric vector of the right-censoring indicator (event = 1, censored = 0)}

\item{group}{A numeric vector of the group indicator (treatment = 1, control = 0)}

\item{tau}{the upper end of the range used in defining the summary measures. Must be user-specified.}

\item{alpha}{Significance level. The default value is 0.05.}

\item{tie}{The default is \code{TRUE}. Add very tiny values to the observations when sorting them to avoid ties.}

\item{bound}{A boundary (\code{-bound}, \code{bound}) for estimating the parameters in the short-term and long-term hazard ratio model (Yang and Prentice, 2005). These parameters are \eqn{\beta_1} and \eqn{\beta_2} in their notations. The default boundary is \eqn{(-50, 50)}.}

\item{repnum}{the number of replications for the resampling method in obtaining the limiting variance estimators of the measures. The default value is 2000.}
}
\description{
\code{ypsummary} provides estimates of several summary measures of the treatment effect for design and analysis of clinical trials with survival outcomes, introduced in Yang (2018). The function utilizes the short-term and long-term hazard ratio model proposed in Yang and Prentice (2005), which can acccomodate various nonproportional hazard scenarios. The asymptotic properties of the summary measures are also discussed in Yang and Prentice (2011), Yang (2013), and Yang (2018).
}
\details{
The function \code{ypsummary} provides five summary measures of the treatment effect (see, Yang 2018), which can be utilized for various nonproportional haards scenarios:
\itemize{
\item the average hazard ratio (AHR): \eqn{\int  hr(x) dx}{∫  hr(x)dx}
\item the weighted average hazard ratio (WAHR): \eqn{\int hr(x) dw(x)}{∫  hr(x)dw(x)} where \eqn{dw(x) = dF_c(x)/F_c(\tau)}
\item the restricted superiority probability ratio (RSPR): \eqn{\frac{\int Sc(x) dF_t(x)}{\int S_t(x) dF_c(x)}}{∫  S_c(x)dF_t(x) / ∫  S_t(x)dF_c(x)}
\item the restricted mean survival difference (RMSD): \eqn{\int S_t(x) dx - \int S_c(x) dx}{∫  S_t(x)dx - ∫  S_c(x)dx}
\item the ratio of restricted mean times lost (RRMTL): \eqn{\frac{\tau - \int S_t(x) dx}{\tau - \int S_c(x) dx}}{(\tau - ∫  S_t(x)dx) / (\tau - ∫  S_c(x)dx)}
}

where \eqn{hr(x)} is the hazard ratio of the treatment group over the control group, \eqn{F_t(x)} and \eqn{F_c(x)} are the distribution functions for the treatment group and control group, respectively, and \eqn{S_t(x) = 1 - F_t(x)} and \eqn{S_c(x) = 1 - F_c(x)}. Note that all integrals are taken from 0 to \eqn{\tau}. 
}

\value{
  \item{Estimate}{The point estimate for the corresponding summary meausure}
  \item{CI}{The confidence interval constructed by a re-sampling method. If the measure is a ratio, the z-value is the standardized log of the estimate. If the measure is a difference, the z-value is the standardized estimate.}
  \item{z-value}{Normally distributed value derived from the asymptotic results}
  \item{p-value}{the (two-sided) p-value using \code{z-value}}
}
\references{
Yang, S. (2018). Improving testing and description of treatment effect in clinical trials with survival outcomes. Statistics in medicine.

Yang S, and Ross L. Prentice (2005). Semiparametric analysis of short-term and long-term hazard ratios with two-sample survival data. Biometrika, 92.1:1-17.

Yang S, and Ross L. Prentice (2011). Estimation of the 2-sample hazard ratio function using a semiparametric model. Biostatistics, 12.2:354-368.

Yang S. (2013). Semiparametric inference on the absolute risk reduction and the restricted mean survival difference in clinical trials. Special issue on risk assessment. Lifetime Data analysis, 19:219-241.
}

\examples{
library(ClinicalTrialSummary)
data(ggas)
time <- ggas$time
event <- ggas$event
group <- ggas$group
result <- ypsummary(time, event, group, tau=8.2) # tau must be supplied.
result
summary(result)
}
