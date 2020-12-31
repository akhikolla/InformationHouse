\name{ClinicalTrialSummary-package}
\alias{ClinicalTrialSummary-package}
\docType{package}
\title{
	Summary Measures for Clinical Trials with Survival Outcomes
}
\description{
ClinicalTrialSummary provides estimates of several summary measures of the treatment effect for design and analysis of clinical trials with survival outcomes, introduced in Yang (2018). These estimates are obtained under the short-term and long-term hazard ratio model (Yang and Prentice, 2005) which allows a range of time-varying hazard ratio shapes including crossing hazards situations.

Let \eqn{hr(x) = \lambda_{t}(x)/\lambda_{c}(x)}{hr(x) = λ_t(x)/λ_c(x)} be the hazard ratio function, where \eqn{\lambda_{t}(x)}{λ_t(x)} and \eqn{\lambda_{c}(x)}{λ_c(x)} are the hazard functions for the treatment group and the control group, respectively. The main function of the package, \code{\link{ypsummary}}, provides following five summary measures:
\itemize{
\item the average hazard ratio (AHR): \eqn{\int  hr(x) dx}{∫  hr(x)dx}
\item the weighted average hazard ratio (WAHR): \eqn{\int hr(x) dw(x)}{∫  hr(x)dw(x)} where \eqn{dw(x) = dF_c(x)/F_c(\tau)}
\item the restricted superiority probability ratio (RSPR): \eqn{\frac{\int S_c(x) dF_t(x)}{\int St(x) dF_c(x)}}{∫  S_c(x)dF_t(x) / ∫  S_t(x)dF_c(x)}
\item the restricted mean survival difference (RMSD): \eqn{\int S_t(x) dx - \int S_c(x) dx}{∫  S_t(x)dx - ∫  S_c(x)dx}
\item the ratio of restricted mean times lost (RRMTL): \eqn{\frac{\tau - \int S_t(x) dx}{\tau - \int S_c(x) dx}}{(\tau - ∫  S_t(x)dx) / (\tau - ∫  S_c(x)dx)}
}

Note that  all integrals are taken over \eqn{(0, \tau)} for \eqn{\tau} less than or equal to the maximum follow-up duration of the trial. The asymptotic results for the average hazard ratio and the restricted mean survival difference were established in Yang and Prentice (2011) and Yang (2013), respectively. For other measures, the asymptotic results were established in Yang (2018).

The data used as an example in this package was from Gastrointestinal Tumor Study Group (1982) and the object returned by \code{\link{ypsummary}} can be formatted to a table using the function \code{\link{summary}}.
}
\details{
\tabular{ll}{
Package: \tab ClinicalTrialSummary\cr
Type: \tab Package\cr
Version: \tab 1.0.0\cr
Date: \tab 2018-12-19\cr
License: \tab GPL (>= 3)\cr
}
}
\author{
Daewoo Pak and Song Yang
}
\references{
Yang, S. (2018). Improving testing and description of treatment effect in clinical trials with survival outcomes. Statistics in medicine.

Yang S, and Ross L. Prentice (2005). Semiparametric analysis of short-term and long-term hazard ratios with two-sample survival data. Biometrika, 92.1:1-17.

Yang S, and Ross L. Prentice (2011). Estimation of the 2-sample hazard ratio function using a semiparametric model. Biostatistics, 12.2:354-368.

Yang S. (2013). Semiparametric inference on the absolute risk reduction and the restricted mean survival difference in clinical trials. Special issue on risk assessment. Lifetime Data analysis, 19:219-241.

Gastrointestinal Tumor Study Group (1982). A comparison of combination chemotherapy and combined modality therapy for locally advanced gastric carcinoma. Cancer.
}

\keyword{Survival analysis, Clinical trials, Hazard ratio, Yang and Prentice model}
\seealso{
\code{\link{ypsummary}}
}
\examples{
	library(ClinicalTrialSummary)
	data(ggas)
	result <- ypsummary(time=ggas$time, event=ggas$event, group=ggas$group, tau = 8.2)
	print(result)
	summary(result)
}
