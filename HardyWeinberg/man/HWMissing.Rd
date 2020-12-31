\name{HWMissing}
\alias{HWMissing}
\title{
  Test a bi-allelic marker for Hardy-Weinberg equilibrium in the
  presence of missing genotype information.
}
\description{
  Function \code{HWMissing} imputes missing genotype data with a
  multinomial logit model that uses information from allele intensities
  and/or neighbouring markers. Multiple imputation algorithms
  implemented in the Mice package are used to obtain imputed data sets.
  Inference for HWE is carried out by estimating the inbreeding
  coefficient or exact p-values for each imputed data set, and
  by combining all estimates
  using Rubin's pooling rules.
}
\usage{
HWMissing(X, imputecolumn = 1, m = 50, coding = c(0,1,2), verbose = FALSE, alpha = 0.05,
    varest = "oneovern", statistic = "chisquare",  alternative =
"two.sided", ...)
}
\arguments{
  \item{X}{An input data frame. By default, the first column should
    contain the SNP with missing values.}
  \item{imputecolumn}{Indicates which column of the supplied data frame
    is to be imputed (by default, the first colum, \code{imputecolumn=1}}
  \item{m}{The number of imputations (50 by default)}
  \item{coding}{Indicates how the genotype data is coded (e.g. 0 for AA,
1 for AB, and 2 for BB).}
  \item{verbose}{ \code{verbose = TRUE} prints results, \code{verbose = FALSE} is silent.}
   \item{alpha}{significance level (0.05 by default) used when computing
    confidence intervals}
  \item{varest}{Estimator for the variance of the inbreeding
    coefficient. \code{varest="oneovern"} is the default and sets the
    variance under the null (1/n). \code{varest="bailey"} uses an
    approximation (see details).}
  \item{statistic}{If \code{statistic = "chisquare"} then inbreeding
    coefficients (equivalent to chisquare statistics) will be computed
    for each imputed data set and then combined. If \code{statistic =
      "exact"} then one-sided exact tests will be computed for each
    imputed data set and the resulting p-values will be combined.}
  \item{alternative}{\code{two.sided} (default) will perform a two-sided
    test where both an excess and a dearth of heterozygotes count as
    evidence against HWE. \code{less} is a one-sided test where only
    dearth of heterozygotes counts a evidence against HWE,
    \code{greater} is a one-sided test where only excess of
    heterozygotes counts as evidence against HWE.}
  \item{\ldots}{additional options for function \code{mice} of the Mice package}
}
\details{
   The function \code{HWMissing} tests one genetic marker (e.g. a SNP)
   with missings for HWE. By default, this marker is supposed to be the
   first column of dataframe \code{X}. The other columns of \code{X}
   contain covariates to be used in the imputation model. Covariates
   will typically be other, correlated markers or allele intensities of
   the SNP to be imputed. Covariate markers should be coded as factor
   variables whereas allele intensities should be numerical
   variables. By default, a polytomous regression model will be used to
   impute the missings. If the covariates also contain missings, an
   imputation method for each column of \code{X} can be specified by
   using the \code{method} of mice (see example below).

   If there are no covariates, missings can be imputed under the MCAR
   assumption. In that case, missings are imputed by taking a random
   sample from the observed data. This is what \code{HWMissing} will do
   if no covariates are supplied, \code{X} being a single factor
   variable.

   Several estimators for the variance of the inbreeding coefficient
   have been described in the literature. The asymptotic variance of the
   inbreeding coefficient under the null hypothesis is 1/n, and is used
   if \code{varest = "oneovern"} is used. This is the recommended
   option. Alternatively, the approximation described in Weir (p. 66) can be used
   with \code{varest = "bailey"}.
   }
\value{
  \item{Res}{A vector with the inbreeding coefficient, a confidence interval
  for the inbreeding coefficient, a p-value for a HWE test and missing
  data statistics.}
  \item{Xmat}{A matrix with the genotypic composition of each of the
   \code{m} imputed data sets.}
}

\references{
Little, R. J. A. and Rubin, D. B. (2002) Statistical analysis with
missing data. Second edition, New York, John Wiley & sons.

Graffelman, J., S\'anchez, M., Cook, S. and Moreno, V. (2013) Statistical
inference for Hardy-Weinberg proportions in the presence of missing
genotype information. PLoS ONE 8(12):
e83316. \url{http://dx.doi.org/10.1371/journal.pone.0083316}

Graffelman, J. (2015) Exploring Diallelic Genetic Markers: The HardyWeinberg Package.
\emph{Journal of Statistical Software} 64(3): 1-23. \url{https://www.jstatsoft.org/v64/i03/}.
}
\author{ Jan Graffelman \email{jan.graffelman@upc.edu} }
\seealso{ \code{\link{HWChisq}} }
\examples{
data(Markers)
\dontrun{
set.seed(123)
Results <- HWMissing(Markers[,1],m=50,verbose=TRUE)$Res # no covariates, imputation assuming MCAR.
set.seed(123)
Results <- HWMissing(Markers[,1:3],m=50,verbose=TRUE)$Res # impute with two allele intensities.
set.seed(123)
Results <- HWMissing(Markers[,c(1,4,5)],m=50,verbose=TRUE)$Res # impute with two covariate SNPs
}
}
\keyword{misc}

