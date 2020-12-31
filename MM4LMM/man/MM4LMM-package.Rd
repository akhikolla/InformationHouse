\name{MM4LMM-package}
\alias{MM4LMM-package}
\alias{MM4LMM}
\docType{package}
\title{
  Min-Max algorithms for Variance Component Mixed Model Inference.
}
\description{
This package provides a function to perform either ML or ReML estimation in a Variance Component Mixed Model. The optimization of the (possibly Restricted) log-likelihood is perfomed using the Min-Max algorithm described in Hunter et al. (2004). Depending on the number of variance components, different computational tricks are used to speed up inference. Additionally, the \code{\link{AnovaTest}} function provides  Type I, Type III and Type III Kenward Roger approximation test series for fixed effects. The nullity of a specific linear combination can also be tested.  
}
\details{
  Variance Component Mixed Models are mixed models of the form
  \deqn{Y = X \beta +  \sum_{k=1}^K Z_k u_k}
  where Y is the response vector, X and \eqn{\beta} are respectively the incidence matrix and the coefficient vector associated with the fixed effects, \eqn{u_k} is the kth vector of random effects and corresponds to its associated incidence matrix. All random effect are assumed to be Gaussian with mean 0 and covariance \eqn{\sigma_k^2 R_k }, where \eqn{R_k} is a known correlation matrix and  \eqn{\sigma_k^2} is an unknown variance parameter. All random effects are assumed to be independent. In many applications the last random component corresponds to the error and therefore both \eqn{Z_k} and \eqn{R_k} correspond to the identity matrix. 
  
  In such models the inference of both the unknown variance components \eqn{\sigma_k^2},...,\eqn{\sigma_K^2} and the fixed effect \eqn{\beta} can be achieved through Maximum Likelihood (ML) or Restricted Maximum Likelihood (ReML) estimation. Neither ML nor ReML yield close form expressions of the estimates, consequently the maximization of the (restricted) log-likelihood has to be performed numerically. This package provides the user with Min-Max algorithms for the optimization.  Efficient tricks such as profiling, MME trick and MM acceleration are implemented for computational efficiency (see Johnson et al. (1995), Varadhan et al. (2008) for details). The main function \code{\link{MMEst}} handles parallel inference of multiple models sharing the same set of correlation matrices - but possibly different fixed effects, an usual situation in GWAS analysis for instance.
}
\author{
Fabien Laporte and Tristan Mary-Huard

Maintainer: Fabien Laporte <fabien.laporte@pasteur.fr>
}
\references{
Johnson, D. L., & Thompson, R. (1995). Restricted maximum likelihood estimation of variance components for univariate animal models using sparse matrix techniques and average information. Journal of dairy science, 78(2), 449-456.

Hunter, D. R., & Lange, K. (2004). A tutorial on MM algorithms. The American Statistician, 58(1), 30-37.

Kenward, M. G., & Roger, J. H. (1997). Small sample inference for fixed effects from restricted maximum likelihood. Biometrics, 983-997.

Varadhan, R., & Roland, C. (2008). Simple and globally convergent methods for accelerating the convergence of any EM algorithm. Scandinavian Journal of Statistics, 35(2), 335-353.

Zhou, H., Hu, L., Zhou, J., & Lange, K. (2015). MM algorithms for variance components models. arXiv preprint arXiv:1509.07426.
}
\keyword{ package }
%%\seealso{
%%  Optional links to other man pages
%%}
%%\examples{
%%  \dontrun{
%%     ## Optional simple examples of the most important functions
%%     ## These can be in \dontrun{} and \donttest{} blocks.   
%%  }
%%}
