%% File Name: mnlfa.Rd
%% File Version: 0.336

\name{mnlfa}
\alias{mnlfa}
\alias{summary.mnlfa}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{
Moderated Nonlinear Factor Analysis
}

\description{
General function for conducting moderated nonlinear factor analysis
(Curran et al., 2014). Item slopes and item intercepts can be modeled as functions
of person covariates.

Parameter regularization is allowed. For categorical covariates,
group lasso can be used for regularization.
}


\usage{
mnlfa(dat, items, item_type="2PL", formula_int=~1, formula_slo=~1, formula_mean=~0,
   formula_sd=~0, theta=NULL, parm_list_init=NULL, parm_trait_init=NULL, prior_init=NULL,
   regular_lam=c(0, 0), regular_type=c("none", "none"), maxit=1000, msteps=4, conv=1e-05,
   conv_mstep=1e-04, h=1e-04, parms_regular_types=NULL, parms_regular_lam=NULL,
   parms_iterations=NULL, center_parms=NULL, center_max_iter=6, L_max=.07,
   verbose=TRUE)

\method{summary}{mnlfa}(object, file=NULL, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
Data frame with item responses
}
  \item{items}{
Vector containing item names
}
  \item{item_type}{
String or vector of item types. Currently, only item types \code{"1PL"} or
\code{"2PL"} can be chosen.
}
  \item{formula_int}{
String or list with formula for item intercepts
}
  \item{formula_slo}{
String or list with formula for item slopes
}
  \item{formula_mean}{
Formula for mean of the trait distribution
}
  \item{formula_sd}{
Formula for standard deviation of the trait distribution
}
  \item{theta}{
Grid of \eqn{\theta} values used for approximation of normally
distributed trait
}
  \item{parm_list_init}{
Optional list of initial item parameters
}
\item{parm_trait_init}{Optional list of initial parameters for trait distribution}
\item{prior_init}{Optional matrix of prior distribution for persons}
  \item{regular_lam}{
Vector of length two containing two general regularization
parameters for item intercepts and item slopes
}
  \item{regular_type}{
Type of regularization method. Can be \code{"none"},
\code{"lasso"}, \code{"scad"} or \code{"mcp"}.
}
  \item{maxit}{
Maximum number of iterations
}
  \item{msteps}{
Maximum number of M-steps
}
  \item{conv}{
Convergence criterion with respect to parameters
}
  \item{conv_mstep}{
Convergence criterion in M-step
}
  \item{h}{
Numerical differentiation parameter
}
  \item{parms_regular_types}{
Optional list containing parameter specific regularization types
}
  \item{parms_regular_lam}{
Optional list containing parameter specific regularization parameters
}
  \item{parms_iterations}{
Optional list containing sequence of parameter indices used for updating
}
\item{center_parms}{Optional list indicating which parameters
should be centered during initial iterations.}
\item{center_max_iter}{Maximum number of iterations in which parameters
  should be centered.}
\item{L_max}{Majorization parameter used in regularization}
\item{verbose}{
Logical indicating whether output should be printed
}
\item{object}{Object of class \code{mnlfa}}
\item{file}{Optional file name}
\item{\dots}{Further arguments to be passed}
}

\details{
The moderated factor analysis model for dichotomous responses defined as
\deqn{P(X_{pi}=1 | \theta_p )=invlogit( a_{pi} \theta_p - b_{pi} ) }
The trait distribution \eqn{\theta_p \sim N( \mu_p, \sigma_p^2)}
allows a latent regression of person covariates on the mean
with \eqn{\mu_p=\bold{X}_p \bold{\gamma}} (to be specified in \code{formula_mean})
and the logarithm of the standard deviation \eqn{\log \sigma_p=\bold{Z}_p \bold{\delta} }
(to be specified in \code{formula_sd}).
Item intercepts and item slopes can be moderated by person covariates, i.e.
\eqn{a_{pi}=\bold{W}_{pi} \bold{\alpha}_i } and
\eqn{b_{pi}=\bold{V}_{pi} \bold{\beta}_i }. Regularization on (some of) the
\eqn{\bold{\alpha}_i} or \eqn{\bold{\beta}_i} parameters is allowed.

The model is estimated using an EM algorithm with the coordinate descent
method during the M-step (Sun et al., 2016).
}


\value{
List with model results including

\item{item}{Summary table for item parameters}
\item{trait}{Summary table for trait parameters}
}


\references{
Curran, P. J., McGinley, J. S., Bauer, D. J., Hussong, A. M., Burns, A.,
Chassin, L., Sher, K., & Zucker, R. (2014).
A moderated nonlinear factor model for the development of commensurate measures
in integrative data analysis. \emph{Multivariate Behavioral Research, 49}(3),
214-231.
http://dx.doi.org/10.1080/00273171.2014.889594

Sun, J., Chen, Y., Liu, J., Ying, Z., & Xin, T. (2016). Latent variable selection for 
multidimensional item response theory models via L1 regularization. 
\emph{Psychometrika, 81}(4), 921-939.
 https://doi.org/10.1007/s11336-016-9529-6
}


%\author{
%%  ~~who you are~~
%}

%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See also the \pkg{aMNLFA} package for automatized moderated nonlinear factor analysis
which provides convenient wrapper functions for automized analysis
in the \emph{Mplus} software.

See the \pkg{GPCMlasso} package for the regularized generalized partial credit model.
}

\examples{
#############################################################################
# EXAMPLE 1: Dichotomous data, 1PL model
#############################################################################

data(data.mnlfa01, package="mnlfa")

dat <- data.mnlfa01
# extract items from dataset
items <- grep("I[0-9]", colnames(dat), value=TRUE)
I <- length(items)

# maximum number of iterations (use only few iterations for the only purpose of
# providing CRAN checks)
maxit <- 10

#***** Model 1: 1PL model without moderating parameters and without covariates for traits

# no covariates for trait
formula_mean <- ~0
formula_sd <- ~1
# no item covariates
formula_int <- ~1
formula_slo <- ~1

mod1 <- mnlfa::mnlfa( dat=dat, items, item_type="1PL", formula_int=formula_int,
             formula_slo=formula_slo, formula_mean=formula_mean, formula_sd=formula_sd,
             maxit=maxit )
summary(mod1)

\donttest{
#***** Model 2: 1PL model without moderating parameters and with covariates for traits

# covariates for trait
formula_mean <- ~female + age
formula_sd <- ~1

mod2 <- mnlfa::mnlfa( dat=dat, items, item_type="1PL", formula_int=formula_int,
             formula_slo=formula_slo, formula_mean=formula_mean, formula_sd=formula_sd)
summary(mod2)

#***** Model 3: 1PL model with moderating parameters and with covariates for traits
#***   Regularization method 'mcp'

# covariates for trait
formula_mean <- ~female + age
formula_sd <- ~1
# moderation effects for items
formula_int <- ~1+female+age
formula_slo <- ~1

# center parameters for female and age in initial iterations for improving convergence
center_parms <- list( rep(2,I), rep(3,I) )

# regularization parameters for item intercept and item slope, respectively
regular_lam <- c(.06, .25)
regular_type <- c("mcp","none")

mod3 <- mnlfa::mnlfa( dat=dat, items, item_type="1PL", formula_int=formula_int,
            formula_slo=formula_slo, formula_mean=formula_mean, formula_sd=formula_sd,
            center_parms=center_parms, regular_lam=regular_lam, regular_type=regular_type )
summary(mod3)


#***** Model 4: 1PL model with selected moderated item parameters

#* trait distribution
formula_mean <- ~0+female+age
formula_sd <- ~1

#* formulas for item intercepts
formula_int <- ~1
formula_int <- mnlfa::mnlfa_expand_to_list(x=formula_int, names_list=items)
mod_items <- c(4,5,6,7)
for (ii in mod_items){
    formula_int[[ii]] <- ~1+female+age
}
formula_slo <- ~1

mod4 <- mnlfa::mnlfa( dat=dat, items, item_type="1PL", formula_int=formula_int,
              formula_slo=formula_slo, formula_mean=formula_mean, formula_sd=formula_sd)
mod4$item
mod4$trait
summary(mod4)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
