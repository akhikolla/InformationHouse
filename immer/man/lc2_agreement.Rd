%% File Name: lc2_agreement.Rd
%% File Version: 0.21

\name{lc2_agreement}
\alias{lc2_agreement}
\alias{summary.lc2_agreement}
\alias{logLik.lc2_agreement}
\alias{anova.lc2_agreement}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A Latent Class Model for Agreement of Two Raters
}

\description{
Estimates a latent class model for agreement of two raters
(Schuster & Smith, 2006). See Details for the description of
the model.
}

\usage{
lc2_agreement(y, w=rep(1, nrow(y)), type="homo", method="BFGS", ...)

\method{summary}{lc2_agreement}(object, digits=3,...)

\method{logLik}{lc2_agreement}(object, ...)

\method{anova}{lc2_agreement}(object, ...)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
A data frame containing the values of two raters in columns
}
  \item{w}{
Optional vector of weights
}
\item{type}{
Type of model specification. Can be \code{"unif"}, \code{"equal"},
  \code{"homo"} or \code{"hete"}. See Details.
}
  \item{method}{
Optimization method used in \code{\link[stats:optim]{stats::optim}}
}
  \item{\dots}{
Further arguments passed to \code{\link[stats:optim]{stats::optim}}
}
\item{object}{Object of class \code{l2_agreement}}
\item{digits}{Number of digits for rounding}
}

\details{
The latent class model for two raters decomposes a portion of ratings
which conform to true agreement and another portion of ratings which
conform to a random rating of a category. Let \eqn{X_r} denote the rating of
rater \eqn{r}, then for \eqn{ i \neq j}, it is assumed that
\deqn{P(X_1=i, X_2=j)=\phi_{1i}  \phi_{2j} ( 1 - \gamma )}
For \eqn{i=j} it is assumed that
\deqn{P(X_1=i, X_2=i)=\tau_i \gamma + \phi_{1i} \phi_{2i} ( 1 - \gamma )}
where \eqn{\gamma} denotes the proportion of true ratings.

All \eqn{\tau_i} and \eqn{\phi_{ri}} parameters are estimated
using \code{type="hete"}. If the \eqn{\phi} parameters are assumed
as invariant across the two raters (i.e. \eqn{\phi_{1i}=\phi_{2i}=\phi_{i}}),
then \code{type="homo"} must be specified. The constraint \eqn{\tau_i=\phi_i}
is imposed by \code{type="equal"}. All \eqn{\phi_i} parameters
are set equal to each other using \code{type="unif"}.
}

\value{
\item{model_output}{Output of the fitted model}
\item{saturated_output}{Output of the saturated model}
\item{LRT_output}{Output of the likelihood ratio test of model fit}
\item{partable}{Parameter table}
\item{parmsummary}{Parameter summary}
\item{agree_true}{True agreement index shich is the \eqn{\gamma} parameter}
\item{agree_chance}{Agreement by chance}
\item{rel_agree}{Conditional reliability of agreement}
\item{optim_output}{Output of \code{optim} from the fitted model}
\item{nobs}{Number of observations}
\item{type}{Model type}
\item{ic}{Information criteria}
\item{loglike}{Log-likelihood}
\item{npars}{Number of parameters}
\item{y}{Used dataset}
\item{w}{Used weights}
}

\references{
Schuster, C., & Smith, D. A. (2006). Estimating with a latent
class model the reliability of nominal judgments upon which
two raters agree. \emph{Educational and Psychological Measurement, 66}(5), 739-747.
}


%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}

\examples{
#############################################################################
# EXAMPLE 1: Dataset in Schuster and Smith (2006)
#############################################################################

data(data.immer08)
dat <- data.immer08

# select ratings and frequency weights
y <- dat[,1:2]
w <- dat[,3]

#*** Model 1: Uniform distribution phi parameters
mod1 <- immer::lc2_agreement( y=y, w=w, type="unif")
summary(mod1)

#*** Model 2: Equal phi and tau parameters
mod2 <- immer::lc2_agreement( y=y, w=w, type="equal")
summary(mod2)

\dontrun{
#*** Model 3: Homogeneous rater model
mod3 <- immer::lc2_agreement( y=y, w=w, type="homo")
summary(mod3)

#*** Model 4: Heterogeneous rater model
mod4 <- immer::lc2_agreement( y=y, w=w, type="hete")
summary(mod4)

#--- some model comparisons
anova(mod3,mod4)
IRT.compareModels(mod1,mod2,mod3,mod4)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
