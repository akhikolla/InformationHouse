%% File Name: immer_cml.Rd
%% File Version: 0.42

\name{immer_cml}
\alias{immer_cml}
\alias{summary.immer_cml}
\alias{logLik.immer_cml}
\alias{anova.immer_cml}
\alias{coef.immer_cml}
\alias{vcov.immer_cml}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Conditional Maximum Likelihood Estimation for the Linear Logistic
Partial Credit Model
}

\description{
Conditional maximum likelihood estimation for the linear logistic
partial credit model (Molenaar, 1995; Andersen, 1995; Fischer, 1995).
The \code{immer_cml} function allows for known
integer discrimination parameters like in the one-parameter logistic
model (Verhelst & Glas, 1995).
}

\usage{
immer_cml(dat, weights=NULL, W=NULL, b_const=NULL, par_init=NULL,
    a=NULL, irtmodel=NULL, normalization="first", nullcats="zeroprob",
    diff=FALSE, use_rcpp=FALSE, ...)

\method{summary}{immer_cml}(object, digits=3, file=NULL, ...)

\method{logLik}{immer_cml}(object,...)

\method{anova}{immer_cml}(object,...)

\method{coef}{immer_cml}(object,...)

\method{vcov}{immer_cml}(object,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
Data frame with item responses
}
  \item{weights}{
Optional vector of sample weights
}
  \item{W}{
Design matrix \eqn{\bold{W}} for linear logistic partial credit model. Every
row corresponds to a parameter for item \eqn{i} in category \eqn{h}
}
  \item{b_const}{
Optional vector of parameter constants \eqn{b_{0ih}} which can be
used for parameter fixings.
}
  \item{par_init}{
Optional vector of initial parameter estimates
}
  \item{a}{
Optional vector of integer item discriminations
}
  \item{irtmodel}{
Type of item response model. \code{irtmodel="PCM"} and \code{irtmodel="PCM2"}
follow the conventions of the \pkg{TAM} package.
}
  \item{normalization}{
The type of normalization in partial credit models. Can be \code{"first"}
for the first item or \code{"sum"} for a sum constraint.
}
  \item{nullcats}{
A string indicating whether categories with zero frequencies should
have a probability of zero (by fixing the constant parameter to a large
value of \code{99}).
}
  \item{diff}{
Logical indicating whether the difference algorithm should be used. See
\code{\link[psychotools:elementary_symmetric_functions]{psychotools::elementary_symmetric_functions}}
for details.
}
\item{use_rcpp}{Logical indicating whether \pkg{Rcpp} package should be
used for computation.}
\item{\dots}{
Further arguments to be passed to \code{\link[stats:optim]{stats::optim}}.
}
\item{object}{Object of class \code{immer_cml}}
\item{digits}{Number of digits after decimal to be rounded.}
\item{file}{Name of a file in which the output should be sunk}
}

\details{
The partial credit model can be written as
\deqn{P(X_{pi}=h ) \propto \exp( a_i h \theta_p - b_{ih}) }
where the item-category parameters \eqn{b_{ih}} are linearly
decomposed according to
\deqn{ b_{ih}=\sum_{v} w_{ihv} \beta_v + b_{0ih}}
with unknown basis parameters \eqn{\beta_v} and fixed values \eqn{w_{ihv}}
of the design matrix \eqn{\bold{W}} (specified in \code{W})
and constants \eqn{b_{0ih}} (specified in \code{b_const}).
}

\value{
List with following entries:

\item{item}{Data frame with item-category parameters}
\item{b}{Item-category parameters \eqn{b_{ih}}}
\item{coefficients}{Estimated basis parameters \eqn{\beta_{v}}}
\item{vcov}{Covariance matrix of basis parameters \eqn{\beta_{v}}}
\item{par_summary}{Summary for basis parameters}
\item{loglike}{Value of conditional log-likelihood}
\item{deviance}{Deviance}
\item{result_optim}{Result from optimization in \code{\link[stats:optim]{stats::optim}}}
\item{W}{Used design matrix \eqn{\bold{W}}}
\item{b_const}{Used constant vector \eqn{b_{0ih}}}
\item{par_init}{Used initial parameters}
\item{suffstat}{Sufficient statistics}
\item{score_freq}{Score frequencies}
\item{dat}{Used dataset}
\item{used_persons}{Used persons}
\item{NP}{Number of missing data patterns}
\item{N}{Number of persons}
\item{I}{Number of items}
\item{maxK}{Maximum number of categories per item}
\item{K}{Maximum score of all items}
\item{npars}{Number of estimated parameters}
\item{pars_info}{Information of definition of item-category parameters \eqn{b_{ih}}}
\item{parm_index}{Parameter indices}
\item{item_index}{Item indices}
\item{score}{Raw score for each person}
}



\references{
Andersen, E. B. (1995). Polytomous Rasch models and their estimation.
In G. H. Fischer & I. W. Molenaar (Eds.). \emph{Rasch Models} (pp. 39--52).
New York: Springer.

Fischer, G. H. (1995). The linear logistic test model.
In G. H. Fischer & I. W. Molenaar (Eds.). \emph{Rasch Models} (pp. 131--156).
New York: Springer.

Molenaar, I. W. (1995). Estimation of item parameters.
In G. H. Fischer & I. W. Molenaar (Eds.). \emph{Rasch Models} (pp. 39--52).
New York: Springer.

Verhelst, N. D. &, Glas, C. A. W. (1995). The one-parameter logistic model.
In G. H. Fischer & I. W. Molenaar (Eds.). \emph{Rasch Models} (pp. 215--238).
New York: Springer.
}




%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
For CML estimation see also the \pkg{eRm} and \pkg{psychotools} packages and the
functions \code{\link[eRm:RM]{eRm::RM}} and
\code{\link[psychotools:raschmodel]{psychotools::raschmodel}} for the Rasch model
and \code{\link[eRm:PCM]{eRm::PCM}} and
\code{\link[psychotools:pcmodel]{psychotools::pcmodel}} for the partial
credit model.

See \code{\link[eRm:LLTM]{eRm::LLTM}} for the linear logistic test model
and \code{\link[eRm:LPCM]{eRm::LPCM}} for the linear logistic partial
credit model in the \pkg{eRm} package for CML implementations.

The \code{immer_cml} function makes use of
\code{\link[psychotools:elementary_symmetric_functions]{psychotools::elementary_symmetric_functions}}.

For CML estimation with sample weights see also the
\pkg{RM.weights} package.
}


\examples{
library(sirt)
library(psychotools)
library(TAM)
library(CDM)
library(eRm)

#############################################################################
# EXAMPLE 1: Dichotomous data data.read
#############################################################################

data(data.read, package="sirt")
dat <- data.read
I <- ncol(dat)

#----------------------------------------------------------------
#--- Model 1: Rasch model, setting first item difficulty to zero
mod1a <- immer::immer_cml( dat=dat)
summary(mod1a)
logLik(mod1a) # extract log likelihood
coef(mod1a)   # extract coefficients

\dontrun{
# estimate model in psychotools package
mod1b <- psychotools::raschmodel(dat)
summary(mod1b)
logLik(mod1b)

# estimate model in eRm package
mod1c <- eRm::RM(dat, sum0=FALSE)
summary(mod1c)
mod1c$etapar

# compare estimates of three packages
cbind( coef(mod1a), coef(mod1b), mod1c$etapar )

#----------------------------------------------------------------
#-- Model 2: Rasch model sum normalization
mod2a <- immer::immer_cml( dat=dat, normalization="sum")
summary(mod2a)

# compare estimation in TAM
mod2b <- tam.mml( dat, constraint="items"  )
summary(mod2b)
mod2b$A[,2,]

#----------------------------------------------------------------
#--- Model 3: some fixed item parameters
# fix item difficulties of items 1,4,8
# define fixed parameters in constant parameter vector
b_const <- rep(0,I)
fix_items <- c(1,4,8)
b_const[ fix_items ] <- c( -2.1, .195, -.95 )
# design matrix
W <- matrix( 0, nrow=12, ncol=9)
W[ cbind( setdiff( 1:12, fix_items ), 1:9 ) ] <- 1
colnames(W) <- colnames(dat)[ - fix_items ]
# estimate model
mod3 <- immer::immer_cml( dat=dat, W=W, b_const=b_const)
summary(mod3)

#----------------------------------------------------------------
#--- Model 4: One parameter logistic model
# estimate non-integer item discriminations with 2PL model
I <- ncol(dat)
mod4a <- sirt::rasch.mml2( dat, est.a=1:I )
summary(mod4a)
a <- mod4a$item$a     # extract (non-integer) item discriminations
# estimate integer item discriminations ranging from 1 to 3
a_integer <- immer::immer_opcat( a, hmean=2, min=1, max=3 )
# estimate one-parameter model with fixed integer item discriminations
mod4 <- immer::immer_cml( dat=dat, a=a_integer )
summary(mod4)

#----------------------------------------------------------------
#--- Model 5: Linear logistic test model

# define design matrix
W <- matrix( 0, nrow=12, ncol=5 )
colnames(W) <- c("B","C", paste0("Pos", 2:4))
rownames(W) <- colnames(dat)
W[ 5:8, "B" ] <- 1
W[ 9:12, "C" ] <- 1
W[ c(2,6,10), "Pos2" ] <- 1
W[ c(3,7,11), "Pos3" ] <- 1
W[ c(4,8,12), "Pos4" ] <- 1

# estimation with immer_cml
mod5a <- immer::immer_cml( dat, W=W )
summary(mod5a)

# estimation in eRm package
mod5b <- eRm::LLTM( dat, W=W )
summary(mod5b)

# compare models 1 and 5 by a likelihood ratio test
anova( mod1a, mod5a )

#############################################################################
# EXAMPLE 2: Polytomous data | data.Students
#############################################################################

data(data.Students,package="CDM")
dat <- data.Students
dat <- dat[, grep("act", colnames(dat) ) ]
dat <- dat[1:400,]  # select a subdataset
dat <- dat[ rowSums( 1 - is.na(dat) ) > 1, ]
    # remove persons with less than two valid responses

#----------------------------------------------------------------
#--- Model 1: Partial credit model with constraint on first parameter
mod1a <- immer::immer_cml( dat=dat )
summary(mod1a)
# compare pcmodel function from psychotools package
mod1b <- psychotools::pcmodel( dat )
summary(mod1b)
# estimation in eRm package
mod1c <- eRm::PCM( dat, sum0=FALSE )
  # -> subjects with only one valid response must be removed
summary(mod1c)

#----------------------------------------------------------------
#-- Model 2: Partial credit model with sum constraint on item difficulties
mod2a <- immer::immer_cml( dat=dat, irtmodel="PCM2", normalization="sum")
summary(mod2a)
# compare with estimation in TAM
mod2b <- TAM::tam.mml( dat, irtmodel="PCM2", constraint="items")
summary(mod2b)

#----------------------------------------------------------------
#-- Model 3: Partial credit model with fixed integer item discriminations
mod3 <- immer::immer_cml( dat=dat, normalization="first", a=c(2,2,1,3,1) )
summary(mod3)

#############################################################################
# EXAMPLE 3: Polytomous data | Extracting the structure of W matrix
#############################################################################

data(data.mixed1, package="sirt")
dat <- data.mixed1

# use non-exported function "lpcm_data_prep" to extract the meaning
# of the rows in W which are contained in value "pars_info"
res <- immer:::lpcm_data_prep( dat, weights=NULL, a=NULL )
pi2 <- res$pars_info

# create design matrix with some restrictions on item parameters
W <- matrix( 0, nrow=nrow(pi2), ncol=2 )
colnames(W) <- c( "P2", "P3" )
rownames(W) <- res$parnames

# joint item parameter for items I19 and I20 fixed at zero
# item parameter items I21 and I22
W[ 3:10, 1 ] <- pi2$cat[ 3:10 ]
# item parameters I23, I24 and I25
W[ 11:13, 2] <- 1

# estimate model with design matrix W
mod <- immer::immer_cml( dat, W=W)
summary(mod)

#############################################################################
# EXAMPLE 4: Partial credit model with raters
#############################################################################

data(data.immer07)
dat <- data.immer07

#*** reshape dataset for one variable
dfr1 <- immer::immer_reshape_wideformat( dat$I1, rater=dat$rater, pid=dat$pid )

#-- extract structure of design matrix
res <- immer:::lpcm_data_prep( dat=dfr1[,-1], weights=NULL, a=NULL)
pars_info <- res$pars_info

# specify design matrix for partial credit model and main rater effects
# -> set sum of all rater effects to zero
W <- matrix( 0, nrow=nrow(pars_info), ncol=3+2 )
rownames(W) <- rownames(pars_info)
colnames(W) <- c( "Cat1", "Cat2", "Cat3", "R1", "R2" )
# define item parameters
W[ cbind( pars_info$index, pars_info$cat ) ] <- 1
# define rater parameters
W[ paste(pars_info$item)=="R1", "R1" ] <- 1
W[ paste(pars_info$item)=="R2", "R2" ] <- 1
W[ paste(pars_info$item)=="R3", c("R1","R2") ] <- -1
# set parameter of first category to zero for identification constraints
W <- W[,-1]

# estimate model
mod <- immer::immer_cml( dfr1[,-1], W=W)
summary(mod)

#############################################################################
# EXAMPLE 5: Multi-faceted Rasch model | Estimation with a design matrix
#############################################################################

data(data.immer07)
dat <- data.immer07

#*** reshape dataset
dfr1 <- immer::immer_reshape_wideformat( dat[, paste0("I",1:4) ], rater=dat$rater,
                pid=dat$pid )

#-- structure of design matrix
res <- immer:::lpcm_data_prep( dat=dfr1[,-1], weights=NULL, a=NULL)
pars_info <- res$pars_info

#--- define design matrix for multi-faceted Rasch model with only main effects
W <- matrix( 0, nrow=nrow(pars_info), ncol=3+2+2 )
parnames <- rownames(W) <- rownames(pars_info)
colnames(W) <- c( paste0("I",1:3), paste0("Cat",1:2), paste0("R",1:2) )
#+ define item effects
for (ii in c("I1","I2","I3") ){
    ind <- grep( ii, parnames )
    W[ ind, ii ] <- pars_info$cat[ind ]
                }
ind <- grep( "I4", parnames )
W[ ind, c("I1","I2","I3") ] <- -pars_info$cat[ind ]
#+ define step parameters
for (cc in 1:2 ){
    ind <- which( pars_info$cat==cc )
    W[ ind, paste0("Cat",1:cc) ] <- 1
                }
#+ define rater effects
for (ii in c("R1","R2") ){
    ind <- grep( ii, parnames )
    W[ ind, ii ] <- pars_info$cat[ind ]
                }
ind <- grep( "R3", parnames )
W[ ind, c("R1","R2") ] <- -pars_info$cat[ind ]

#--- estimate model with immer_cml
mod1 <- immer::immer_cml( dfr1[,-1], W=W, par_init=rep(0,ncol(W) ) )
summary(mod1)

#--- comparison with estimation in TAM
resp <- dfr1[,-1]
mod2 <- TAM::tam.mml.mfr( resp=dat[,-c(1:2)], facets=dat[, "rater", drop=FALSE ],
            pid=dat$pid, formulaA=~ item + step + rater )
summary(mod2)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Conditional maximum likelihood estimation}
\keyword{summary}
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
