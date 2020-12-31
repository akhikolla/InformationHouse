% $Id: stoned.Rd 221 2020-06-09 21:18:37Z lao $
\name{stoned}
\alias{stoned}

\title{Convex nonparametric least squares}

\description{
    Convex nonparametric least squares here for convex (Cost) function 
    function or concave (Production) function with multiplicative or 
    additive error term. the StoNED estimator combines the
    axiomatic and non-parametric frontier (the DEA aspect) with
    a stochastic noise term (the SFA aspect)}

\usage{
stoned(X, Y, RTS = "vrs", COST = 0, MULT = 0, METHOD = "MM")
}

\arguments{
  \item{X}{Inputs (right hand side) of firms to be evaluated, a K x m matrix
  of observations of K firms with m inputs (firm x input).}

  \item{Y}{Output or cost (left hand side) of firms to be evaluated, a K x 1 matrix
  of observations of K firms with 1 output or cost (firm x input).}

  \item{RTS}{RTS determines returns to scale assumption: RTS="vrs", 
  "drs", "crs"  and "irs" are possible for constant or variable returns
  to scale; see \code{\link{dea}} for a verbal description and numberring scheme.}

  \item{COST}{COST specifies whether a cost function needs is 
  estimated (COST=1) or a production function (COST=0).}

  \item{MULT}{MULT determines if multiplicative (MULT=1) or 
  additive (MULT=0) model is estimated.}

  \item{METHOD}{METHOD specifies the way efficiency is estimated: 
  MM for Method of Momente and PSL for pseudo likelihood estimation.}
}

\details{Convex nonparametric least squares here for convex (cost)
function with multiplicative error term: Y=b*X*exp(e) or additive
error term: Y=b*X + e.}

\value{The results are returned in a list with the components:
    \item{residualNorm}{Norm of residual}
    \item{solutionNorm}{Norm of solution}
    \item{error}{Is there an error in the solution?}
    %%\item{type}{Type of method used for finding the solution in the QP problem}
    \item{coef}{beta_matrix, estimated coefficients as a Kxm matrix; if there 
        is an intercept the first collumn is the intercept, and the matrix is Kx(1+m)}
    \item{residuals}{Residuals}
    \item{fit}{Fitted values}
    \item{eff}{Efficinecy score}
    \item{front}{Points on the frontier}
    \item{sigma_u}{sigma_u}
}

\references{
Kuosmanen and Kortelainen, "Stochastic non-smooth envelopment of data: 
semi-parametric frontier estimation subject to shape constraints", 
\emph{Journal of Productivity Analysis} 2012
}

\author{
Stefan Seifert \email{s.seifert@ilr.uni-bonn.de} and Lars Otto \email{larsot23@gmail.com}
}

\note{
Convex nonparametric least squares here for convex (Cost) function
with multiplicative error term:
\code{Y=b*X*exp(e)} or additive error term: \code{Y=b*X + e}.

The intercept is absent for the constant returns to scale assumption; all 
other technology assumptions do have an in tercept.

Note that the method \code{stoned} is a rather slow method and probably only 
works in a reasonable time for less than 3-400 units.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%%}

\examples{
#### Example: Single Input Production Function
n=10

x1 <- runif(n,10,20)
v <- rnorm(n,0,0.01)
u <- abs(rnorm(n,0,0.04))

y <- (x1^0.8)*exp(-u)*exp(v)

sol_MM <- stoned(x1, y)
sol_PSL <- stoned(x1, y, METHOD="PSL")

plot(x1,y)
curve(x^0.8, add=TRUE)
points(x1,sol_MM$front, col="red")
points(x1,sol_PSL$front, col="blue", pch=16, cex=.6)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{efficiency}% use one of  RShowDoc("KEYWORDS")
\keyword{models}% __ONLY ONE__ keyword per line
