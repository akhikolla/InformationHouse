\name{trcontrol}
\alias{trcontrol}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Set the Control Parameters Used by \code{transport}.
}
\description{
  Set the control parameters for the algorithm used by the function \code{\link{transport}}.
}
\usage{
trcontrol(method = c("networkflow", "revsimplex", "shortsimplex", "primaldual", "aha",
          "shielding", "auction", "auctionbf"), para = list(),
          start = c("auto",          "modrowmin", "nwcorner", "russell"), 
          nscales = 1, scmult = 2,returncoarse = FALSE, a = NULL,
          b = NULL, M = NULL, N = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{method}{
    The algorithm to be used to compute the optimal transference plan. See details for the function \code{\link{transport.pgrid}}.
}
  \item{para}{
    A list of parameters that are specific to the chosen method. See the table on the help page of the function \code{\link{transport}}.
}
  \item{start}{
    If \code{method == "revsimplex"}, the method for computing a starting solution.
}
  \item{nscales, scmult, returncoarse}{
    The parameters for the multiscale versions of certain algorithms. See the help on \code{\link{transport}}.
}
  \item{a,b,M,N}{
    The two objects \code{a} and \code{b} for which the transportation problem is to be solved \emph{or} the sizes \code{M}
    and \code{N} of these objects. Based on the information available here, \code{trcontrol} tries hard to find reasonable values
    for the control parameters of the algorithm not specified directly.
}
}
\details{
  For further details about the parameters of the individual algorithms see the help page for \code{\link{transport}}.
}
\value{
  A list with components \code{method}, \code{para}, \code{start}, \code{nscales}, \code{scmult}, \code{returncoarse} as
  entered or adapted/computed based on the arguments \code{method}, \code{a}, \code{b}, \code{M}, \code{N}.
}
\author{
  Dominic Schuhmacher \email{dschuhm1@uni-goettingen.de}
}
\note{
  This function is typically only called by the user to check what the parameter settings used
  by the function \code{\link{transport}} are for a given problem.
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{transport}}
}
\keyword{control}
\keyword{transport parameters}