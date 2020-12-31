% $Id: cost.opt.Rd 229 2020-07-04 13:39:18Z lao $
\name{cost.opt}
\alias{cost.opt}
\alias{print.cost.opt}
\alias{summary.cost.opt}

\alias{revenue.opt}
\alias{print.revenue.opt}
\alias{summary.revenue.opt}

\alias{profit.opt}
\alias{print.profit.opt}
\alias{summary.profit.opt}


\title{DEA optimal cost, revenue, and profit}

\description{Estimates the input and/or output vector(s) that minimize
cost, maximize revenue or maximize profit in the context of a DEA
technology}


\usage{
cost.opt(XREF, YREF, W, YOBS=NULL, RTS="vrs", param=NULL,
         TRANSPOSE=FALSE, LP=FALSE, CONTROL=NULL, LPK = NULL)  

revenue.opt(XREF, YREF, P, XOBS=NULL, RTS="vrs",  param=NULL,
            TRANSPOSE = FALSE, LP = FALSE, CONTROL=NULL, LPK = NULL)

profit.opt(XREF, YREF, W, P, RTS = "vrs",  param=NULL,
           TRANSPOSE = FALSE, LP = FALSE, CONTROL=NULL, LPK = NULL)

}

\arguments{

\item{}{Input and output matrices are in the same form as for the
  method \code{\link{dea}}.}

\item{XREF}{Input of the firms defining the technology, a K x m matrix
  of observations of K firms with m inputs (firm x input).  In case
  \code{TRANSPOSE=TRUE} the input matrix is transposed as input x
  firm. }

\item{YREF}{output of the firms defining the technology, a K x n
  matrix of observations of K firms with n outputs (firm x input).  In
  case \code{TRANSPOSE=TRUE} the output matrix is transposed as output
  x firm. }

\item{W}{Input prices as a matrix. Either same prices for all firms or
  individual prices for all firms; i.e. either a 1 x m or a K x m
  matrix for K firms and m inputs }

\item{P}{Output prices as a matrix. Either same prices for all firms
  or individual prices for all firms; i.e. either a 1 x n or K x n
  matrix for K firms and n outputs }

\item{XOBS}{The input for which an optimal, revenue maximizing, output
  vector is to be calculated. Defaults is \code{XREF}.  Same form as
  \code{XREF}}

\item{YOBS}{The output for which an optimal, cost minimizing input
  vector is to be calculated. Defaults is \code{YREF}.  Same form as
  \code{YREF}}

\item{RTS}{

  A text string or a number defining the underlying DEA technology /
  returns to scale assumption.

\tabular{rll}{
0 \tab fdh \tab Free disposability hull, no convexity assumption \cr
1 \tab vrs \tab Variable returns to scale, convexity and free disposability \cr
2 \tab drs \tab Decreasing returns to scale, convexity, down-scaling and free disposability \cr
3 \tab crs \tab Constant returns to scale, convexity and free disposability \cr
4 \tab irs \tab Increasing returns to scale, (up-scaling, but not down-scaling), convexity and free disposability \cr
5 \tab add \tab Additivity (scaling up and down, but only with integers), and free disposability \cr
6 \tab fdh+ \tab A combination of free disposability and restricted
 or local constant return to scale
}

}

\item{param}{Possible parameters. At the moment only used for
    RTS="fdh+" to set low and high values for restrictions on lambda;
    see the section details and examples in \code{\link{dea}} for its
    use. Future versions might also use \code{param} for other
    purposes.}


\item{TRANSPOSE}{Input and output matrices are treated as firms times
    goods for the default value \code{TRANSPOSE=FALSE} corresponding
    to the standard in R for statistical models.  When \code{TRUE}
    data matrices, quantities and prices, are transposed to goods times
    firms matrices.}


\item{LP}{Only for debugging. If \code{LP=TRUE} then input and output
    for the LP program are written to standard output for each unit.}

\item{CONTROL}{Possible controls to \pkg{lpSolveAPI}, see the
    documentation for that package. For examples of use see the
    function \code{\link{dea}}.}

\item{LPK}{When \code{LPK=k} then a mps file is written for firm
    \code{k}; it can be used as input to an alternative LP solver
    to check the results.}

}



\details{ The LP optimization problem is formulated in Bogetoft and
  Otto (2011, pp 35 and 102) and is solved by the LP method in the
  package \pkg{lpSolveAPI}.

  The methods \code{print} and \code{summary} are working for
  \code{cost.opt}, \code{revenue.opt}, and \code{profit.opt}

}



\value{ The values returned are the optimal input, and/or optimal
        output. When saved in an object the following components are
        available:

\item{xopt}{The optimal input, returned as a matrix by
     \code{cost.opt} and \code{profit.cost}.}

\item{yopt}{The optimal output, returned as a matrix by
     \code{revenue.opt} and \code{profit.cost}.}

\item{cost}{The optimal/minimal cost.}

\item{revenue}{The optimal/maximal revenue}

\item{profit}{The optimal/maximal profit}

\item{lambda}{The peer weights that determines the technology, a
    matrix. Each row is the lambdas for the firm corresponding to that
    row; for the vrs technology the rows sum to 1. A column shows for
    a given firm how other firms are compared to this firm; i.e. peers
    are firms with a positive element in their columns.}

}


\references{
Bogetoft and Otto; \emph{Benchmarking with DEA, SFA, and R}; Springer
   2011
}


\author{ Peter Bogetoft and Lars Otto \email{larsot23@gmail.com}}


\note{

The index for peer units can be returned by the method \code{peers}
and the weights are returned in \code{lambda}.  Note that the peers
now are the firms for the optimal input and/or output allocation, not
just the technical efficient firms.

}


\seealso{

Paul W. Wilson (2008), \dQuote{FEAR 1.0: A Software Package for Frontier
Efficiency Analysis with R,} \emph{Socio-Economic Planning Sciences}
42, 247--254

}

\examples{

x <- matrix(c(2,12, 2,8, 5,5, 10,4, 10,6, 3,13), ncol=2, byrow=TRUE)
y <- matrix(1,nrow=dim(x)[1],ncol=1)
w <- matrix(c(1.5, 1),ncol=2)

txt <- LETTERS[1:dim(x)[1]]
dea.plot(x[,1],x[,2], ORIENTATION="in",  cex=1.25)
text(x[,1],x[,2],txt,adj=c(-.7,-.2),cex=1.25)

# technical efficiency
te <- dea(x,y,RTS="vrs")
xopt <- cost.opt(x,y,w,RTS=1)
cobs <- x \%*\% t(w)
copt <- xopt$x \%*\% t(w)
# cost efficiency
ce <- copt/cobs
# allocaltive efficiency
ae <- ce/te$eff
data.frame("ce"=ce,"te"=te$eff,"ae"=ae)
print(cbind("ce"=c(ce),"te"=te$eff,"ae"=c(ae)),digits=2)

# isocost line in the technology plot
abline(a=copt[1]/w[2], b=-w[1]/w[2], lty="dashed")
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
\keyword{efficiency}
\keyword{DEA}
