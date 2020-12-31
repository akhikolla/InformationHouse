% $Id: dea.dual.Rd 118 2011-07-22 23:28:53Z Lars $
\name{dea.dual}
\alias{dea.dual}
% \alias{dea.dual.mat}
% \encoding{UTF-8}

\title{ Dual DEA models and assurance regions }

\description{ Solution of dual DEA models, possibly with partial value
information given as restrictions on the ratios (assurance regions) }

\usage{
dea.dual(X, Y, RTS = "vrs", ORIENTATION = "in", 
         XREF = NULL, YREF = NULL, 
         FRONT.IDX = NULL, DUAL = NULL, DIRECT=NULL,
         TRANSPOSE = FALSE, LP = FALSE, CONTROL=NULL, LPK=NULL)
}


% dea.dual.mat(X,Y, RTS="vrs", ORIENTATION="in", XREF=NULL,YREF=NULL,
%         FRONT.IDX=NULL, DUAL=NULL, # SLACK=FALSE, delta=0, DELTA=FALSE, 
%         TRANSPOSE=FALSE, LP=FALSE, ...)


\arguments{

\item{X}{Inputs of firms to be evaluated, a K x m matrix of
  observations of K firms with m inputs (firm x input).  In case
  \code{TRANSPOSE=TRUE} the input matrix is transposed to input x
  firm. }

\item{Y}{Outputs of firms to be evaluated, a K x n matrix of
  observations of K firms with n outputs (firm x input).  In case
  \code{TRANSPOSE=TRUE} the output matrix is transposed to output x
  firm. }

\item{RTS}{ A text string or a number defining the underlying DEA
 technology / returns to scale assumption.

\tabular{rll}{
% 0 \tab fdh \tab  \emph{Not available} for \code{dea.dual}.\cr
1 \tab vrs \tab Variable returns to scale, convexity and free disposability \cr
2 \tab drs \tab Decreasing returns to scale, convexity, down-scaling and free disposability \cr
3 \tab crs \tab Constant returns to scale, convexity and free disposability \cr
4 \tab irs \tab Increasing returns to scale, (up-scaling, but not down-scaling), convexity and free disposability.
}

}

\item{ORIENTATION}{Input efficiency "in" (1), output efficiency "out"
      (2), and graph efficiency "graph" (3) (not yet implemented). For
      use with \code{DIRECT} an additional option is "in-out" (0). In
      this case, "graph" is not feasible}

\item{XREF}{Input of the firms determining the technology, defaults to
\code{X}}

\item{YREF}{Output of the firms determining the technology, defaults
to \code{Y}}

\item{FRONT.IDX}{Index for firms determining the technology}


\item{DUAL}{

 Matrix of order \dQuote{number of input plus number of outputs minus
 2} times 2.  The first column is the lower bound and the second column is
 the upper bound for the restrictions on the multiplier ratios. The ratios are relative to the first input and the first
 output, respectively. This implyies that there is no restriction for neither the first
 input nor the first output so that the number of restrictions is two less
 than the total number of inputs and outputs.

}

\item{DIRECT}{Directional efficiency, \code{DIRECT} is either a
     scalar, an array, or a matrix with non-negative elements.

     NB \emph{Not yet implemented}}


\item{TRANSPOSE}{Input and output matrices are treated as firms
      times goods for the default value \code{TRANSPOSE=FALSE}
      corresponding to the standard in R for statistical models.  When
      \code{TRUE} data matrices shall be transposed to good times firms
      matrices as is normally used in LP formulation of the problem.}


 \item{LP}{Only for debugging. If \code{LP=TRUE} then input and
      output for the LP program are written to standard output for
      each unit.}

\item{CONTROL}{Possible controls to lpSolveAPI, see the documentation
  for that package.}

 \item{LPK}{When \code{LPK=k} then a mps file is written for firm
   \code{k}; it can be used as input to an alternative LP solver just
   to check the our results.}

}

\details{ Solved as an LP program using the package lpSolveAPI. The
method \code{dea.dual.dea} calls the method \code{dea} with the option
\code{DUAL=TRUE}.  }



\value{

\item{eff}{ The efficiencies }

\item{objval}{ The objective value as returned from the LP problem, normally the same as eff }

\item{RTS}{The return to scale assumption as in the option \code{RTS}
     in the call}

\item{ORIENTATION}{ The efficiency orientation as in the call }

\item{TRANSPOSE}{ As in the call }

\item{u}{Dual values, prices, for inputs }

\item{v}{Dual values, prices, for outputs }

\item{gamma}{ The values of gamma, the shadow price(s) for returns to
scale restriction }

\item{sol}{Solution of all variables as one component, sol=c(u,v,gamma).}

}


\references{ Bogetoft and Otto; \emph{Benchmarking with DEA, SFA, and
 R}; Springer 2011. Sect.  5.10: Partial value information }


\author{ Peter Bogetoft and Lars Otto \email{larsot23@gmail.com}}


\note{ 

Note that the dual values are not unique for extreme points in the
technology set. In this case the value of the calculated dual variable can
depend on the order of the complete efficient firms.

}



\seealso{
\code{\link{dea}}
}

\examples{

x <- matrix(c(2,5 , 1,2 , 2,2 , 3,2 , 3,1 , 4,1), ncol=2,byrow=TRUE)
y <- matrix(1,nrow=dim(x)[1])
dea.plot.isoquant(x[,1],x[,2],txt=1:dim(x)[1])
segments(0,0, x[,1], x[,2], lty="dotted")


e <- dea(x,y,RTS="crs",SLACK=TRUE)
ed <- dea.dual(x,y,RTS="crs")
print(cbind("e"=e$eff,"ed"=ed$eff, peers(e), lambda(e), 
            e$sx, e$sy, ed$u, ed$v), digits=3)

dual <- matrix(c(.5, 2.5), nrow=dim(x)[2]+dim(y)[2]-2, ncol=2, byrow=TRUE)
er <- dea.dual(x,y,RTS="crs", DUAL=dual)
print(cbind("e"=e$eff,"ar"=er$eff, lambda(e), e$sx, e$sy, er$u, 
            "ratio"=er$u[,2]/er$u[,1],er$v),digits=3)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
