% $Id: addModel.Rd 118 2011-07-22 23:28:53Z Lars $
\name{dea.add}
\alias{dea.add}


\title{Additive DEA model}

\description{ Calculates additive efficiency as sum of input and
output slacks within different DEA models }


\usage{

dea.add(X, Y, RTS="vrs", XREF=NULL, YREF=NULL, 
        FRONT.IDX=NULL, param=NULL, TRANSPOSE=FALSE, LP=FALSE)

}


\arguments{

\item{X}{Inputs of firms to be evaluated, a K x m matrix
  of observations of K firms with m inputs (firm x input).  In case
  \code{TRANSPOSE=TRUE} the input matrix is transposed to input x
  firm. }

\item{Y}{Outputs of firms to be evaluated, a K x n matrix
  of observations of K firms with n outputs (firm x input).  In case
  \code{TRANSPOSE=TRUE} the output matrix is transposed to output x
  firm. }

\item{RTS}{ Text string or a number defining the underlying DEA
technology / returns to scale assumption.

\tabular{rll}{
0 \tab fdh \tab Free disposability hull, no convexity assumption \cr
1 \tab vrs \tab Variable returns to scale, convexity and free disposability \cr
2 \tab drs \tab Decreasing returns to scale, convexity, down-scaling and free disposability \cr
3 \tab crs \tab Constant returns to scale, convexity and free disposability \cr
4 \tab irs \tab Increasing returns to scale, (up-scaling, but not down-scaling), convexity and free disposability \cr
5 \tab add \tab Additivity (scaling up and down, but only with integers), and free disposability
}

}

\item{XREF}{Inputs of the firms determining the technology, defaults
    to \code{X}}

\item{YREF}{Outputs of the firms determining the technology, defaults
    to \code{Y}}

\item{FRONT.IDX}{Index for firms determining the technology}

\item{param}{Possible parameters. At the moment only used for
    RTS="fdh+" to set low and high values for restrictions on lambda;
    see the section details and examples for its use. Future versions
    might also use \code{param} for other purposes.}

\item{TRANSPOSE}{Input and output matrices are treated as firms
      times goods matrices for the default value \code{TRANSPOSE=FALSE}
      corresponding to the standard in R for statistical models.  When
      \code{TRUE} data matrices are transposed to good times firms
      matrices as is normally used in LP formulation of the problem.}

\item{LP}{Only for debugging. If \code{LP=TRUE} then input and output
   for the LP program are written to standard output for each unit.}

}


\value{

\item{sum}{Sum of all slacks for each firm,
     \code{sum=sum(sx)+sum(sy)}.}

\item{slack}{A non-NULL vector of logical variables, \code{TRUE} if
     there is slack for the corresponding firm, and \code{FALSE} if
     the there is no slack, i.e. the sum of slacks is zero.}

\item{sx}{A matrix of input slacks for each firm}

\item{sy}{A matrix of output slack for each firm}

\item{lambda}{The lambdas, i.e. the weights of the peers for each firm}

}


\details{ The sum of the slacks is maximized in a LP formulation of
  the DEA technology. The sum of the slacks can be seen as distance to
  the frontier when you only move parallel to the axes of inputs and
  outputs, i.e. not a usual Euclidean distance, but what is also known
  as an L1 norm.

  Since it is the sum of slacks that is calculated, there is no
  exogenous ORIENTATION in the problem. Rather, there is generally
  both an input and an output direction in the slacks. The model
  considers the input excess and output shortfall simultaneously and
  finds a point on the frontier that is most distant to the point
  being evaluated.

}


\note{ This is neither a Farrell nor a Shephard like efficiency.

The value of the slacks depends on the scaling of the different
inputs and outputs.  Therefore the values are not independent of how
the input and output are measured.}


\author{ Peter Bogetoft and Lars Otto \email{larsot23@gmail.com}}


\source{
   Corresponds to Eqs. 4.34-4.38 in Cooper et al. (2007)
}


\references{

Bogetoft and Otto; \emph{Benchmarking with DEA, SFA, and R}; Springer
   2011

Cooper, Seiford, and Tone; \emph{Data Envelopment Analysis: A
	Comprehensive Text with Models, Applications, References and
	DEA-Solver Software}; Second edition, Springer 2007

}


\examples{
x <- matrix(c(2,3,2,4,6,5,6,8),ncol=1)
y <- matrix(c(1,3,2,3,5,2,3,5),ncol=1)
dea.plot.frontier(x,y,txt=1:dim(x)[1])

sb <- dea.add(x,y,RTS="vrs")
data.frame("sx"=sb$sx,"sy"=sb$sy,"sum"=sb$sum,"slack"=sb$slack)
}

\keyword{models}

