% $Id: dea.direct.Rd 129 2014-06-16 13:31:01Z b002961 $
\name{dea.direct}
\alias{dea.direct}

\title{ Directional efficiency }

\description{ Directional efficiency rescaled to an interpretation a la
   Farrell efficiency and the corresponding peer importance
   (lambda). }


\usage{
dea.direct(X, Y, DIRECT, RTS = "vrs", ORIENTATION = "in", 
          XREF = NULL, YREF = NULL, FRONT.IDX = NULL, 
          SLACK = FALSE, param=NULL, TRANSPOSE = FALSE)
}


\arguments{

  \item{X}{Inputs of firms to be evaluated, a K x m matrix of
      observations of K firms with m inputs (firm x input) }

  \item{Y}{Outputs of firms to be evaluated, a K x n matrix of
      observations of K firms with n outputs (firm x input). }

  \item{DIRECT}{ Directional efficiency, \code{DIRECT} is either a
      scalar, an array, or a matrix with non-negative elements.

      If the argument is a scalar, the direction is (1,1,...,1) times
      the scalar; the value of the efficiency depends on the scalar as
      well as on the unit of measurements.

      If the argument an array, this is used for the direction for
      every firm; the length of the array must correspond to the
      number of inputs and/or outputs depending on the
      \code{ORIENTATION}.

      If the argument is a matrix then different directions are used
      for each firm. The dimensions depends on the \code{ORIENTATION}
      (and \code{TRANSPOSE}), the number of firms must correspond to
      the number of firms in \code{X} and \code{Y}.

      \code{DIRECT} must not be used in connection with
      \code{DIRECTION="graph"}.
      
  }

  \item{RTS}{Text string or a number defining the underlying DEA
      technology / returns to scale assumption.

\tabular{rll}{
0 \tab fdh \tab Free disposability hull, no convexity assumption \cr
1 \tab vrs \tab Variable returns to scale, convexity and free disposability \cr
2 \tab drs \tab Decreasing returns to scale (down-scaling, but not up-scaling), convexity, and free disposability \cr
3 \tab crs \tab Constant returns to scale, convexity and free disposability \cr
4 \tab irs \tab Increasing returns to scale (up-scaling, but not down-scaling), convexity, and free disposability \cr
6 \tab add \tab Additivity (scaling up and down, but only with integers), and free disposability\cr
7 \tab fdh+ \tab A combination of free disposability and restricted
 or local constant return to scale
}

}

  \item{ORIENTATION}{Input efficiency "in" (1), output
      efficiency "out" (2), and graph efficiency "graph" (3). For use
      with \code{DIRECT}, an additional option is "in-out" (0). }

  \item{XREF}{Inputs of the firms determining the technology, defaults
      to \code{X}. }

  \item{YREF}{Outputs of the firms determining the technology,
      defaults to \code{Y}. }

  \item{FRONT.IDX}{Index for firms determining the technology. }

  \item{SLACK}{See \code{\link{dea}} and \code{\link{slack}}. }

  \item{param}{Possible parameters. At the moment only used for
    RTS="fdh+" to set low and high values for restrictions on lambda;
    see the section details and examples in \code{\link{dea}} for its
    use. Future versions might also use \code{param} for other
    purposes.}

  \item{TRANSPOSE}{ see \code{\link{dea}}}

}


\details{

  When the argument \code{DIRECT=d} is used then component
  \code{objval} of the returned object for input orientation is the
  maximum value of e where for input orientation \eqn{x-e d}, and for
  output orientation \eqn{y+e d} are in the generated technology
  set. The returned component \code{eff} is for input \eqn{1-e d/X}
  and for output \eqn{1+e d /Y} to make the interpretation as for a
  Farrell efficiency.  Note that when the direction is not
  proportional to \code{X} or \code{Y} the returned \code{eff} are
  different for different inputs or outputs and \code{eff} is a matrix
  and not just an array.  The directional efficiency can be restricted
  to inputs (\code{ORIENTATION="in"}), restricted to outputs
  (\code{ORIENTATION="out"}), or both include inputs and output
  directions (\code{ORIENTATION="in-out"}). Directional efficiency is
  discussed on pages 31--35 and 121--127 in Bogetoft and Otto (2011).

  The Farrell efficiency interpretation is the ratio by which a firm
  can proportionally reduce all inputs (or expand all outputs) without
  producing less outputs (using more inputs). The directional
  efficiecies have the same interpretation expect that the direction is
  not proportional to the inputs (or outputs) and therefore the
  different inputs may have different reduction ratios, the
  efficiency is an array and not just a number.

}



\value{The results are returned in a Farrell object with the following
      components.  The method \code{slack} only returns the three
      components in the list relevant for slacks.

\item{eff}{The Farrell efficiencies.  Note that the the efficiencies
    are calculated to have the same interpretations as Farrell
    efficiencies. \code{eff} is a matrix if there are more than 1
    good.}

\item{lambda}{The lambdas, i.e. the weight of the peers, for each firm}

\item{objval}{The objective value as returned from the LP program; the
     \code{objval} are excess values in DIRECT units of measurement.}

\item{RTS}{The return to scale assumption as in the option \code{RTS}
     in the call}

\item{ORIENTATION}{The efficiency orientation as in the call}

\item{TRANSPOSE}{As in the call}

\item{slack}{A vector with sums of the slacks for each firm. Only
     calculated in dea when option \code{SLACK=TRUE}}

\item{sx}{A matrix for input slacks for each firm, only calculated if
     the option \code{SLACK} is \code{TRUE} or returned from the
     method \code{slack}}

\item{sy}{A matrix for output slack, see \code{sx}}
}


\references{ Directional efficiency is discussed on pages 31--35 and
121--127 in Bogetoft and Otto (2011).

Peter Bogetoft and Lars Otto; \emph{Benchmarking with DEA, SFA, and
R}; Springer 2011 }


\author{ Peter Bogetoft and Lars Otto \email{larsot23@gmail.com} }


\note{ To handle fixed, non-discretionary inputs, one can let it
  appear as negative output in an input-based mode, and reversely for
  fixed, non-discretionary outputs. Fixed inputs (outputs) can also be
  handled by directional efficiency; set the direction, the argument
  \code{DIRECT}, equal to the variable, discretionary inputs (outputs)
  and 0 for the fixed inputs (outputs).

  When the the argument \code{DIRECT=X} is used the then the returned
  efficiency is equal to 1 minus the Farrell efficiency for input
  orientation and equal to the Farrell efficiency minus 1 for output
  orientation.

}


\seealso{
\code{\link{dea}}
}


\examples{
# Directional efficiency
x <- matrix(c(2,5 , 1,2 , 2,2 , 3,2 , 3,1 , 4,1), ncol=2,byrow=TRUE)
y <- matrix(1,nrow=dim(x)[1])
dea.plot.isoquant(x[,1], x[,2],txt=1:dim(x)[1])

E <- dea(x,y)
z <- c(1,1)
e <- dea.direct(x,y,DIRECT=z)
data.frame(Farrell=E$eff, Perform=e$eff, objval=e$objval)
# The direction
arrows(x[,1], x[,2], (x-z)[,1], (x-z)[,2], lty="dashed")
# The efficiency (e$objval) along the direction
segments(x[,1], x[,2], (x-e$objval*z)[,1], (x-e$objval*z)[,2], lwd=2)



# Different directions
x1 <- c(.5, 1, 2, 4, 3, 1)
x2 <- c(4,  2, 1,.5, 2, 4)
x <- cbind(x1,x2)
y <- matrix(1,nrow=dim(x)[1])
dir1 <- c(1,.25)
dir2 <- c(.25, 4)
dir3 <- c(1,4)
e <- dea(x,y)
e1 <- dea.direct(x,y,DIRECT=dir1)
e2 <- dea.direct(x,y,DIRECT=dir2)
e3 <- dea.direct(x,y,DIRECT=dir3)
data.frame(e=eff(e),e1=e1$eff,e2=e2$eff,e3=e3$eff)[6,]

# Technology and directions for all firms
dea.plot.isoquant(x[,1], x[,2],txt=1:dim(x)[1])
arrows(x[,1], x[,2],  x[,1]-dir1[1], x[,2]-dir1[2],lty="dashed")
segments(x[,1], x[,2],  
    x[,1]-e1$objval*dir1[1], x[,2]-e1$objval*dir1[2],lwd=2)
# slack for direction 1
dsl1 <- slack(x,y,e1)
cbind(E=e$eff,e1$eff,dsl1$sx,dsl1$sy, sum=dsl1$sum)



# Technology and directions for firm 6, 
# Figure 2.6 page 32 in Bogetoft & Otto (2011)
dea.plot.isoquant(x1,x2,lwd=1.5, txt=TRUE)
arrows(x[6,1], x[6,2],  x[6,1]-dir1[1], x[6,2]-dir1[2],lty="dashed")
arrows(x[6,1], x[6,2],  x[6,1]-dir2[1], x[6,2]-dir2[2],lty="dashed")
arrows(x[6,1], x[6,2],  x[6,1]-dir3[1], x[6,2]-dir3[2],lty="dashed")
segments(x[6,1], x[6,2],  
    x[6,1]-e1$objval[6]*dir1[1], x[6,2]-e1$objval[6]*dir1[2],lwd=2)
segments(x[6,1], x[6,2],  
    x[6,1]-e2$objval[6]*dir2[1], x[6,2]-e2$objval[6]*dir2[2],lwd=2)
segments(x[6,1], x[6,2],  
    x[6,1]-e3$objval[6]*dir3[1], x[6,2]-e3$objval[6]*dir3[2],lwd=2)
}


\keyword{ models }
\keyword{ efficiency }
