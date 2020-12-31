% $Id: dea.Rd 229 2020-07-04 13:39:18Z lao $
\name{dea}
\alias{dea}
% \alias{dea.mat}
% \alias{efficiencies}
% \alias{eff}
% \alias{efficiencies.Farrell}
% \alias{eff.Farrell}
\alias{print.Farrell}
\alias{summary.Farrell}
% \encoding{UTF-8}

\title{DEA efficiency}

\description{Estimates a DEA frontier and calculates efficiency
measures a la Farrell.}

\usage{
dea(X, Y, RTS="vrs", ORIENTATION="in", XREF=NULL, YREF=NULL,
    FRONT.IDX=NULL, SLACK=FALSE, DUAL=FALSE, DIRECT=NULL, param=NULL,
    TRANSPOSE=FALSE, FAST=FALSE, LP=FALSE, CONTROL=NULL, LPK=NULL)

\method{print}{Farrell}(x, digits=4, ...) 
\method{summary}{Farrell}(object, digits=4, ...)
}



%- maybe also 'usage' for other objects documented here.
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

1 \tab vrs \tab Variable returns to scale, convexity and free
disposability \cr

2 \tab drs \tab Decreasing returns to scale, convexity, down-scaling
and free disposability \cr

3 \tab crs \tab Constant returns to scale, convexity and free
disposability \cr

4 \tab irs \tab Increasing returns to scale, (up-scaling, but not
down-scaling), convexity and free disposability \cr

5 \tab irs2 \tab Increasing returns to scale (up-scaling, but not
down-scaling), additivity, and free disposability \cr

6 \tab add \tab Additivity (scaling up and down, but only with
integers), and free disposability; also known af replicability and 
free disposability, the free disposability and replicability hull 
(frh) -- no convexity assumption \cr

7 \tab fdh+ \tab A combination of free disposability and restricted
 or local constant return to scale \cr
 
%8 \tab fdh++ \tab Should be the same as fdh+, but without use of Special 
%Ordered Set in the solution of the LP problem; unfortunately this is not
%always the case, and the problem is with fdh++ \cr

%9 \tab fdh0 \tab Same as fdh, but always solved as a LP problem 

10 \tab vrs+ \tab As vrs, but with restrictions on the individual lambdas 
via \code{param}

}

}

\item{ORIENTATION}{Input efficiency "in" (1), output efficiency "out"
    (2), and graph efficiency "graph" (3). For use with \code{DIRECT},
    an additional option is "in-out" (0).}

\item{XREF}{Inputs of the firms determining the technology, defaults
    to \code{X}}

\item{YREF}{Outputs of the firms determining the technology, defaults
    to \code{Y}}

\item{FRONT.IDX}{Index for firms determining the technology}

\item{SLACK}{Calculate slack in a phase II calculation by an intern 
    call of the function \code{\link{slack}}. Note that the precision
    for calculating slacks for orientation graph is low. }
     
\item{DUAL}{Calculate dual variables, i.e. shadow prices; not
    calculated for orientation graph as that is not an LP problem.}

\item{DIRECT}{Directional efficiency, \code{DIRECT} is either a
    scalar, an array, or a matrix with non-negative elements.

    If the argument is a scalar, the direction is (1,1,...,1) times
    the scalar; the value of the efficiency depends on the scalar as
    well as on the unit of measurements.

    If the argument is an array, this is used for the direction for every
    firm; the length of the array must correspond to the number of
    inputs and/or outputs depending on the \code{ORIENTATION}.

    If the argument is a matrix then different directions are used for
    each firm. The dimensions depends on the \code{ORIENTATION} (and
    \code{TRANSPOSE}), the number of firms must correspond to the
    number of firms in \code{X} and \code{Y}.

    \code{DIRECT} must not be used in connection with
    \code{ORIENTATION="graph"}.

   }

\item{param}{Possible parameters. At the moment only used for
    RTS="fdh+" to set low and high values for restrictions on lambda;
    see the section details and examples for its use. Future versions
    might also use \code{param} for other purposes.}

\item{TRANSPOSE}{Input and output matrices are treated as firms
      times goods matrices for the default value \code{TRANSPOSE=FALSE}
      corresponding to the standard in R for statistical models.  When
      \code{TRUE} data matrices are transposed to good times firms
      matrices as is normally used in LP formulation of the problem.}

\item{LP}{Only for debugging. If \code{LP=TRUE} then input and
     output for the LP program are written to standard output for each
     unit.}

\item{FAST}{Only calculate efficiencies and just return them as a
     vector, i.e. no lambda or other output. The return when using
     FAST cannot be used as input for \code{slack} and \code{peers}.}

\item{CONTROL}{Possible controls to \pkg{lpSolveAPI}, see the
      documentation for that package.}

\item{\dots}{Optional parameters for the print and summary methods.}

\item{object, x}{An object of class Farrell (returned by the
    function \code{dea}) -- R code uses \sQuote{\code{object}} and
    \sQuote{\code{x}} alternating for generic methods.}

\item{digits}{digits in printed output, handled by format in print.}

\item{LPK}{when \code{LPK=k} then a mps file is written for firm
   \code{k}; it can be used as input to an alternative LP solver
   to check the results.}

}


\details{ 

The return from \code{dea} and \code{sdea} is an object of class
Farrell.  The efficiency in \code{dea} is calculated by the LP method
in the package \pkg{lpSolveAPI}.  Slacks can be calculated either in
the call of \code{dea} using the option \code{SLACK=TRUE} or in a
following call to the function \code{\link{slack}}.

The directional efficiency when the argument \code{DIRECT} is used,
depends on the unit of measurement and is not restricted to be less
than 1 (or greater than 1 for output efficiency) and is therefore
completely different from the Farrell efficiency.

The crs factor in RTS="fdh+" that sets the lower and upper bound can
be changed by the argument \code{param} that will set the lower and
upper bound to 1-param and 1+param; the default value is
param=.15. The value must be greater than or equal to 0 and strictly
less than 1. A value of 0 corresponds to RTS="fdh".  To get an
asymmetric interval set param to a 2 dimentional array with values for
the low and high end for interval, for instance
\code{param=c(.8,1.15)}.  The FDH+ technology set is described in
Bogetoft and Otto (2011) pages 73--74.

The technology RTS="vrs+" uses the parameter \code{param} to set 
restrictions on lambda, the convexity parameters. The elements of param 
are \code{param=(low, high, sum_low, sum_high)} where "low" and "high"
are restrictions on the individual lambda and "sum_low" and "sum_high" are 
restrictions on the sum of lambdas. The individual lambda must be in the 
interval from low to high or be zero. With one parameter the restrictions 
set are \code{(param, 1+1-(param),1,1)}, with two parameters 
\code{(param[1], param[2],1,1)}, and with four parameters 
\code{(param[1], param[2],param[3], param[4])}. The resulting technology set is 
not necessarily convex.

% The method \code{dea.mat} uses \code{lpSolve} as the LP engine instead
% of the \code{lpSolveAOI} LP engine. The result is that \code{dea.mat}
% is much slower than \code{dea} except for data with very small number
% of firms, like less than 10 firms.

The graph orientated efficiency is calculated by bisection between
feasible and infeasible values of G. The precision in the result is
less than for the other orientations.

When the argument \code{DIRECT=d} is used then the returned value
\code{e} for input orientation is the exces input measured in \code{d}
units of measurements, i.e. \eqn{x-e d}, and for output orientation
\eqn{y+e d}.  The directional efficency can be restricted to inputs
(\code{ORIENTAION="in"}), restricted to outputs
(\code{ORIENTAION="out"}), or both include inputs and output
directions (\code{ORIENTAION="in-out"}). Dirctional efficiency is
discussed on pages 31--35 and 121--127 in Bogetoft and Otto (2011).

}


\value{ The results are returned in a Farrell object with the
    following components.  The last three components in the list are
    only part of the object when \code{SLACK=TRUE}.

\item{eff}{The efficiencies.  Note when DIRECT is used then the
    efficencies are not Farrell efficiencies but rather exces values
    in DIRECT units of measurement }

\item{lambda}{The lambdas, i.e. the weight of the peers, for each firm}

\item{objval}{The objective value as returned from the LP program;
    normally the same as eff, but for \code{slack} it is the the sum
    of the slacks. }

\item{RTS}{The return to scale assumption as in the option \code{RTS}
     in the call}

\item{ORIENTATION}{The efficiency orientation as in the call}

\item{TRANSPOSE}{As in the call}

\item{slack}{A logical vector where the component for a firm is
     \code{TRUE} if the sums of slacks for the corresponding firm is
     positive. Only calculated in dea when option \code{SLACK=TRUE}}

\item{sum}{A vector with sums of the slacks for each firm. Only
     calculated in dea when option \code{SLACK=TRUE}}

\item{sx}{A matrix for input slacks for each firm, only calculated if
     the option \code{SLACK} is \code{TRUE} or returned from the
     method \code{slack}}

\item{sy}{A matrix for output slack, see \code{sx}}

\item{ux}{Dual variable for input, only calculated if \code{DUAL} is
     \code{TRUE}.}

\item{vy}{Dual variable for output, only calculated if \code{DUAL} is
    \code{TRUE}.}

}


\references{ Peter Bogetoft and Lars Otto; \emph{Benchmarking with
DEA, SFA, and R}; Springer 2011 }

\author{ Peter Bogetoft and Lars Otto \email{larsot23@gmail.com}}


\note{ 

The arguments \code{X}, \code{Y}, \code{XREF}, and \code{YREF} are
supposed to be matrices or numerical data frames that in the function
will be converted to matrices. When subsetting a matrix or data frame
to just one column then the class of the resulting object/variable is
no longer a matrix or a data frame, but just a numeric (array,
vector). Therefore, in this case a numeric input that is not a matrix
nor a data frame is transformed to a 1 column matrix, and here the use
of the argument \code{TRANSPOSE=TRUE} gives an error.

The dual values are not unique for extreme points (firms on the
boundary with an efficiency of 1) and therefore the calculated dual
values for these firms can depend on the order of firms in the
reference technology. The same lack of uniqueness also make the peers
for some firms depend on the order of firms in the reference
technology.

To calucalte slack use the argument \code{SLACK=TRUE} or use the
function \code{\link{slack}} directly.

When there is slack, and slack is not taken into consideration, then
the peers for a firm with slack might depend on the order of firms in
the data set; this is a property of the LP algorithm used to solve the
problem.

To handle fixed, non-discretionary inputs, one can let it appear as
negative output in an input-based mode, and reversely for fixed,
non-discretionary outputs. Fixed inputs (outputs) can also be handled
by directional efficiency; set the direction, the argument
\code{DIRECT}, equal to the variable, discretionary inputs (outputs)
and 0 for the fixed inputs (outputs).

When the the argument \code{DIRECT=X} is used the then the returned
effiency is equal to 1 minus the Farrell efficiency for input
orientation and to the Farrell effiency minus 1 for output
orientation.

To use matrices \code{X} and \code{Y} prepared for the methods in the
package \pkg{FEAR} (Wilson 2008) set the options
\code{TRANSPOSE=TRUE}; for consistency with \pkg{FEAR} the options
\code{RTS} and \code{ORIENTATION} also accepts numbers as in
\pkg{FEAR}.

The tolerance that lambda is zero or one is 1e-7, the default value of
'epsint' in the package \code{lpSolveAPI}, i.e. values closer than
1e-7 from zero or one are set to respective integer value. The 'epsint'
is the tolerance that is used to determine whether a floating-point
number is in fact an in teger. The same tolerance is used for
efficiency value near one.

Some scaling is done in the function, but this does not always work 
satisfactory, i.e. sometime, a solution cannot always be found -- the 
program prints a warning and the efficiency for the firm is set to NA. 
Often this is due to a bad scaling of the data. Either the user can try 
a different scaling of data when calling the function or one can use the
option \code{CONTROL} to try a different scaling I the program. For 
instance one can try \code{CONTROL=list(scaling=c("geometric", "equilibrate")}
or \code{CONTROL=list(scaling=c("curtisreid", "equilibrate")} in the 
option list for the function call.

If a numerical problem occurs, status=5, the best solution is probably to
scale the input \code{X} and output \code{Y} yourself. The best results are
obtained when the variables are close to 1. If the some variable are in the
millions, then let the unit of measure be a million units.

}



\seealso{
Paul W. Wilson (2008), \dQuote{FEAR 1.0: A Software Package for Frontier
Efficiency Analysis with R,} \emph{Socio-Economic Planning Sciences}
42, 247--254
}


\examples{
x <- matrix(c(100,200,300,500,100,200,600),ncol=1)
y <- matrix(c(75,100,300,400,25,50,400),ncol=1)
dea.plot.frontier(x,y,txt=TRUE)

e <- dea(x,y)
eff(e)
print(e)
summary(e)
lambda(e)

# Input savings potential for each firm
(1-eff(e)) * x
(1-e$eff) * x

# calculate slacks
el <- dea(x,y,SLACK=TRUE)
data.frame(e$eff,el$eff,el$slack,el$sx,el$sy)

# Fully efficient units, eff==1 and no slack
which(eff(e) == 1 & !el$slack)

# fdh+ with limits in the interval [.7, 1.2]
dea(x,y,RTS="fdh+", param=c(.7,1.2))
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
\keyword{ efficiency }% __ONLY ONE__ keyword per line
