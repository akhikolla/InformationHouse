\name{partial.BoostMLR}
\alias{partial.BoostMLR}
\title{Partial plot analysis}
\description{
  Partial dependence plot of x and time against adjusted predicted y.
}
\usage{
\method{partial}{BoostMLR}(Object,
        xvar.name,
        n.x = 10,
        n.tm = 10,
        x.unq = NULL,
        tm.unq = NULL,
        Mopt,
        plot.it = TRUE,
        path_saveplot = NULL,
        Verbose = TRUE,
        ...)
        }
\arguments{
  \item{Object}{A boosting object of class \code{(BoostMLR, grow)}.}
  \item{xvar.name}{Name of the x-variable to be used for partial plot.}
  \item{n.x}{Maximum number of unique points used for \code{xvar.name}. Reduce this value if
    plotting is slow.}
  \item{n.tm}{Maximum number of unique points used for \code{tm}. Reduce this value if
    plotting is slow.}
  \item{x.unq}{Unique values used for the partial plot for variable \code{xvar.name}. 
  Default is NULL in which case
  unique values are obtained uniformaly based on the range of variable.}
  \item{tm.unq}{Unique time points used for the partial plots of x against y.
  Default is NULL in which case
  unique values are obtained uniformaly based on the range of \code{tm}.
  }
  \item{Mopt}{The optimal number of boosting iteration. If missing, the value from 
the \code{Object} will be used.} 
  \item{plot.it}{Should partial plot be displayed?}
  \item{path_saveplot}{
        Provide the location where plot should be saved. By default the plot will be
        saved at temporary folder.
}
\item{Verbose}{
       Display the path where the plot is saved?
}
  \item{...}{Further arguments passed to or from other methods.}
}
\details{
  Partial dependence plot (Friedman, 2001) of x values specified by
  \code{xvar.name} against the adjusted predicted y-values over a set
  of time points specified by \code{tm.unq}.
}

\value{
\item{x.unq}{Unique values used for the partial plot for variable \code{xvar.name}}
\item{tm.unq}{Unique time points used for the partial plots of x against y.}
\item{pList}{List with number of elements equal to number of multivariate response.
Each element of the list is a matrix with number of rows equal to length of \code{x.unq},
and number of columns equal to length of \code{tm.unq}. Values in the matrix represent
predicted partial values.}
\item{sList}{List with number of elements equal to number of multivariate response.
Each element is a matrix with the same dimension as described in \code{pList}.
Values are calculated using the local smoother (loess) for \code{tm.unq} and the i'th row
of the matrix from \code{pList}. Users are encouraged to use \code{pList} to
genenrate their own \code{sList} so that they will have more control over the
different arguments of local smoother.}
}
\author{
  Amol Pande and Hemant Ishwaran
}
\references{
  Friedman J.H. Greedy function approximation: a gradient
  boosting machine, \emph{Ann. of Statist.}, 5:1189-1232, 2001.
}
\examples{
\donttest{
##------------------------------------------------------------
## Generate partial plot for covariate x1
##-------------------------------------------------------------

dta <- simLong(n = 100, N = 5, rho =.80, model = 1, q_x = 0, 
                                  q_y = 0,type = "corCompSym")$dtaL

# Boosting call: Raw values of covariates, B-spline for time, 
# no shrinkage, no estimate of rho and phi

boost.grow <- BoostMLR(x = dta$features, tm = dta$time, id = dta$id, 
                            y = dta$y, M = 100, VarFlag = FALSE)

Partial_Plot_x1 <- partial.BoostMLR(Object = boost.grow, xvar.name = "x1",plot.it = FALSE)

}
}
\keyword{plot}
