% $Id: slack.Rd 229 2020-07-04 13:39:18Z lao $
\name{slack}
\alias{slack}
\alias{print.slack}
\alias{summary.slack}


\title{
Calculate slack in an efficiency analysis
}


\description{
Slacks are calculated after taking the efficiency into consideration.
}


\usage{
slack(X, Y, e, XREF = NULL, YREF = NULL, FRONT.IDX = NULL, LP = FALSE, CONTROL=NULL)
}



\arguments{

  \item{X}{Inputs of firms to be evaluated, a K x m matrix of
      observations of K firms with m inputs (firm x input). }

  \item{Y}{Outputs of firms to be evaluated, a K x n matrix of
      observations of K firms with n outputs (firm x input). }

  \item{e}{A Farrell object as returned from \code{\link{dea}} et al.}

  \item{XREF}{Inputs of the firms determining the technology, defaults
      to \code{X}}

  \item{YREF}{Outputs of the firms determining the technology,
      defaults to \code{Y}}

  \item{FRONT.IDX}{Index for firms determining the technology}

  \item{LP}{Set \code{TRUE} for debugging.}
  
    \item{CONTROL}{Possible controls to \pkg{lpSolveAPI}, see the
        documentation for that package. For examples of use see the
        function \code{\link{dea}}.}
}


\details{ Slacks are calculated in a LP problem where the sum of all
   slacks are maximied after correction for efficiency. 
    The for calculating slacks for orientation graph is low because of the low precision in the calculated graph efficiency. }


\value{The result is returned as the Farrell object used as the
   argument in the call of the function with the following added
   components:

  \item{slack}{A logical vector where the component for a firm is
      \code{TRUE} if the sums of slacks for the corresponding firm is
      positive. Only calculated in dea when option \code{SLACK=TRUE}}

  \item{sum}{A vector with sums of the slacks for each firm. Only
      calculated in dea when option \code{SLACK=TRUE}}


  \item{sx}{A matrix for input slacks for each firm, only calculated
      if the option \code{SLACK} is \code{TRUE} or returned from the
      method \code{slack}}

  \item{sy}{A matrix for output slack, see \code{sx}}

}


\references{ Peter Bogetoft and Lars Otto; \emph{Benchmarking with
DEA, SFA, and R}; Springer 2011. Sect. 5.6 page 127.  

WW Cooper, LM Seiford, and K Tone; \emph{Data Envelopment Analysis: A
Comprehensive Text with Models, Applications, References and
DEA-Solver Software}, 2nd edn. Springer 2007 .  }


\author{ Peter Bogetoft and Lars Otto \email{larsot23@gmail.com}}


\examples{
x <- matrix(c(100,200,300,500,100,200,600),ncol=1)
y <- matrix(c(75,100,300,400,25,50,400),ncol=1)
dea.plot.frontier(x,y,txt=1:dim(x)[1])

e <- dea(x,y)
eff(e)

# calculate slacks
sl <- slack(x,y,e)
data.frame(e$eff,sl$slack,sl$sx,sl$sy)
}

\keyword{models}
\keyword{efficiency}
\keyword{slack}% __ONLY ONE__ keyword per line
