\name{mphase1.methods}
\alias{postsignal}
\alias{print.mphase1}
\alias{plot.mphase1}
\alias{postsignal.mphase1}
\title{
  Methods for objects of class \code{mphase1}
}
\description{
  Methods \code{print} and \code{plot} allow to write
  to the console and plot (optionally changing the layout)
  the result of the Phase I analysis performed with function
  \code{mphase1}.
  
  Method \code{postsignal} implements the post-signal Phase I analysis
  based on the adaptive LASSO described in Capizzi and Masarotto (2016).
  It uses the p-value and the results on the forward search
  contained in its first argument. Hence, it is
  useful for re-running the analysis with different values
  of \code{alpha} and/or \code{gamma}.
}
\usage{
\method{print}{mphase1}(x,...)
\method{plot}{mphase1}(x,layout=c(1,p),...)
\method{postsignal}{mphase1}(x, plot = TRUE, alpha = 0.05, gamma = 0.5,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
    an object returned by function \code{\link{mphase1}}.
  }
  \item{layout}{
    an integer vector describing the multi-panel (and
    possible multi-page) layout. 
  }
  \item{plot}{
    logical; if \code{TRUE} the diagnostic plot is displayed.
  }
  \item{alpha}{
    real; the acceptable false alarm probability; if the
    observed p-value is greater than \code{alpha}, then
    the estimated mean function is a constant.
  }
  \item{gamma}{
    real; the extra penalization for the extended BIC criteria.
  }
  \item{...}{ignored.}
}
\value{
  An object of class \code{mphase1}. See \code{\link{mphase1}}
  for the description.
}
\references{
  G. Capizzi and G. Masarotto (2016), Phase I Distribution-Free
  Analysis of Multivariate Data, \emph{Technometrics},
  to appear, \url{http://dx.doi.org/10.1080/00401706.2016.1272494}.
}
\author{
\packageAuthor{dfphase1}.
}


\examples{
  data(gravel)
  u <- mphase1(gravel,plot=FALSE)
  print(u)
  plot(u,layout=c(2,1))
  postsignal(u,plot=FALSE,gamma=1)
}


