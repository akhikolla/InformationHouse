% $Id: outlier.Rd 226 2020-06-27 21:48:14Z lao $
\name{outlier.ap}
\alias{outlier.ap}
\alias{outlier.ap.plot}
\alias{outlierC.ap}


\title{
Detection of outliers in benchmark models
}

\description{ The functions implements the Wilson (1993) outlier
detection method. One written entirely in R and another written in C++.}

\usage{
outlier.ap (X, Y, NDEL = 3, NLEN = 25, TRANSPOSE = FALSE)
outlierC.ap(X, Y, NDEL = 3, NLEN = 25, TRANSPOSE = FALSE)

outlier.ap.plot(ratio, NLEN = 25, xlab = "Number of firms deleted", 
                ylab = "Log ratio", ..., ylim)
}


\arguments{

  \item{X}{ Input as a firms times goods matrix, see \code{TRANSPOSE}.
}

  \item{Y}{ Output as a firms times goods matrix, see
\code{TRANSPOSE}.  }

  \item{NDEL}{ The maximum number of firms to be considered as a group
of outliers, i.e. the maximum number of firms to be deleted.  }

  \item{NLEN}{ The number of ratios to save for each level
or removal, the number of rows in \code{ratio} used.  }

\item{TRANSPOSE}{Input and output matrices are treated as firms
      times goods matrices for the default value \code{TRANSPOSE=FALSE}
      corresponding to the standard in R for statistical models.  When
      \code{TRUE} data matrices are transposed to good times firms
      matrices as is normally used in LP formulation of the problem.}


  \item{ratio}{ The \code{ratio} component from the list as output
from \code{outlier.ap}.} 

  \item{xlab}{ Label for the x-axis.  }

  \item{ylab}{ Label for the y-axis }

  \item{ylim}{ The y limits \code{(y1, y2)} of the plot, an
array/vector of length 2.  }

 \item{\dots}{Usual options for the methods \code{plot} and
 \code{lines}.}
 
}


\details{ 

An implementation of the method in Wilson (1993) using only R
functions and especially the function \code{det} to calculate
\eqn{R^{(i)}_{\min}}. The alternative method \code{outlierC.ap} 
is written completely in C++ and is much faster, but still not
as fast at the method in \pkg{FEAR}.

An elementary presentation of the method is found in Bogetoft and Otto
(2011), Sect. 5.13 on outliers.

For a data set with 10 firms and considering at the most 3 outliers there are
175 combinations of firms to delete. For 100 firms there are 166,750
combinations and for at most 5 outliers there are 79,375,495 combinatins, for
at most 8 outliers there are  203,366,882,995 combinations. For 200 firms whith
respectively 3,5 and 8 outliers there are  1,333,500, and 2,601,668,490, and 
a number we do not know what to call
57,467,902,686,615 combinations. Thus the number of combinations are increasing
exponentialy in both number of firms and number of firms to be deleted and so
is the computational time. Thus you should limit the numbers \code{NDEL} to a
very small number like at the most 3 or perhabs 5 depending of the number
of firms. Or you should use the extremely fast method \code{ap} from the
package \pkg{FEAR} mendtioned in the references.
}


\value{

\item{ratio}{A \code{min(NLEN,K) x NDEL} matrix with the log-ratios to
 be plotted.}

\item{imat}{A \code{NDEL x NDEL} matrix with indicies for deleted firms.}

\item{r0}{A \code{NDEL} array with the minimum value \eqn{R^{i}} of
the for each number of deleted firms. }

}


\references{ Bogetoft and Otto; \emph{Benchmarking with DEA, SFA, and
R}; Springer 2011

Wilson (1993), \dQuote{Detecing outliers in deterministic
nonparametric frontier models with multiple outputs,} \emph{Journal of
Business and Economic Statistics} 11, 319-323.

Wilson (2008), \dQuote{FEAR 1.0: A Software Package for Frontier
Efficiency Analysis with R,} \emph{Socio-Economic Planning Sciences}
42, 247--254

}


\author{ Peter Bogetoft and Lars Otto \email{larsot23@gmail.com}}


\note{ 

The function \code{outlier.ap} is extremely slow and for \code{NDEL}
larger than 3 or 4 it might be advisable to use the function \code{ap}
from the package \pkg{FEAR}.

The name of the returned components are the same as for \code{ap} in
the package \pkg{FEAR}.

}


\seealso{
The function \code{ap} in the package \pkg{FEAR}.
}


\examples{
n <- 25
x <- matrix(rnorm(n))
y <- .5 + 2.5*x + 2*rnorm(25)
tap <- outlier.ap(x,y, NDEL=2)
print(cbind(tap$imat,tap$rmin), na.print="", digit=2)
outlier.ap.plot(tap$ratio)
}

