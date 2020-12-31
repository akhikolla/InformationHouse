% $Id: excess.Rd 206 2019-07-11 14:15:27Z lao $
\name{excess}
\alias{excess}

\title{
Excess input compared over frontier input
}

\description{ Excess input compared over frontier input and/or less
output than frontier/transformation/optimal output.  }

\usage{
excess(object, X = NULL, Y = NULL)
}


\arguments{ 

  \item{object}{ A Farrell object as returned from functions like
      \link{dea}, \link{dea.direct}, link{sdea}, and \link{mea}. }

  \item{X}{ Input matrix, only neccesary for ordinary input Farrell
      efficiency }

  \item{Y}{ Ouput matrix , only neccesary for ordinary output Farrell
      efficiency }

}

\details{

  For Farrell input efficiency E the exess input is \eqn{(1-E) X} and
  for Farrell ouput efficiency F the missing output is \eqn{(F-1) Y}.
  
  Notice that the excess calculated does not include any slackvalues. In 
  case slacks are present and calculated it might be more appropriate 
  to add slack; i.e. to use \code{excess(object, X, Y) + slack(X, Y, object)}.

  For directional efficiency e in the direction D the excess input is
  \eqn{e D}.

  If a firm is outside the technology set, as could be the case when
  calculating super-efficiencies, the Farrell input efficiency is
  larger then 1 and then the excess values are negative.

}


\value{
Return a matrix with exces input and/or less output.
}

\references{ Peter Bogetoft and Lars Otto; \emph{Benchmarking with
DEA, SFA, and R}; Springer 2011 }

\author{
Peter Bogeroft and Lars Otto \email{larsot23@gmail.com}
}




\examples{
x <- matrix(c(100,200,300,500,100,200,600),ncol=1)
y <- matrix(c(75,100,300,400,25,50,400),ncol=1)

e <- dea(x,y)
excess(e,x)
x - eff(e) * x

e <- dea(x,y, ORIENTATION="graph")
excess(e, x, y)
x - eff(e) * x
1/eff(e) * y -y

me <- mea(x,y)
excess(me)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
\keyword{ efficiency }% __ONLY ONE__ keyword per line
