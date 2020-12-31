% $Id: dea.merge.Rd 229 2020-07-04 13:39:18Z lao $
\name{dea.merge}
\alias{dea.merge}

\title{ Estimate potential merger gains and their decompositions }

\description{ Calculate and decompose potential gains from mergers of
similar firms (horizontal integration).  }

\usage{

dea.merge(X, Y, M, RTS = "vrs", ORIENTATION = "in",
          XREF = NULL, YREF = NULL, FRONT.IDX = NULL, TRANSPOSE=FALSE,
          CONTROL=NULL)
}

\arguments{\item{}{ Most of the arguments correspond to the arguments in
\code{\link{dea}}, with K firms, m inputs, and n outputs.}

  \item{X}{ K times m matrix as in \code{dea} }

  \item{Y}{ K times n matrix as in \code{dea} }

  \item{M}{Kg times K matrix where each row defines a merger by the
  firms (collums) included; matrix as returned from method
  \code{\link{make.merge}}}

  \item{RTS}{as in \code{\link{dea}} }

  \item{ORIENTATION}{ as in \code{\link{dea}} }

  \item{XREF}{ as in \code{dea} }

  \item{YREF}{ as in \code{dea} }

  \item{FRONT.IDX}{ as in \code{dea} }

  \item{TRANSPOSE}{ as in \code{dea} }

  \item{CONTROL}{Possible controls to \pkg{lpSolveAPI}, see the
      documentation for that package. For examples of use see the
      function \code{\link{dea}}.}

}


\details{ 
The K firms are merged into Kg new, merged firms.

The decomposition is summarized on page 275 and in table 9.1 page 276
in Bogetoft and Otto (2011) and is based on Bogetoft and Wang (2005) }


\value{

\item{Eff}{Overall efficiencies of mergers, Kg vector}

\item{Estar}{Adjusted overall efficiencies of mergers after the
removal of individual learning, Kg vector}

\item{learning}{Learning effects, Kg vector}

\item{harmony}{Harmony (scope) effects, Kg vector}

\item{size}{Size (scale) effects, Kg vector}

}


\references{ Bogetoft and Otto; \emph{Benchmarking with DEA, SFA, and
R}; chapter 9; Springer 2011.

Bogetoft and Wang; \dQuote{Estimating the Potential Gains from
Mergers}; \emph{Journal of Productivity Ana-lysis}, 23, pp. 145-171,
2005. }


\author{ Peter Bogetoft and Lars Otto \email{larsot23@gmail.com}}


% \note{
%%  ~~further notes~~
% }


\seealso{
    \code{\link{dea}} and \code{\link{make.merge}}
}



\examples{
x <- matrix(c(100,200,300,500),ncol=1,dimnames=list(LETTERS[1:4],"x"))
y <- matrix(c(75,100,300,400),ncol=1,dimnames=list(LETTERS[1:4], "y"))

dea.plot.frontier(x,y,RTS="vrs",txt=LETTERS[1:length(x)],
    xlim=c(0,1000),ylim=c(0,1000) )
dea.plot.frontier(x,y,RTS="drs", add=TRUE, lty="dashed", lwd=2)
dea.plot.frontier(x,y,RTS="crs", add=TRUE, lty="dotted")

dea(x,y,RTS="crs")
M <- make.merge(list(c(1,2), c(3,4)), X=x)
xmer <- M \%*\% x
ymer <- M \%*\% y
points(xmer,ymer,pch=8)
text(xmer,ymer,labels=c("A+B","C+D"),pos=4)
dea.merge(x,y,M, RTS="vrs")
dea.merge(x,y,M, RTS="crs")
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
\keyword{ efficiency }
\keyword{ merger }
