% $Id: dea.plot.Rd 121 2011-10-10 21:12:44Z Lars $
\name{dea.plot}
\alias{dea.plot}
\alias{dea.plot.frontier}
\alias{dea.plot.isoquant}
\alias{dea.plot.transform}


\title{Plot of DEA technologies}

\description{Draw a graph of a DEA technology. Designed for two goods
   illustrations, either isoquant (2 inputs), transformation curve (2
   outputs), or a production function (1 input and 1 output).  If the
   number of good is larger than 2 then aggregation occur, either
   simple or weighted.}

\usage{

dea.plot(x, y, RTS="vrs", ORIENTATION="in-out", txt=NULL, add=FALSE, 
         wx=NULL, wy=NULL, TRANSPOSE=FALSE, fex=1, GRID=FALSE,
         RANGE=FALSE, param=NULL, ..., xlim, ylim, xlab, ylab)

dea.plot.frontier(x, y, RTS="vrs",...)

dea.plot.isoquant(x1, x2, RTS="vrs",...)

dea.plot.transform(y1, y2, RTS="vrs",...)

}


\arguments{

  \item{x}{The good illustrated on the first axis. If there are more
      than 1 input then inputs are just summed or, if \code{wx} is
      present, a weighted sum of inputs is used.}

  \item{y}{The good illustrated on the second axis. If there are more
      than 1 output then outputs are just summed or, if \code{wy} is
      present, a weighted sum of outputs is used.}

  \item{x1, y1}{The good illustrated on the first axis}
  \item{x2, y2}{The good illustrated on the  second axis}

  \item{RTS}{Underlying DEA model / assumptions about returns to
      scale: "fdh" (0), "vrs" (1), "drs" (2), "crs" (3), "irs" (4),
      "irs2" (5) (irs without convexity), "add" (6), and "fdh+" (7).
      Numbers in parenthesis can also be used as values for
      \code{RTS}}

  \item{ORIENTATION}{Input-output graph of 1 input and 1 output is
      "in-out" (0), graph of 2 inputs is "in" (1), and graph of 2
      outputs is "out" (2). }

  \item{txt}{\code{txt} is an array to label the observations. If
      \code{txt=TRUE} the observations are labeled by the observation
      number or rownames if there are any. }

  \item{add}{For \code{add=T} the technology is drawn on top of an
      existing graph. With the default \code{add=F}, a new graph is
      made.}

  \item{wx}{Weight to aggregate the first axis if there are more than
      1 good behind the first axis.}
  
  \item{wy}{Weights to aggregate for the second axis if there are more
      than 1 good behind the second the second axis.}

  \item{TRANSPOSE}{Only relevant for more than 1 good for each axis,
      see \code{\link{dea}} for a description of this option.}

 \item{GRID}{If \code{GRIF=TRUE} a gray grid is put on the plot.}

  \item{\dots}{Usual options for the methods \code{plot}, \code{lines}, 
      and \code{abline} etc.}

  \item{fex}{Relative size of the text/labels on observations;
      corresponds to \code{cex}, but only changes the size of the
      text.}

  \item{RANGE}{A logical variable, if \code{RANGE=TRUE} the limits for
      the graph is the range of the variables; zero is always
      included. Default is \code{RANGE=FALSE} when the range is from
      zero to the max values.  Relevant if some values are negative.}

  \item{param}{Possible parameters. At the moment only used for
      RTS="fdh+"; see the section details and examples for its
      use. Future versions might also use \code{param} for other
      purposes.}

  \item{xlim}{Possible limits \code{c(x1,x2)} for the first axis}

  \item{ylim}{Possible limits  \code{c(y1,y2)} for the second axis}

  \item{xlab}{Possible label for the x-axis}

  \item{ylab}{Possible label for the y-axis}

}


\details{ The method \code{dea.plot} is the general plotting
   method. The the 3 others are specialized versions for frontiers (1
   input and 1 output), isoquant curves (2 inputs for given outputs),
   and transformation curves (2 outputs for given inputs) obtained by
   using the argument \code{ORIENTATION}.

   The crs factor in RTS="fdh+" that sets the lower and upper bound
   can be changed by the argument \code{param} that will set the lower
   and upper bound to 1-param and 1+param; the default value is
   param=.15. The value must be greater than or equal to 0 and
   strictly less than 1. A value of 0 corresponds to RTS="fdh". The
   FDH+ technology set is described in Bogetoft and Otto (2011) pages
   72--73.}

\note{If there are more than 1 good for the arguments \code{x} and
    \code{y} then the goods are just summed or, if \code{wx} or
    \code{wy} are present, weighted sum of goods are used. In this
    case the use of the command \code{identify} must be calles as
    \code{dea.plot(rowSums(x),rowSums(y))}. 

    \emph{Warning} If you use this facility to plot multi input and
    multi output then the plot may deceive you as fully multi
    efficient firms are not necessarely placed on the two dimensional
    frontier. 

    Note that \code{RTS="add"} and \code{RTS="fdh+"} only works for
    \code{ORIENTATION="in-out"} (0).  }


\value{No return, uses the original graphing system.  }


\references{ Peter Bogetoft and Lars Otto; \emph{Benchmarking with
   DEA, SFA, and R}; Springer 2011 

Paul Murrell; \emph{R Graphics}; Chapman & Hall 2006}


\author{ Peter Bogetoft and Lars Otto \email{larsot23@gmail.com}}


\seealso{The documentation for the function \code{plot} and Murrell
(2006) for further options and on customizing plots.  }


\examples{

x <- matrix(c(100,200,300,500,600,100),ncol=1)
y <- matrix(c(75,100,300,400,400,50),ncol=1)

dea.plot(x,y,RTS="vrs",ORIENTATION="in-out",txt=LETTERS[1:length(x)])
dea.plot(x,y,RTS="crs",ORIENTATION="in-out",add=TRUE,lty="dashed")

dea.plot.frontier(x,y,txt=1:dim(x)[1])

n <- 10
x <- matrix(1:n,,1)
y <- matrix(x^(1.6) + abs(rnorm(n)),,1)
dea.plot.frontier(x,y,RTS="irs",txt=1:n)
dea.plot.frontier(x,y,RTS="irs2",add=TRUE,lty="dotted")

# Two different forms of irs: irs and irs2, and two different ways to
# make a frontier
id <- sample(1:n,30,replace=TRUE)
dea.plot(x[id],y[id],RTS="irs",ORIENTATION="in-out")
dea.plot.frontier(x[id],y[id],RTS="irs2")

# Difference between the FDH technology and the additive 
# FRH technology
x <- matrix(c(100,220,300,520,600,100),ncol=1)
y <- matrix(c(75,100,300,400,400,50),ncol=1)
dea.plot(x,y,RTS="fdh",ORIENTATION="in-out",txt=LETTERS[1:length(x)])
dea.plot(x,y,RTS="add",ORIENTATION="in-out",add=TRUE,lty="dashed",lwd=2)
dea.plot(x,y,RTS="fdh+",ORIENTATION="in-out",add=TRUE,
                            lty="dotted",lwd=3,col="red")

# Use of parameter in FDH+
dea.plot(x,y,RTS="fdh",ORIENTATION="in-out",txt=LETTERS[1:length(x)])
dea.plot(x,y,RTS="fdh+",ORIENTATION="in-out",add=TRUE,lty="dashed")
dea.plot(x,y,RTS="fdh+",ORIENTATION="in-out",add=TRUE,lty="dotted",param=.5)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
\keyword{plot}
