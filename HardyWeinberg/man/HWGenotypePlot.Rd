\name{HWGenotypePlot}
\alias{HWGenotypePlot}
\title{
Scatter plot of the genotype frequencies 
}
\description{
\code{HWGenotypePlot} makes a scatterplots of the AB or BB frequency
versus the AA frequency and represents a blue curve indicating the
Hardy-Weinberg equilibrium condition.}
\usage{
HWGenotypePlot(X, plottype = 1, xlab = expression(f[AA]), ylab =
ifelse(plottype == 1, expression(f[AB]), expression(f[BB])), asp = 1,
pch = 19, xlim = c(0, 1), ylim = c(0, 1), cex = 1, cex.axis = 2, cex.lab = 2, ...)
}
\arguments{
    \item{X}{A matrix of genotype counts or frequencies with three columns (AA, AB, BB)}
    \item{plottype}{\code{plottype=1} produces a plot of AB versus AA,
      \code{plottype=2} produced a plot of BB versus AA.}
  \item{xlab}{A label for the x axis}
  \item{ylab}{A label for the y axis}
  \item{asp}{Aspec ratio (1 by default)}
  \item{pch}{Plotting charachter (19 by default)}
  \item{xlim}{Limits for the x axis (0-1 by default)}
  \item{ylim}{Limits for the y axis (0-1 by default)}
  \item{cex}{Character expansion factor (1 by default)}
  \item{cex.axis}{Character expansion factor for the axes (2 by default)}
  \item{cex.lab}{Character expansion factor for labels of axis (2 by default)}
  \item{\dots}{Additional arguments for the \code{plot} function}
}
\value{NULL}
\author{ Jan Graffelman \email{jan.graffelman@upc.edu} }
\seealso{ \code{\link{HWTernaryPlot}} }
\examples{
n <- 100 # sample size
m <- 100 # number of markers
Xc <- HWClo(HWData(n,m))
HWGenotypePlot(Xc,plottype=1,main="Heterozygote-homozygote scatterplot")
}
\keyword{aplot}

