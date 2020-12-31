\name{randomGVARmodel}
\alias{randomGVARmodel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulate a graphical VAR model
}
\description{
Simulates an contemporaneous and temporal network using the method described by Yin and Li (2001)
}
\usage{
randomGVARmodel(Nvar, probKappaEdge = 0.1, probKappaPositive = 0.5, probBetaEdge = 0.1, 
      probBetaPositive = 0.5, maxtry = 10, kappaConstant = 1.1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Nvar}{
Number of variables
}
  \item{probKappaEdge}{
Probability of an edge in contemporaneous network
}
  \item{probKappaPositive}{
Proportion of positive edges in contemporaneous network
}
  \item{probBetaEdge}{
Probability of an edge in temporal network
}
  \item{probBetaPositive}{
  Propotion of positive edges in temporal network
}
  \item{maxtry}{
Maximum number of attempts to create a stationairy VAR model
}
  \item{kappaConstant}{
The constant used in making kappa positive definite. See Yin and Li (2001)
}
}
\details{
The resulting simulated networks can be plotted using the plot method.
}
\value{
A list containing:
\item{kappa}{True kappa structure (residual inverse variance-covariance matrix)}
\item{beta}{True beta structure}
\item{PCC}{True partial contemporaneous correlations}
\item{PDC}{True partial temporal correlations}
}
\references{
Yin, J., & Li, H. (2011). A sparse conditional gaussian graphical model for analysis of genetical genomics data. The annals of applied statistics, 5(4), 2630-2650.
}
\author{
Sacha Epskamp
}
