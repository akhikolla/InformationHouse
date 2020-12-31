\name{preston_plot}
\alias{preston_plot}

\title{
  Barplot in Preston style of an abundance dataset
}
\description{
  This function first sorts abundances into octaves, and then plots the resulting distribution.
}
\usage{
  preston_plot(abund, expected, ...)
}

\arguments{
  \item{abund}{vector containing the number of individuals per species}
  
  \item{expected}{vector containing the expected number of species per octave}
  
  \item{\dots}{further graphical arguments that can be passed to \code{barplot()}}
}

\author{
  Thijs Janzen
}

\examples{
  theta = 10
  m = 0.1
  J = 1000
  I = m * (J - 1) / (1 - m)
  
  abund <- generate.ESF(theta, I, J)
  par(mfrow = c(1,2))
  preston_plot(abund)
  abund.expect <- expected.SAD(theta, m, J)
  preston_plot(abund, abund.expect)
}