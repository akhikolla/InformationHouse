\name{KruskalStress}
\alias{KruskalStress}
\title{
Kruskal stress calculation
}
\description{
Calculates the stress as defined by Kruskal for 2 distance matrices
}
\usage{
KruskalStress(InputDistances, OutputDistances)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{InputDistances}{
Distance matrix of the original Data
}
  \item{OutputDistances}{
Distance matrix of the projected Data
}
}
\value{
A single numerical representing the Kruskal stress of the distance matrices.
}
\author{
Felix Pape
}
 \details{
 An short overview of different types of quality measures can be found in [Thrun, 2018, p.68, Fig. 6.3] (\url{https://doi.org/10.1007/978-3-658-20540-9}).
}