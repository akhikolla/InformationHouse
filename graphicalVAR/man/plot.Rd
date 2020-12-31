\name{plot.graphicalVAR}
\alias{plot.graphicalVAR}
\title{
Plot method for graphicalVAR objects
}
\description{
Sends the estimated PCC and PDC networks to \code{\link[qgraph]{qgraph}}.
}
\usage{
\method{plot}{graphicalVAR}(x, include = c("PCC", "PDC"), repulsion = 1, 
                horizontal = TRUE, titles = TRUE, sameLayout = TRUE, 
                unweightedLayout = FALSE, ...)
}
\arguments{
  \item{x}{
A \code{graphicalVAR} object
}
  \item{include}{
A vector of at most two containing \code{"PCC"} and \code{"PDC"} indicating which networks should be plotted and in what order.
}
  \item{repulsion}{
The repulsion argument used in \code{\link[qgraph]{qgraph}}
}
  \item{horizontal}{
Logical, should the networks be plotted horizontal or vertical?
}
  \item{titles}{
Logical, should titles be added to the plots?
}
  \item{sameLayout}{
Logical, should both networks be plotted in the same layout?
}
 \item{unweightedLayout}{
Logical, should the layout be based on the unweighted network instead of the weighted network?
}
  \item{\dots}{
Arguments sent to \code{\link[qgraph]{qgraph}}
}
}
\author{
Sacha Epskamp <mail@sachaepskamp.com>
}