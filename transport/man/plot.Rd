\name{plot}
\alias{plot.pgrid}
\alias{plot.pp}
\alias{plot.wpp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Methods for Plotting Objects of Class pgrid, pp and wpp
}
\description{
  Methods for plotting objects of class pgrid, pp and wpp,
    possibly together with a transference plan.
}
\usage{
  \method{plot}{pgrid}(x, y = NULL, tplan = NULL, mass = c("colour", "thickness"),
    length = 0.1, acol, bcol=4, lwd, rot = FALSE, overlay = FALSE, static.mass =TRUE, ...)
  \method{plot}{pp}(x, y = NULL, tplan = NULL, cols = c(4, 2), cex = 0.8,
    acol = grey(0.3), lwd = 1, overlay = TRUE, ...)
  \method{plot}{wpp}(x, y = NULL, tplan = NULL, pmass=TRUE, tmass=TRUE, cols = c(4, 2),
    cex = 0.8, aglevel = 0.4, acol = grey(0.3), lwd = 1, overlay = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x,y}{
    one or two objects of class \code{pgrid} or class \code{pp} to be plotted.
  }
  \item{tplan}{
    a transference plan between the two objects \code{x} and \code{y},
    typically an optimal transference plan obtained by a call to \code{\link{transport}}.
  }
  \item{mass, pmass, tmass}{
    for \code{pgrid} objects with a \code{tplan}:
    if \code{mass == "colour"}, the mass transferred is depicted by heatmap colours;
    if \code{mass == "thickness"}, it is depicted by the line widths of the arrows.\cr
    For \code{wpp} objects:
    \code{pmass, tmass} are logicals controlling whether the \emph{amount} of mass
    associated with the points and the mass transferred should be depicted in
    the plot.
  }
  \item{length}{
    the length of the arrow heads in inches.
  }
  \item{aglevel}{
  	for \code{wpp} objects with \code{tmass = TRUE}: the grey level
  	chosen for depicting the transport of an average amount of mass.
  }
  \item{acol}{
    the colour of the arrows/lines of the transference plan. Ignored for \code{pgrid}
    objects if \code{mass = "colour"} and for \code{wpp} objects if \code{tmass}
    is \code{TRUE}.  
  }   
  \item{bcol}{
    the colour of the cell boundaries for a semidiscrete transport plan. Ignored
    in all other instances.  
  }    
  \item{cols}{
  	for \code{pp} objects:
    A vector of size 2 specifying the colours of the two \code{pp} objects.
  }
  \item{cex,lwd,\dots}{
  	further graphic parameters used by plot. Note that for pgrid objects
  	\code{acol} is ignored for \code{mass == "colour"}, and \code{lwd} is ignored
  	for \code{mass == "thickness"}. Setting any of these parameters is optional.
  }
  \item{rot}{
  	logical. Whether the mass matrices of pgrid objects should be rotated before calling
  	\code{\link{image}} so that the orientation of the ploted pixelgrid and the orientation
  	of the mass matrix are the same. Otherwise plotting follows the usual convention of
  	\code{\link{image}}.
  }
  \item{overlay}{
  	in the case of two objects \code{x} and \code{y} whether they should be plotted on top
  	of one another (for pgrid objects the difference \code{x-y} is plotted) or not. In the
  	presence of a transference plan \code{overlay} is forced to be true.
  }
  \item{static.mass}{
  	for a transference plan that explicitly lists the \dQuote{static mass transports} (i.e.
  	mass that stays at the same site), should these transports also be plotted as disks
  	with colours/sizes corresponding to the amount of mass that stays?\cr
  	Note that it is wrong to assume that an optimal transference plan obtained by one of
  	the algorithms will automatically list static mass transports. It is not the
  	case for \eqn{p=1}, where static mass transport at site \eqn{i} is trivially equal
  	to the minimum of source mass and target mass, and it is currently not the case
  	for results obtained by \code{method="aha"}.
  }
}
\value{
  Used for its side effect.  
}
\author{
  Dominic Schuhmacher \email{dschuhm1@uni-goettingen.de}
}

