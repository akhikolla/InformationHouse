\name{PlotProjectedPoints}
\alias{PlotProjectedPoints}
\title{
Plot Projected Points
}
\description{
plots XY data colored by Cls with ggplot2
}
\usage{
PlotProjectedPoints(Points,Cls,BMUorProjected=F,PlotLegend=FALSE,

xlab='X',ylab='Y',main="Projected Points",PointSize=2.5)
}
\arguments{
  \item{Points}{
  	[1:n,1:2] xy cartesian coordinates of a projection
}

  \item{Cls}{ numeric vector, given Classification in numbers: every element is the cluster number of a certain corresponding element of data.
}
  \item{BMUorProjected}{
		Default ==FALSE, If TRUE assuming BestMatches of ESOM instead of Projected Points
}
  \item{PlotLegend}{
  ...
}
  \item{xlab}{
  Optional: Label of the x axis
}
  \item{ylab}{
  Optional: Label of the y axis
}
  \item{main}{
  Optional: title
}
  \item{PointSize}{
  Optional: size of  points
}
  
}
\value{
ggobject of ggplot2
}
\author{
Michael Thrun
}