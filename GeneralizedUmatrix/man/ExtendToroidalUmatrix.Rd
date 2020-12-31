\name{ExtendToroidalUmatrix}
\alias{ExtendToroidalUmatrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Extend Toroidal Umatrix
}
\description{
Extends Umatrix by toroidal continuation of the given Umatrix defined by \code{ExtendBorders} in all four directions.
}
\usage{
ExtendToroidalUmatrix(Umatrix, Bestmatches, ExtendBorders)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Umatrix}{[1:Lines,1:Columns] Matrix of Umatrix Heights}
  \item{Bestmatches}{[1:n, 1:2] Matrix with positions of Bestmatches for n datapoints, first columns is the position in \code{Lines} and second column in \code{Columns}}
  \item{ExtendBorders}{
number of lines and columns the umatrix should be extended with
}
}
\details{
Function assumes that U-matrix is not planaer (has no borders), i.e. is toroidal, and not tiled. Bestmatches are moved to new positions accordingly. Example is shown in conference talk of [Thrun et al., 2020].
}
\value{
  \item{Umatrix}{[1:Lines+2*ExtendBorders,1:Columns+2*ExtendBorders] Matrix of  U-Heights}
  \item{Bestmatches}{Array with positions of Bestmatches }
}
\references{
[Thrun et al., 2020]  Thrun, M. C., Pape, F., & Ultsch, A.: Interactive Machine Learning Tool for Clustering in Visual Analytics, 7th IEEE International Conference on Data Science and Advanced Analytics (DSAA 2020), Vol. accepted, pp. 1-9, IEEE, Sydney, Australia, 2020.
}
\author{
Michael Thrun
}
\note{
Currently can be only used if untiled U-Matrix (the default) is presented, but 4-tiled U-matrix does not work.
}

\examples{
#ToDO
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{ExtendToroidalUmatrix}
\keyword{Umatrix}
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
\concept{Extend Toroidal Umatrix}
\concept{Extend Toroidal U-matrix}
% Use only one concept per line.
