\name{upscaleUmatrix}
\alias{upscaleUmatrix}
\title{
Upscale a Umatrix grid
}
\description{
Use linear interpolation to increase the size of a umatrix. This can be used to produce nicer ggplot plots in \code{\link{plotTopographicMap}} and is going to be used for further normalization of the umatrix.
}
\usage{
upscaleUmatrix(Umatrix, Factor = 2,BestMatches, Imx)
}
\arguments{
  \item{Umatrix}{
    The umatrix which should be upscaled
  }
    \item{BestMatches}{
    The BestMatches which should be upscaled
  }
  \item{Factor}{
    Optional: The factor by which the axes will be scaled. Be aware that the size of the matrix will grow by Factor squared. Default: 2
  }
  \item{Imx}{
    Optional: Island cutout of the umatrix. Should also be scaled to the new size of the umatrix.
  }
}
\value{ A List consisting of:
\item{Umatrix}{A matrix representing the upscaled umatrix.}
\item{BestMatches}{If BestMatches was given as parameter: The rescaled BestMatches for an island cutout. Otherwise: \code{NULL}}
\item{Imx}{If Imx was given as parameter: The rescaled matrix for an island cutout. Otherwise: \code{NULL}}
}
\author{
Felix Pape
}
\concept{topographic map}
\concept{Generalized U-matrix}
\keyword{U-matrix}