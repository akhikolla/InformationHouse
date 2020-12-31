\name{PCA}
\alias{PCA}
\title{Principal Component Analysis (PCA)}
\usage{
PCA(Data,OutputDimension=2,Scale=FALSE,Center=FALSE,PlotIt=FALSE,Cls)
}
\arguments{
\item{Data}{numerical matrix of data: n cases in rows, d variables in columns}
\item{OutputDimension}{Number of dimensions in the Outputspace, default=2}
  \item{Scale}{
	  a logical value indicating whether the variables should be scaled to have unit variance before the analysis takes place. The default is FALSE for consistency with S, but in general scaling is advisable. Alternatively, a vector of length equal the number of columns of x can be supplied. The value is passed to scale.}
	  \item{Center}{
		a logical value indicating whether the variables should be shifted to be zero centered. Alternately, a vector of length equal the number of columns of x can be supplied. The value is passed to scale
	}
	
	\item{PlotIt}{
	Default: FALSE, If TRUE: Plots the projection as a 2d visualization. 
	 OutputDimension>2: only the first two dimensions will be shown
	}
	\item{Cls}{[1:n,1] Optional,: only relevant if PlotIt=TRUE. Numeric vector, given Classification in numbers: every element is the cluster number of a certain corresponding element of data.
}
}
\value{

  \item{ProjectedPoints}{
 [1:n,OutputDimension], n by OutputDimension matrix containing coordinates of the Projectio
}

  \item{Rotation}{
the matrix of variable loadings (i.e., a matrix whose columns contain the eigenvectors)
}
  \item{sDev}{
the standard deviations of the principal components (i.e., the square roots of 
                                the eigenvalues of the covariance/correlation matrix, 
                                though the calculation is actually done with the singular values of the data matrix)
}
  \item{TransformedData}{
matrix  with PCA transformed Data
}
  \item{Center}{
   the centering used, or FALSE
}
  \item{Scale}{
 the scaling used, or FALSE
}
}
\description{
 Performs a principal components analysis on the given data matrix projection=SammonsMapping(Data)}
\author{Michael Thrun}
\note{
A wrapper for \code{\link[stats]{prcomp}}

 You can use the standard \code{ShepardScatterPlot} or the better approach through the \code{ShepardDensityPlot} of the CRAN package \code{DataVisualizations}.
}
 \details{
 An short overview of different types of projection methods can be found in [Thrun, 2018, p.42, Fig. 4.1] (\url{https://doi.org/10.1007/978-3-658-20540-9}).
 }
 
  \examples{
data('Hepta')
Data=Hepta$Data

Proj=PCA(Data)

\dontrun{
PlotProjectedPoints(Proj$ProjectedPoints,Hepta$Cls)
}
}

\keyword{PCA}

\concept{Principal component analysis}
\concept{Projection Method}
\concept{Dimensionality Reduction}
\keyword{DR}