\name{ProjectionPursuit}
\alias{ProjectionPursuit}
\title{ Projection Pursuit}
\arguments{
\item{Data}{array of data: n cases in rows, d variables in columns, matrix is not symmetric or distance matrix, in this case matrix has to be symmetric}
  \item{OutputDimension}{Number of dimensions in the Outputspace, default=2}
	
	\item{Indexfunction}{
	Criterium for Minimization:
	
	default: 'logcosh' G(u)=1/a*log cosh(a*u) (ICA)
                            'exp': G(u)=-exp(u^2/2)
                            'kernel'  1/(1* pi )*exp(r/2)
	}
	\item{Alpha}{
	constant with 1<=alpha<=2 used in approximation to neg-entropy when fun == "logcosh"
	}
	\item{Iterations}{
	 maximum number of iterations to perform.
	}
	\item{PlotIt}{
	Default: FALSE, If TRUE: Plots the projection as a 2d visualization. 
	 OutputDimension>2: only the first two dimensions will be shown
	}
	\item{Cls}{[1:n,1] Optional,: only relevant if PlotIt=TRUE. Numeric vector, given Classification in numbers: every element is the cluster number of a certain corresponding element of data.
}
}
\usage{
ProjectionPursuit(Data,OutputDimension=2,Indexfunction="logcosh",

Alpha=1,Iterations=200,PlotIt=FALSE,Cls)
}
\value{  \item{ProjectedPoints}{
 [1:n,OutputDimension], n by OutputDimension matrix containing coordinates of the Projectio
}
}
\description{
 In the absence of a generative model for the data the algorithm can be used to find the projection pursuit directions. Projection pursuit is a technique for finding 'interesting' directions in multidimensional  datasets
}
\author{
Michael Thrun}

\note{
 You can use the standard \code{ShepardScatterPlot} or the better approach through the \code{ShepardDensityPlot} of the CRAN package \code{DataVisualizations}.
}

 \details{
 An short overview of different types of projection methods can be found in [Thrun, 2018, p.42, Fig. 4.1] (\url{https://doi.org/10.1007/978-3-658-20540-9}).
 }
 
 \keyword{ProjectionPursuit}
  \concept{Projection Pursuit}
 \concept{Projection Method}
\concept{Dimensionality Reduction}
\keyword{DR}