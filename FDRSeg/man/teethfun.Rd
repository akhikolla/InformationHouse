\name{teethfun}
\alias{teethfun}

\title{Teeth function}
\description{Creat the teeth function with specified lengths and number of change-points. 
}

\usage{
teethfun(n, K, h = 3)
}

\arguments{
\item{n}{length of the vector (values of the teeth function)}
\item{K}{number of change-points}
\item{h}{height of the jump}
}

\value{
A numeric vector gives values of the teeth function. 
}

\references{
Fryzlewicz, P. (2014). Wild binary segmentation for multiple change-point detection.
\emph{Ann. Statist.}, 42(6): 2243--1572.
} 

\examples{
# create teeth function
u <- teethfun(100, 6)

# plot 
plot(u, type = "s")
}

\keyword{nonparametric}


