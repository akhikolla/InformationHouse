\name{LinTransform}
\alias{LinTransform}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Transform parameters for linear transformations on response catagories
}
\description{
This function is mainly used to translate parameters estimated with response options set to 0 and 1 to a model in which the response options are -1 and 1, but can be used for any linear transformation of response options.
}
\usage{
LinTransform(graph, thresholds, from = c(0L, 1L), to = c(-1L, 1L), a, b)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{graph}{
A matrix containing an Ising graph
}
  \item{thresholds}{
A vector containing thresholds
}
  \item{from}{
The original response encoding
}
  \item{to}{
The response encoding to transform to
}
  \item{a}{
The slope of the transformation. Overwrites \code{to}.
}
  \item{b}{
The intercept of the transformation. Overwrites \code{to}.
}
}

\author{
Sacha Epskamp <sacha.epskamp@gmail.com>
}

\examples{
N <- 4 # Number of nodes

# Ising parameters:
Graph <- matrix(sample(0:1,N^2,TRUE,prob = c(0.7, 0.3)),N,N) * rnorm(N^2)
Graph <- pmax(Graph,t(Graph)) / N
diag(Graph) <- 0
Thresh <- -(rnorm(N)^2)
Beta <- 1

p1 <- IsingLikelihood(Graph, Thresh, Beta, c(0,1))

a <- 2
b <- -1

# p2 <- IsingLikelihood(Graph/(a^2), Thresh/a - (b*rowSums(Graph))/a^2, Beta, c(-1,1))

p2 <- IsingLikelihood(LinTransform(Graph,Thresh)$graph, 
                      LinTransform(Graph,Thresh)$thresholds , 
                      Beta, c(-1,1))
LinTransform

round(cbind(p1[,1],p2[,1]),5)

plot(p1[,1],p2[,1])
abline(0,1)
}

