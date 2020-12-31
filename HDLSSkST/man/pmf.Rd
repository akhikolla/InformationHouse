\name{pmf}
\alias{pmf}
\title{Generalized Hypergeometric Probability
}
\description{A function that provides the probability of observing an \eqn{r\times c} contingency table using generalized hypergeometric probability.
}
\usage{
pmf(M)
}
\arguments{
  \item{M}{\eqn{r\times c} contingency table
}
}

\value{
a single value between \eqn{0} and \eqn{1}
}
\references{
Biplab Paul, Shyamal K De and Anil K Ghosh (2020).  Distribution-free Exact k-sample Tests for High Dimensional Low Sample Size data Based on Clustering.

Cyrus R Mehta and Nitin R Patel (1983). A network algorithm for performing Fisher's exact test in rxc contingency tables, \emph{Journal of the American Statistical Association}, 78(382):427-434, doi:10.2307/2288652.
}
\author{ Biplab Paul, Shyamal K. De and Anil K. Ghosh

Maintainer: Biplab Paul<biplab.paul@niser.ac.in>
}
\examples{
   # Generalized hypergeometric probability of rxc Contingency Table:
   mat <- matrix(1:20,5,4, byrow = TRUE)
   pmf(mat)

   ## outputs:
   #[1] 4.556478e-09
}
\keyword{ pmf function}