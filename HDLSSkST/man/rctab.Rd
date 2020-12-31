\name{rctab}
\alias{rctab}
\title{Generates an \eqn{r\times c} Contingency Table
}
\description{A function that generates an \eqn{r\times c} contingency table with the same marginal totals as given \eqn{r\times c} contingency table.
}
\usage{
rctab(M)
}
\arguments{
  \item{M}{\eqn{r\times c} contingency table
}
}

\value{
generated \eqn{r\times c} contingency table
}
\references{
Cyrus R Mehta and Nitin R Patel (1983). A network algorithm for performing Fisher's exact test in rxc contingency tables, \emph{Journal of the American Statistical Association}, 78(382):427-434, doi:10.2307/2288652.
}
\author{ Biplab Paul, Shyamal K. De and Anil K. Ghosh

Maintainer: Biplab Paul<biplab.paul@niser.ac.in>
}
\examples{
   # Generation of rxc Contingency Table:
   set.seed(151)
   mat <- matrix(1:20,5,4, byrow = TRUE)
   rctab(mat)

   ## outputs:
   #      [,1] [,2] [,3] [,4]
   #[1,]    3    4    0    3
   #[2,]    4    5   10    7
   #[3,]    8    7   12   15
   #[4,]   18   16   13   11
   #[5,]   12   18   20   24
}
\keyword{ rctab function}