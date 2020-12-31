\name{HWIlr}
\alias{HWIlr}
\title{
Compute isometric log ratio coordinates.
}
\description{
\code{HWIlr} computes isometric log ratio coordinates for genotypic
compositions (AA, AB, BB)
}
\usage{
HWIlr(X, zeroadj = 0.5)
}
\arguments{
  \item{X}{A matrix of genotype counts, markers in rows, counts for AA,
    AB and BB in three columns
}
  \item{zeroadj}{Adjustment for zeros (0.5 by defaults)
}
}
\value{
A matrix of log ratio coordinates.
}
\references{Egozcue, J.J., Pawlowsky-Glahn, V., Mateu-Figueras, G. and
  Barcelo-Vidal, C. (2003) Isometric Logratio Transformations for Compositional
  Data Analysis. Mathematical Geology 35(3), pp. 279-300.
  
  Graffelman, J. and Egozcue, J. J. (2011) Hardy-Weinberg
  equilibrium: a non-parametric compositional approach. In: Vera
  Pawlowsky-Glahn and Antonella Buccianti (eds.) Compositional Data
  Analysis: Theory and Applications, John Wiley \& Sons, Ltd, pp. 207-215
  }
\author{ Jan Graffelman (jan.graffelman@upc.edu) }
\seealso{
\code{\link{HWAlr},\link{HWClr}}
}
\examples{
X <- HWData(100,100)
Y <- HWIlr(X)
}
\keyword{misc}

