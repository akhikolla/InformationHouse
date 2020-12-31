\name{mpboost}

\alias{mpboost}

\title{
Compute Allocation Sequences by the Maximal Procedure
}

\description{
This function produces a treatment assignment sequence generated
according to the maximal procedure of Berger et al. (2003). It is an
implementation of the algorithm proposed by Salama et al. (2008).
}

\usage{
mpboost(N1, N2, MTI = 2L)
}

\arguments{
\item{N1}{An integer specifying the size of the sample assigned to treatment 1.
}
\item{N2}{An integer specifying the size of the sample assigned to treatment 2.
}
\item{MTI}{An integer specifying the maximum tolerated imbalance
  (MTI). The default is 2.
}
}

\value{A vector of N1 1's and N2 2's representing the allocation sequence.}

\references{Berger, V. W., Ivanova, A., Knoll, M. D. (2003). Minimizing
  predictability while retaining balance through the use of less
  restrictive randomization procedures. Statistics in Medicine, 22: 3017-3028.
  \url{https://doi.org/10.1002/sim.1538}.
  
  Salama, I., Ivanova, A., Qaqish, B. (2008). Efficient generation of
  constrained block allocation sequences. Statistics in Medicine, 27,
  1421-1428. \url{https://doi.org/10.1002/sim.3014}.
}

\note{See the package's vignette for more details and further examples of the use of the function.
} 

\examples{
#\donttest{
mpboost(N1 = 25, N2 = 25, MTI = 3)
#}
}

\keyword{clinical trial}

\keyword{constrained block randomization}
