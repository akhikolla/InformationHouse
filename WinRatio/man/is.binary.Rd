\name{is.binary}
\alias{is.binary}
\title{Binary Vectors}
\usage{
is.binary(x)
}
\arguments{
\item{x}{a vector to be tested.}
}
\value{
is.binary returns TRUE or FALSE depending on whether a vector is of binary type or not.
}
\description{
is.binary is used to test if a vector is of binary type, i.e. if a vector contains only 0 and/or 1 values.
}
\examples{
x1 <- runif(100, 0, 1)
x2 <- ifelse(x1 < 0.5, 1, 0)
is.binary(x1)
is.binary(x2)
}


