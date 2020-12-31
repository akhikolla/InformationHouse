\name{tri.kern}
\alias{tri.kern}

\title{Trapezoidal shape taper}

\description{
Returns a half of the trapezoidal shape taper symmetric around zero, for average block length selection
}

\usage{
tri.kern(h)
}

\arguments{
  \item{h}{ bandwidth }
}

\value{
A vector of length \code{h+1} containing the trapezoidal shape taper symmetric around zero,
for indices greater than or equal to zero.
}

\author{
Haeran Cho
}

\keyword{internal}
