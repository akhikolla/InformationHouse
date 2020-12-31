\encoding{utf8}

\name{ElegantPairingInt}

\alias{ElegantPairingInt}
\alias{elegant pairing}
\title{Elegant Pairing}

\usage{
ElegantPairingInt(int1, int2)
}
\description{
Unordered Pairing Function creating a new unique integer from two input integers.
}

\arguments{
  \item{int1}{first integer to be paired.}
  \item{int2}{second integer to be paired.}
}
\value{
The function outputs a single non-negative integer that is uniquely associated with that unordered pair.
}
\source{
Szudzik, M. (2006): An Elegant Pairing Function. Wolfram Science Conference NKS 2006.

}
\details{With two of non-negative integers x and y as an input, the pairing is computed as:


\deqn{elegantPairing(x, y) = (x * y) + floor ((|x - y| - 1)^2) / 4)}


The function is commutative. x and y have to be non-negative integers.}

\examples{

ElegantPairingInt(2, 3)

}

\seealso{
  
  \code{\link{ElegantPairingVec}}
}

