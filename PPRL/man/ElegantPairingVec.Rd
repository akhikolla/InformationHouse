\encoding{utf8}

\name{ElegantPairingVec}
\alias{ElegantPairingVec}
\alias{Pairing}
\alias{unordered pairing}

\title{Simple Pairing Function}
\usage{
ElegantPairingVec (ID, data)
}

\description{
Unordered Pairing Function creating a new unique integer from two input integers in a \code{\link{data.frame}}.
}

\arguments{
\item{ID}{A character vector or integer vector containing the IDs of the data.frame.}
  \item{data}{a \code{\link{data.frame}} consisting of two columns containing the integers on which the pairing function is to be applied.}
}
\value{
A data.frame containing IDs and the computed integer.
}
\details{With two of non-negative integers x and y as an input, the pairing is computed as:


\deqn{elegantPairing(x, y) = (x * y) + floor ((|x - y| - 1)^2) / 4)}


The function is commutative. x and y have to be non-negative integers. The function outputs a single
non-negative integer that is uniquely associated with that unordered pair.}

\source{
Szudzik, M. (2006): An Elegant Pairing Function. Wolfram Science Conference NKS 2006.

}


\examples{
# Load test data
testFile <- file.path(path.package("PPRL"), "extdata/testdata.csv")
testData <- read.csv(testFile, head = FALSE, sep = "\t",
  colClasses = "character")

# Create numeric data frame of day and month of birth
dataInt <- data.frame(as.integer(testData$V4), as.integer(testData$V5))

# Use unordered pairing on day and month
res <- ElegantPairingVec(testData$V1, dataInt)

}
\seealso{
  
  \code{\link{ElegantPairingInt}}
}
