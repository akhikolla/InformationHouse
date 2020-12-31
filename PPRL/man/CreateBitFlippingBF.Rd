\encoding{utf8}

\name{CreateBitFlippingBF}

\alias{CreateBitFlippingBF}
\alias{RRT}

\alias{Permanent Randomized Response}

\title{Permanent Randomized Response Bloom Filters}

\usage{
CreateBitFlippingBF(data, password, f)
}
\description{
Applies Permanent Randomized Response to flip bits of the bit vectors given.
}
\arguments{
  \item{data}{a data.frame containing the IDs in the first column and bit vectors created by any Bloom filter-based method in the second column.}
  \item{password}{a string to seed the random bit flipping.}
  \item{f}{a numeric between 0 and 1 giving the probability of flipping a bit.}
}
\details{
The randomized response technique is used on each
bit position B[i] of a Bloom filter B. B[i] is set to one or zero with a probability of \eqn{1/2 * f} for each outcome. The bit position remains unchanged with a probability of \eqn{1 - f}, where  \eqn{ 0 <= f <= 1}.}

\value{
A data.frame containing IDs and the corresponding bit vector.
}
\source{
Schnell, R., Borgs, C. (2016): Randomized Response and Balanced Bloom Filters for Privacy Preserving Record Linkage. IEEE International Conference on Data Mining (ICDM 2016), Barcelona.

}

\examples{
# Load test data
testFile <- file.path(path.package("PPRL"), "extdata/testdata.csv")
testData <-read.csv(testFile, head = FALSE,
  sep = "\t", colClasses = "character")

## Encode data into Bloom Filters
BF <- CreateBF(ID = testData$V1, data = testData$V7,
  k = 20, padding = 1, q = 2, l = 1000,
  password = "(H]$6Uh*-Z204q" )

# Create Permanent Randomized Response Bloom Filter
RR <- CreateBitFlippingBF(BF, password = "l+kfdj1J", f = 0.1)

}
\seealso{
  
  \code{\link{CreateBF}},
  \code{\link{CreateCLK}},
  \code{\link{StandardizeString}}
}
