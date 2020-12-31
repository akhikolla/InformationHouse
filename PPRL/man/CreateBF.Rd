\encoding{utf8}

\name{CreateBF}

\alias{CreateBF}

\alias{BF}
\alias{Bloomfilter}
\alias{Bloom Filter}

\title{Bloom Filter Encoding}
\usage{
CreateBF(ID, data, password, k = 20, padding = 1, qgram = 2, lenBloom = 1000)
}
\description{
Creates Bloom filters for each row of the input data by splitting the input into q-grams which are hashed into a bit vector.
}
\arguments{
\item{ID}{a character vector or integer vector containing the IDs of the data.frame.}
  \item{data}{a character vector containing the data to be encoded. Make sure the input vectors are not factors.}
  \item{password}{a string used as a password for the random hashing of the q-grams.}
  \item{k}{number of bit positions set to one for each bigram.}
  \item{padding}{integer (0 or 1) indicating if string padding is to be used.}
  \item{qgram}{integer (1 or 2) indicating whether to split the input strings into bigrams (q = 2) or unigrams (q = 1).}
  \item{lenBloom}{desired length of the Bloom filter in bits.}
}

\value{
A data.frame containing IDs and the corresponding bit vector.
}
\source{
Schnell, R., Bachteler, T., Reiher, J. (2009): Privacy-preserving record linkage using Bloom filters. BMC Medical Informatics and Decision Making 9: 41.
}

\examples{
# Load test data
testFile <- file.path(path.package("PPRL"), "extdata/testdata.csv")
testData <- read.csv(testFile, head = FALSE, sep = "\t",
  colClasses = "character")

# Encode data
BF <- CreateBF(ID = testData$V1, data = testData$V7,
  k = 20, padding = 1, q = 2, l = 1000,
  password = "(H]$6Uh*-Z204q")

}
\seealso{

  \code{\link{CreateCLK}},
  \code{\link{StandardizeString}}
}
