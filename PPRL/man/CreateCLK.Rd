\encoding{utf8}


\name{CreateCLK}

\alias{CreateCLK}
\alias{CLK}
\alias{CLKs}
\alias{additive BF}
\title{Cryptographic Long-term Keys (CLKs)}
\usage{
CreateCLK(ID, data, password, k = 20,
  padding =  as.integer(c(0)),
  qgram =  as.integer(c(2)), lenBloom = 1000)
}
\description{
Each column of the input data.frame is hashed into a single additive Bloom filter.
}
\arguments{
\item{ID}{A character vector or integer vector containing the IDs of the data.frame.}
  \item{data}{a data.frame containing the data to be encoded. Make sure the input vectors are not factors.}
  \item{password}{a character vector with a password for each column of the input data.frame for the random hashing of the q-grams.}
  \item{k}{number of bit positions set to one for each bigram.}

  \item{padding}{integer vector (0 or 1) indicating if string padding is to be used on the columns of the input. The padding vector must have the same size as the number of columns of the input data.}
  \item{qgram}{integer vector (1 or 2) indicating whether to split the input strings into bigrams (q = 2) or unigrams (q = 1). The qgram vector must have the same size as the number of columns of the input data.}
  \item{lenBloom}{desired length of the final Bloom filter in bits.}
}
\value{
A data.frame containing IDs and the corresponding bit vector.
}

\source{
Schnell, R. (2014): An efficient Privacy-Preserving Record Linkage Technique for Administrative Data and Censuses. Journal of the International Association for Official Statistics (IAOS) 30: 263-270.
}


\examples{
# Load test data
testFile <- file.path(path.package("PPRL"), "extdata/testdata.csv")
testData <-read.csv(testFile, head = FALSE, sep = "\t",
  colClasses = "character")

## Encode data
CLK <- CreateCLK(ID = testData$V1,
  data = testData[, c(2, 3, 7, 8)],
  k = 20, padding = c(0, 0, 1, 1),
  q = c(1, 1, 2, 2), l = 1000,
  password = c("HUh4q", "lkjg", "klh", "Klk5"))

}
\seealso{

  \code{\link{CreateBF}},
  \code{\link{StandardizeString}}
}

