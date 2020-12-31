\encoding{utf8}


\name{CreateEnsembleCLK}
\alias{CreateEnsembleCLK}

\alias{Ensemble CLK}
\alias{majority}
\alias{majority rule}

\title{
Combine multiple independent CLKs using a simple majority rule
}
\description{
Creates multiple CLKs which are combined using a simple majority rule.
}

\usage{
CreateEnsembleCLK(ID, data, password, NumberOfCLK = 1 , k = 20,
  padding = as.integer(c(0)), qgram = as.integer(c(2)),
  lenBloom = 1000)
}

\arguments{
\item{ID}{A character vector or integer vector containing the IDs of the data.frame.}
  \item{data}{a data.frame containing the data to be encoded. Make sure the input vectors are not factors.}
  \item{password}{a character vector with a password for each column of the input data.frame for the random hashing of the q-grams.}
  \item{NumberOfCLK}{number of independent CLKs to be built.}
  \item{k}{number of bit positions set to one for each bigram.}
  \item{padding}{integer vector (0 or 1) indicating if string padding is to be used on the columns of the input. The padding vector must have the same size as the number of columns of the input data.}
  \item{qgram}{integer vector (1 or 2) indicating whether to split the input strings into bigrams (q = 2) or unigrams (q = 1). The qgram vector must have the same size as the number of columns of the input data.}
  \item{lenBloom}{desired length of the final Bloom filter in bits.}
}
\details{
Creates a set number of independent CLKs for each record of the input data.frame and combines them using  a simple majority rule. The bit positions \eqn{B[i]} in the final CLK of the length \eqn{B} are set to \eqn{B[i] = 1} if more than half of the independent CLKs bit positions \eqn{B[i]} have a value of one. Otherwise the final bit position is zero.
}
\value{
A data.frame containing IDs and the corresponding ensemble bit vector.
}
\references{
Kuncheva, L. (2014): Combining Pattern Classifiers: Methods and Algorithms. Wiley.

}

\seealso{

  \code{\link{CreateBF}},
  \code{\link{CreateCLK}},
  \code{\link{StandardizeString}}
}


\examples{
# Load test data
testFile <- file.path(path.package("PPRL"), "extdata/testdata.csv")
testData <-read.csv(testFile, head = FALSE, sep = "\t",
  colClasses = "character")
\dontrun{
# Create Ensemble CLK
EnsembleCLK <- CreateEnsembleCLK(ID = testData$V1,
  data = testData[, c(2, 3, 7, 8)],
  k = 20, padding = c(0, 0, 1, 1),
  q = c(1, 2, 2, 2), l = 1000,
  password = c("HUh4q", "lkjg", "klh", "Klk5"),
  NumberOfCLK = 5)
}
}

