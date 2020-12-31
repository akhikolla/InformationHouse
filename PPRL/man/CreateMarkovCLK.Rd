\encoding{utf8}


\name{CreateMarkovCLK}
\alias{CreateMarkovCLK}
\alias{Markov Chains}

\title{
Create CLKs with Markov Chain-based transition matrix
}
\description{
Builds CLKs encoding additional bigrams based on the transition probabilities as estimated by a Markov Chain.
}
\usage{
CreateMarkovCLK(ID, data, password, markovTable, k1 = 20, k2 = 4,
  padding = as.integer(c(0)), qgram = as.integer(c(2)),
  lenBloom = 1000,  includeOriginalBigram = TRUE, v = FALSE)}
\arguments{
\item{ID}{a character vector or integer vector containing the IDs of the data.frame.}
  \item{data}{a data.frame containing the data to be encoded. Make sure the input vectors are not factors.}
  \item{password}{a character vector with a password for each column of the input data.frame for the random hashing of the q-grams.}
  \item{markovTable}{a numeric matrix containing the transition probabilities for all bigrams possible.}

  \item{k1}{number of bit positions set to one for each bigram.}
  \item{k2}{number of additional bigrams drawn for each original bigram.}

  \item{padding}{integer vector (0 or 1) indicating if string padding is to be used on the columns of the input. The padding vector must have the same size as the number of columns of the input data.}
  \item{qgram}{integer vector (1 or 2) indicating whether to split the input strings into bigrams (q = 2) or unigrams (q = 1). The qgram vector must have the same size as the number of columns of the input data.}
  \item{lenBloom}{desired length of the final Bloom filter in bits.}
  \item{includeOriginalBigram}{by default, the original bigram is encoded together with the additional bigrams. Set this to \code{FALSE} to include only the additional bigrams to further increase the security.}
   \item{v}{verbose.}
}

\details{
A transition matrix for all possible bigrams is built using a function to fit a markov chain distribution using a Laplacian smoother. A transition matrix built for bigrams using the NC Voter Data is included in the package. For each original bigram in the data, \code{k2} new bigrams are drawn according to their follow-up probability as given by the transition matrix. The final bigram set is then encoded following  \code{\link{CreateCLK}}.
}
\value{
A data.frame containing IDs and the corresponding bit vector.
}
\references{
Schnell R., Borgs C. (2017): Using Markov Chains for Hardening Bloom Filter Encryptions against
Cryptographic Attacks in Privacy Preserving Record Linkage. German Record Linkage Center Working Paper.

Schnell, R. (2017): Recent Developments in Bloom Filter-based Methods for Privacy-preserving Record Linkage. Curtin Institute for Computation, Curtin University, Perth, 12.9.2017.

}

\seealso{

  \code{\link{CreateCLK}},
  \code{\link{StandardizeString}}
}


\examples{
# Load test data
testFile <- file.path(path.package("PPRL"), "extdata/testdata.csv")
testData <-read.csv(testFile, head = FALSE, sep = "\t",
  colClasses = "character")

\dontrun{
# Load example Markov chain matrix (created from NC Voter Data)
markovFile <-file.path(path.package("PPRL"), "extdata/TestMatrize.csv")
markovData <-read.csv(markovFile,  sep = " ",
  header = TRUE, check.names = FALSE)
markovData <- as.matrix(markovData)

# Create Markov CLK using
CLKMarkov <- CreateMarkovCLK(ID = testData$V1,
  data = testData[, c(2, 3, 7, 8)],
  markovTable = markovData,
  k1 = 20, k2 = 4, l = 1000,
  padding = c(0, 0, 1, 1),
  q = c(1, 2, 2, 2),
  password = c("(H]$6Uh*-Z204q", "lkjg", "klh", "KJHklk5"),
  includeOriginalBigram = TRUE)
}
}

