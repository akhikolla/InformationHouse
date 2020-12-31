\encoding{utf8}


\name{CreateRecordLevelBF}
\alias{CreateRecordLevelBF}

\alias{RLBF}
\alias{Record Level Bloom Filter}
\alias{Durham}
\title{Record Level Bloom Filter (RLBF) Encoding}

\usage{
CreateRecordLevelBF(ID, data, password, lenRLBF = 1000, k = 20,
                       padding = as.integer(c(0)),
                       qgram = as.integer(c(2)),
                       lenBloom = as.integer(c(500)),
                       method = "StaticUniform",
                       weigths = as.numeric(c(1)))
}
\description{Creates Record Level Bloom filters, combining single Bloom filters into a singular bit vector.
}

\arguments{
\item{ID}{a character vector or integer vector containing the IDs of the data.frame.}
  \item{data}{a character vector containing the data to be encoded. Make sure the input vectors are not factors.}
  \item{password}{a string used as a password for the random hashing of the q-grams and the shuffling.}
  \item{lenRLBF}{length of the final Bloom filter.}
  \item{lenBloom}{an integer vector containing the length of the first level Bloom filters which are to be combined. For the methods "StaticUniform" and "StaticWeighted", a single integer is required, since all original Bloom filters will have the same length.}
  \item{k}{number of bit positions set to one for each q-gram.}
  \item{padding}{integer (0 or 1) indicating if string padding is to be used.}
  \item{qgram}{integer (1 or 2) indicating whether to split the input strings into bigrams (q = 2) or unigrams (q = 1).}
  \item{method}{any of either "StaticUniform", "StaticWeigthed", "DynamicUniform" or "DynamicWeighted" (see details).}
   \item{weigths}{weigths are used for the "StaticWeighted" and "DynamicWeighted" methods. The weights vector must have the same length as number of columns in the input data. The sum of the weights must be  1.}
}
\value{
A data.frame containing IDs and the corresponding bit vector.
}
\details{Single Bloom filters are first built for every variable in the input data.frame. Combining the Bloom filters is done by sampling a set fraction of the original Bloom filters and concatenating the samples. The result is then shuffled. The sampling can be done using four different weighting methods:
\enumerate{
   \item StaticUniform
   \item StaticWeighted
   \item DynamicUniform
   \item DynamicWeighted
}
Details are described in the original publication.}

\source{
Durham, E. A. (2012). A framework for accurate, efficient
private record linkage. Dissertation. Vanderbilt University.
}

\examples{
# Load test data
testFile <- file.path(path.package("PPRL"), "extdata/testdata.csv")
\dontrun{
testData <- read.csv(testFile, head = FALSE, sep = "\t",
  colClasses = "character")

# StaticUniform
RLBF <- CreateRecordLevelBF(ID = testData$V1,
  data = testData[, c(2, 3, 7, 8)],
  lenRLBF = 1000, k = 20,
  padding = c(0, 0, 1, 1), qgram = c(1, 1, 2, 2),
  lenBloom = 500,
  password = c("(H]$6Uh*-Z204q", "lkjg", "klh", "KJHkälk5"),
  method = "StaticUniform")

# StaticWeigthed
RLBF <- CreateRecordLevelBF(ID = testData$V1,
  data = testData[, c(2, 3, 7, 8)],
  lenRLBF = 1000, k = 20,
  padding = c(0, 0, 1, 1), qgram = c(1, 1, 2, 2),
  lenBloom = 500,
  password = c("(H]$6Uh*-Z204q", "lkjg", "klh", "KJHkälk5"),
  method = "StaticWeigthed", weigths = c(0.1, 0.1, 0.5, 0.3))

# DynamicUniform
RLBF <- CreateRecordLevelBF(ID = testData$V1,
  data = testData[, c(2, 3, 7, 8)],
  lenRLBF = 1000, k = 20,
  padding = c(0, 0, 1, 1), qgram = c(1, 1, 2, 2),
  lenBloom = c(300, 400, 550, 500),
  password = c("(H]$6Uh*-Z204q", "lkjg", "klh", "KJHkälk5"),
  method = "DynamicUniform")

# DynamicWeigthed
RLBF <- CreateRecordLevelBF(ID = testData$V1,
  data = testData[, c(2, 3, 7, 8)],
  lenRLBF = 1000, k = 20,
  padding = c(0, 0, 1, 1), qgram = c(1, 1, 2, 2),
  lenBloom = c(300, 400, 550, 500),
  password = c("(H]$6Uh*-Z204q", "lkjg", "klh", "KJHkälk5"),
  method = "DynamicWeigthed", weigths = c(0.1, 0.1, 0.5, 0.3))
}
}
\seealso{
  
  \code{\link{CreateBF}},
  \code{\link{CreateCLK}},
}

