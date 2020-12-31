\encoding{utf8}

\name{CreateBalancedBF}

\alias{CreateBalancedBF}
\alias{BalancedBloomfilter}
\alias{Balanced Codes}
\alias{Green function}
\alias{Balanced Bloom Filter}

\title{Balanced Bloom Filter Encoding}
\usage{
CreateBalancedBF(ID, data, password)
}
\description{

Creates CLKs with constant Hamming weights by adding a negated copy of the binary input vector which is then permutated.
}
\arguments{
\item{ID}{A character vector or integer vector containing the IDs of the data.frame.}
  \item{data}{Bit vectors as created by any Bloom filter-based method.}
  \item{password}{a string used as a password for the random permutation.}
}

\value{
A data.frame containing IDs and the corresponding Balanced Bloom Filter.
}


\references{


Berger, J. M. (1961): A Note on Error Detection Codes for
Asymmetric Channels. In: Information and Control 4: 68–73.


Knuth, Donald E. (1986): Efficient Balanced Codes. In: IEEE
Transactions on Information Theory IT-32 (1): 51–53.

Schnell, R., Borgs, C. (2016): Randomized Response and Balanced Bloom Filters for Privacy Preserving Record Linkage. IEEE International Conference on Data Mining (ICDM 2016), Barcelona.
}

\examples{
# Load test data
testFile <- file.path(path.package("PPRL"), "extdata/testdata.csv")
testData <- read.csv(testFile, head = FALSE, sep = "\t",
  colClasses = "character")

# Create bit vectors e.g. with CreateBF
testData <- CreateBF(ID = testData$V1,
  testData$V7, k = 20, padding = 1, q = 2,
  l = 1000, password = "(H]$6Uh*-Z204q")

# Create Balanced Bloom Filters
BB <- CreateBalancedBF(ID = testData$ID, data = testData$CLKs,
  password = "hdayfkgh")

}
\seealso{

  \code{\link{CreateBF}},
  \code{\link{CreateBitFlippingBF}},
  \code{\link{CreateCLK}},
  \code{\link{CreateDoubleBalancedBF}},
  \code{\link{CreateEnsembleCLK}},
  \code{\link{CreateMarkovCLK}},
  \code{\link{CreateRecordLevelBF}},
  \code{\link{StandardizeString}}
}
