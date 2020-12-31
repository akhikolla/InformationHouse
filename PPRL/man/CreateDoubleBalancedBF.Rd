\encoding{utf8}

\name{CreateDoubleBalancedBF}

\alias{CreateDoubleBalancedBF}
\alias{Double Balanced Bloom Filter}


\title{Double Balanced Bloomfilter Encoding}

\usage{
CreateDoubleBalancedBF(ID, data, password)
}
\description{
Double balanced Bloom Filter are created by first creating balanced Bloom Filters, see \code{\link{CreateBalancedBF}}, negating the whole data set and shuffling each Bloom Filter.
}
\arguments{
  \item{ID}{A character vector containing the ID. The ID vector must have the same size as the number of rows of data.}
  \item{data}{CLKs as created by \code{\link{CreateCLK}} or \code{\link{CreateBF}}.}
  \item{password}{A string to encode the routines.}
}
\value{
A data.frame containing IDs and the corresponding double balanced bit vector.
}


\references{
Schnell, R. (2017): Recent Developments in Bloom Filter-based Methods for Privacy-preserving Record Linkage. Curtin Institute for Computation, Curtin University, Perth, 12.9.2017.
}

\examples{
# Load test data
testFile <- file.path(path.package("PPRL"), "extdata/testdata.csv")
testData <-read.csv(testFile, head = FALSE, sep = "\t",
  colClasses = "character")

# Create bit vectors e.g. with CreateBF
testData <- CreateBF(ID = testData$V1,
  testData$V7, k = 20, padding = 1, q = 2,
  l = 1000, password = "(H]$6Uh*-Z204q")

# Create Double Balanced Bloom Filters
DBB <- CreateDoubleBalancedBF(ID = testData$ID, data = testData$CLKs,
  password = "hdayfkgh")

}
\seealso{

  \code{\link{CreateBalancedBF}},
  \code{\link{CreateBF}},
  \code{\link{CreateCLK}},
  \code{\link{StandardizeString}}
}
