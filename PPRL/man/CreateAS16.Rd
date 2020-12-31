\encoding{utf8}


\name{CreateAS16}

\alias{CreateAS16}

\alias{Armknecht}
\title{Creating Records with Armknechts method create}
\usage{
CreateAS16(ID, data, password)
}
\description{
This method generates a new bit vector out of an existing Bloom Filter. Building and comparisons are both possible with \code{\link{CompareAS16}}.}

\arguments{
  \item{ID}{A character vector or integer vector containing the IDs of the second data.frame.}
  \item{data}{A character vector containing the original bit vectors created by any Bloom Filter-based method.}
  \item{password}{A string containing the password to be used for both "create" and "compare".}
}

\value{
A character vector containing bit vectors created as described in the original publication.
}
\source{
Armknecht, F., Schnell, R. (unpublished): Privacy Preserving Record Linkage Based on
Bloom Filters and Codes. Working Paper.
}


\examples{
# Load test data
testFile <- file.path(path.package("PPRL"), "extdata/testdata.csv")
testData <- read.csv(testFile, head = FALSE, sep = "\t",
  colClasses = "character")

# Create Bloom Filter
testData <- CreateBF(ID = testData$V1, testData$V7, k = 20, padding = 1,
  q = 2, l = 1000, password = "(H]$6Uh*-Z204q")

# Create the new Bloom Filter
testAS <- CreateAS16(testData$ID, testData$CLKs, password = "khäuds")


}
\seealso{
  
  \code{\link{CompareAS16}},
  \code{\link{CreateBF}},
  \code{\link{CreateCLK}}
}
