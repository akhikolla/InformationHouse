\encoding{utf8}


\name{SelectBlockingFunction}
\alias{SelectBlockingFunction}
\alias{blocking}
\alias{indexing}
\alias{partitioning}

\title{Select blocking method prior to linkage}

\description{
Before calling \code{\link{ProbabilisticLinkage}} or \code{\link{DeterministicLinkage}}, a blocking method can be selected. For each blocking variable desired, the function call has to be repeated.
}

\usage{
SelectBlockingFunction(variable1, variable2, method)
}

\arguments{
  \item{variable1}{Column name of blocking variable 1.}
  \item{variable2}{Column name of blocking variable 2.}
  \item{method}{Desired blocking method. Possible values are \code{'exact'} and \code{'exactCL'}.}
}

\details{
\describe{The following methods are available for blocking:

  \item{\code{'exact'}}{Simple exact blocking. All records with the same values for the blocking variable create a block. Searching for links is only done within these blocks.}
  \item{\code{'exactCL'}}{The same as \code{'exact'}. Only works with strings; all caracters are capitalised.}
%  \item{\code{'sortedNeighbourhood'}}{Sorted nearest neighbourhood blocking (see references); to be used exclusively with Bloom filter-based encryptions.}
%  \item{\code{'canopy'}}{Canopy clustering (see references); to be used exclusively with Bloom filter-based encryptions.}
%  \item{\code{'mbt'}}{Multibit tree blocking (see references); to be used exclusively with Bloom filter-based encryptions.}
}

}

\references{
Christen, P. (2012): Data Matching - Concepts and Techniques for Record Linkage, Entity Resolution, and Duplicate Detection. Springer.
}


\seealso{

  \code{\link{DeterministicLinkage}},
  \code{\link{ProbabilisticLinkage}},
  \code{\link{SelectSimilarityFunction}},
  \code{\link{StandardizeString}}
}


\examples{
# load test data
testFile <- file.path(path.package("PPRL"), "extdata/testdata.csv")
testData <- read.csv(testFile, head = FALSE, sep = "\t",
  colClasses = "character")

# define year of birth (V3) as blocking variable
bl <- SelectBlockingFunction("V3", "V3", method = "exact")

# Select first name and last name as linking variables,
# to be linked using the Jaro-Winkler similarity measure (first name)
# and levenshtein distance (last name)
l1 <- SelectSimilarityFunction("V7", "V7", method = "jw")
l2 <- SelectSimilarityFunction("V8", "V8", method = "lv")

# Link the data as specified in bl and l1/l2
# (in this small example data is linked to itself)
res <- ProbabilisticLinkage(testData$V1, testData,
  testData$V1, testData, similarity = c(l1, l2), blocking = bl)

}

