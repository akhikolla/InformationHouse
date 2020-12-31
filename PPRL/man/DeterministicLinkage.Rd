\encoding{utf8}


\name{DeterministicLinkage}
\alias{DeterministicLinkage}
\alias{Record Linkage}

\title{Deterministic Record Linkage}
\usage{
DeterministicLinkage(IDA, dataA, IDB, dataB, blocking = NULL, similarity)
}
\description{
Deterministic Record Linkage of two data sets giving results enabling rule-based methods.
}


\arguments{

\item{IDA}{A character vector or integer vector containing the IDs of the first data.frame.}
\item{dataA}{A data.frame containing the data to be linked and all linking variables as specified in \code{\link{SelectBlockingFunction}} and \code{\link{SelectSimilarityFunction}}.}
\item{IDB}{A character vector or integer vector containing the IDs of the second data.frame.}
\item{dataB}{A data.frame containing the data to be linked and all linking variables as specified in \code{\link{SelectBlockingFunction}} and \code{\link{SelectSimilarityFunction}}.}
\item{blocking}{Optional blocking variables. See \code{\link{SelectBlockingFunction}}.}
\item{similarity}{Variables used for linking and their respective linkage methods as specified in \code{\link{SelectSimilarityFunction}}.}
}

\source{

Christen, P. (2012): Data Matching - Concepts and Techniques for Record Linkage, Entity Resolution, and Duplicate Detection. Springer.

Schnell, R., Bachteler, T., Reiher, J. (2004): A toolbox for record linkage. Austrian Journal of Statistics 33(1-2): 125-133.
}

\details{
To call the Deterministic Linkage function it is necessary to set up linking variables and methods. Using blocking variables is optional. Further options are available in \code{\link{SelectBlockingFunction}} and \code{\link{SelectSimilarityFunction}}.
}

\value{
A data.frame containing ID-pairs and the link status for each linking variable. This way, rules can be put into place allowing the user to classify links and non-links.
}

\examples{
# load test data
testFile <- file.path(path.package("PPRL"), "extdata/testdata.csv")
testData <- read.csv(testFile, head = FALSE, sep = "\t",
  colClasses = "character")

# define year of birth (V3) as blocking variable
bl <- SelectBlockingFunction("V3", "V3", method = "exact")

# Select first name and last name as linking variables,
# to be linked using the soundex phonetic (first name)
# and exact matching (last name)
l1 <- SelectSimilarityFunction("V7", "V7", method = "Soundex")
l2 <- SelectSimilarityFunction("V8", "V8", method = "exact")

# Link the data as specified in bl and l1/l2
# (in this small example data is linked to itself)
res <- DeterministicLinkage(testData$V1, testData,
  testData$V1, testData, similarity = c(l1, l2), blocking = bl)

}
\seealso{

  \code{\link{ProbabilisticLinkage}},
  \code{\link{SelectBlockingFunction}},
  \code{\link{SelectSimilarityFunction}},
  \code{\link{StandardizeString}}
}
