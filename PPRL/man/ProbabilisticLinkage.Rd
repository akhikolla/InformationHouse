\encoding{utf8}

\name{ProbabilisticLinkage}
\alias{ProbabilisticLinkage}
\alias{MergeToolBox}
\alias{Merge Tool Box}
\alias{Probabilistic}

\title{Probabilistic Record Linkage}
\usage{
ProbabilisticLinkage(IDA, dataA, IDB, dataB,  blocking = NULL, similarity)
}

\description{
Probabilistic Record Linkage of two data sets using distance-based or probabilistic methods.
}

\arguments{

\item{IDA}{A character vector or integer vector containing the IDs of the first data.frame.}
\item{dataA}{A data.frame containing the data to be linked and all linking variables as specified in \code{\link{SelectBlockingFunction}} and \code{\link{SelectSimilarityFunction}}.}
\item{IDB}{A character vector or integer vector containing the IDs of the second data.frame.}
\item{dataB}{A data.frame containing the data to be linked and all linking variables as specified in \code{\link{SelectBlockingFunction}} and \code{\link{SelectSimilarityFunction}}.}
\item{blocking}{Optional blocking variables. See \code{\link{SelectBlockingFunction}}.}
\item{similarity}{Variables used for linking and their respective linkage methods as specified in \code{\link{SelectSimilarityFunction}}.}
}

\details{
To call the Probabilistic Linkage function it is necessary to set up linking variables and methods. Using blocking variables is optional. Further options are available in \code{\link{SelectBlockingFunction}} and \code{\link{SelectSimilarityFunction}}. Using this method, the Fellegi-Sunter model is used, with the EM algorithm estimating the weights (Winkler 1988).
}

\value{
A data.frame containing pairs of IDs, their corresponding similarity value and the match status as determined by the linkage procedure.
}

\source{

Christen, P. (2012): Data Matching - Concepts and Techniques for Record Linkage, Entity Resolution, and Duplicate Detection. Springer.

Schnell, R., Bachteler, T., Reiher, J. (2004): A toolbox for record linkage. Austrian Journal of Statistics 33(1-2): 125-133.

Winkler, W. E. (1988): Using the EM algorithm for weight computation in the Fellegi-Sunter model of record linkage. Proceedings of the Section on Survey Research Methods Vol. 667, American Statistical Association: 671.

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
\seealso{

  \code{\link{CreateBF}},
  \code{\link{CreateCLK}},
  \code{\link{DeterministicLinkage}},
  \code{\link{SelectBlockingFunction}},
  \code{\link{SelectSimilarityFunction}},
  \code{\link{StandardizeString}}

}
