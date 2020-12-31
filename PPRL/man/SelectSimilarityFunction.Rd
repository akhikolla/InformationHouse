\encoding{utf8}


\name{SelectSimilarityFunction}
\alias{SelectSimilarityFunction}
\alias{Similarity}
\alias{Linkage}

\title{
Select Similarity Function for Linkage
}

\description{
To call \code{\link{DeterministicLinkage}} or \code{\link{ProbabilisticLinkage}}
 it is mandatory to select a similarity function for each variable. Each element of the setup contains the two variable names and the method. For some methods further informations can be entered.
}

\usage{
SelectSimilarityFunction(variable1, variable2,
  method = "jw",
  ind_c0 = FALSE, ind_c1 = FALSE,
  m = 0.9, u = 0.1, p = 0.05, epsilon = 0.0004,
  lower = 0.0, upper = 0.0,
  threshold = 0.85, jaroWeightFactor = 1.0, lenNgram = 2)
}

\arguments{
  \item{variable1}{name of linking variable 1 in the data.frame. The column must be of type character, numeric or integer,  containing the data to be merged. The data vector must have the same length as the ID vector.}
  \item{variable2}{name of linking variable 2 in the data.frame. The column must be of type character, numeric or integer,  containing the data to be merged. The data vector must have the same length as the ID vector.}
  \item{method}{linking method. Possible values are:
  \itemize{
\item 'exact' = Exact matching
\item 'exactCL' = Exact matching using capital letters
\item 'LCS' = Longest Common Subsequence
\item 'lv' = Levenshtein distance
\item 'dl' = Damerau Levenshtein distance
\item 'jaro' = Jaro similarity
\item 'jw' = Jaro-Winkler similarity
\item 'jw2' = Modified Jaro-Winkler similarity
\item 'ngram' = n-gram similarity
\item 'Gcp' = German census phonetic (Baystat)
\item 'Reth' = Reth-Schek (IBM) phonetic
\item 'Soundex' = Soundex phonetic
\item 'Metaphone' = Metaphone phonetic
\item 'DoubleMetaphone' = Double Metaphone phonetic
}
  }

\item{ind_c0}{Only used for jw2.

Increase the probability of a match when the number of matched
              characters is large. This option allows for a little more
              tolerance when the strings are large. It is not an appropriate
              test when comparing fixed length fields such as phone and
              social security numbers. A nonzero value indicates the option is deactivated.
}

\item{ind_c1}{Only used for jw2.

All lower case characters are converted to upper case prior
              to the comparison. Disabling this feature means that the lower
              case string "code" will not be recognized as the same as the
              upper case string "CODE". Also, the adjustment for similar
              characters section only applies to uppercase characters. A nonzero value indicates the option is deactivated.

}
\item{m}{Initial m value for the EM algorithm. Only used when linking using \code{\link{ProbabilisticLinkage}}. 0.0 < m < 1.0.}
\item{u}{Initial u value for the EM algorithm. Only used when linking using \code{\link{ProbabilisticLinkage}}. 0.0 < u < 1.0.}
\item{p}{Initial p value for the EM algorithm. Only used when linking using \code{\link{ProbabilisticLinkage}}. 0.0 < u < 1.0.}
\item{epsilon}{epsilon is a stop criterum for the EM algorithm. The EM algorithm can be terminated when relative change of likelihood logarithm is less than epsilon. Only used when linking using \code{\link{ProbabilisticLinkage}}.}
\item{lower}{Matches lower than 'lower' are classified as non-match. Everything between 'lower' and 'upper' is classified as possible match. Only used when linking using \code{\link{ProbabilisticLinkage}}.}
\item{upper}{Matches higher than 'upper' are classified as match. Everything between 'lower' and 'upper' is classified as possible match. Only used when linking using \code{\link{ProbabilisticLinkage}}.}
\item{threshold}{If using string similarities: Outputs only matches above the similarity threshold value. If using string distances: Outputs only matches below the set threshold distance.}
\item{jaroWeightFactor}{By the Jaro weight adjustment the matching weight is adjusted according to the degree of similarity between the
variable values. The weight factor which determines the Jaro adjusted matching weight. Only used when linking using \code{\link{ProbabilisticLinkage}}.}
\item{lenNgram}{Length of ngrams. Only used for the method ngram. Length of ngrams must be between 1 and 4.}
}

\value{
Calling the function will not return anything.
}

\references{

Christen, P. (2012): Data Matching - Concepts and Techniques for Record Linkage, Entity Resolution, and Duplicate Detection. Springer.

Schnell, R., Bachteler, T., Reiher, J. (2004): A toolbox for record linkage. Austrian Journal of Statistics 33(1-2): 125-133.

Winkler, W. E. (1988): Using the EM algorithm for weight computation in the Fellegi-Sunter model of record linkage. Proceedings of the Section on Survey Research Methods Vol. 667, American Statistical Association: 671.

}



\seealso{

  \code{\link{DeterministicLinkage}},
  \code{\link{ProbabilisticLinkage}},
  \code{\link{SelectBlockingFunction}},
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
# to be linked using the jaro-winkler (first name)
# and exact matching (last name)
l1 <- SelectSimilarityFunction("V7","V7", method = "jw",
  ind_c0 = FALSE, ind_c1 = FALSE , m = 0.9, u = 0.1,
  lower = 0.0, upper = 0.0)
l2 <- SelectSimilarityFunction("V8","V8", method = "exact")

# Link the data as specified in bl and l1/l2
# (in this small example data is linked to itself)
res <- ProbabilisticLinkage(testData$V1, testData,
  testData$V1, testData, similarity = c(l1, l2), blocking = bl)

}

