\encoding{utf8}


\name{CreateALC}

\alias{CreateALC}
\alias{ALC}
\alias{Soundex}

\title{
Anonymous Linkage Codes (ALCs)
}
\description{
Creates ALCs from clear-text data by creating soundex phonetics for first and last names and concatenating all other identifiers. The resulting code is encrypted using SHA-2. The user can decide on which columns the soundex phonetic is applied.
}
\usage{
CreateALC(ID, data, soundex, password)
}

\arguments{
\item{ID}{A character vector or integer vector containing the IDs of the data.frame.}
  \item{data}{a data.frame containing the data to be encoded.}
  \item{soundex}{a binary vector with one element for each input column, indicating whether soundex is to be used. 1 = soundex is used, 0 = soundex is not used. The soundex vector must have the same length as the number of columns the data.frame.}
  \item{password}{a string used as a password for the HMAC.}
}
\source{
Herzog, T. N., Scheuren, F. J., Winkler, W. E. (2007): Data Quality and Record Linkage Techniques. Springer.

}

\value{
A data.frame containing IDs and the corresponding Anonymous Linkage Codes.
}


\examples{
# Load test data
testFile <- file.path(path.package("PPRL"), "extdata/testdata.csv")
testData <- read.csv(testFile, head = FALSE, sep = "\t",
  colClasses = "character")

# Encrypt data, use Soundex for names
res <- CreateALC(ID = testData$V1,   
  data = testData[, c(2, 3, 7, 8)],
  soundex = c(0, 0, 1, 1),  
  password = "$6Uh*-Z204q")
  
}

\seealso{
  
  \code{\link{Create581}},
  \code{\link{StandardizeString}}
}

