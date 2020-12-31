\encoding{utf8}

\name{Create581}

\alias{Create581}
\alias{ESL}
\alias{581}
\alias{Australian}

\title{Create Encrypted Statistical Linkage Keys}

\usage{
Create581(ID, data, code, password)
}
\description{
Creates ESLs (also known as 581-Keys), which are the hashed combination of the full date of birth and sex and subsets of first and last names. 
}
\arguments{
\item{ID}{a character vector or integer vector containing the IDs of the data.frame.}
  \item{data}{a data.frame containig the data to be encoded.}
  \item{code}{a list indicating  how data is to be encoded for each column. The list must have the same length as the number of columns of the data.frame to be encrypted.}
  \item{password}{a string used as a password for the HMAC.}
}
\details{
The original implementation uses the second and third position of the first name, the second, third and fifth position of the last name and full date of birth and sex as an input for an HMAC. This would be akin to using \code{code = list(c(2, 3), c(2, 3, 5), 0, 0)}. In this implementation, the positions of the subsets can be customized.
}
\source{
Karmel, R., Anderson, P., Gibson, D., Peut, A., Duckett, S., Wells, Y. (2010): Empirical aspects of record linkage across multiple data sets using statistical linkage keys: the experience of the PIAC cohort study. BMC Health Services Research 41(10).
}

\value{
A data.frame containing IDs and the corresponding Encrypted Statistical Linkage Keys.
}

\examples{
# Load test data
testFile <- file.path(path.package("PPRL"), "extdata/testdata.csv")
testData <- read.csv(testFile, head = FALSE, sep = "\t",
  colClasses = "character")

# Encrypt data
res <- Create581(ID = testData$V1,
  data = testData[, c(2, 3, 7, 8)],
  code = list(0, 0, c(2, 3), c(2, 3, 5)),
  password = "(H]$6Uh*-Z204q")

# Code: 0 means the whole string is used, 
# c(2, 3) means the second and third letter of the string is used

}
\seealso{
  
  \code{\link{CreateALC}},
  \code{\link{StandardizeString}}
}
