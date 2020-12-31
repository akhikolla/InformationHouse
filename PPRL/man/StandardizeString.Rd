\encoding{utf8}


\name{StandardizeString}

\alias{StandardizeString}
\alias{Preprocess}
\alias{Standardize}
\alias{Clean}
\title{Standardize String}

\description{
Preprocessing (cleaning) of strings prior to linkage.
}
\usage{
StandardizeString(strings)
}

\arguments{
  \item{strings}{A character vector of strings to be standardized.}
}
\details{
Strings are capitalized, letters are substituted as described below. Leading and trailing blanks are removed. Other non-ASCII characters are deleted.

\itemize{
\item Replace "Æ" with "AE"
\item Replace "æ" with "AE"
\item Replace "Ä" with "AE"
\item Replace "ä" with "AE"
\item Replace "Å" with "A"
\item Replace "å" with "A"
\item Replace "Â" with "A"
\item Replace "â" with "A"
\item Replace "À" with "A"
\item Replace "à" with "A"
\item Replace "Á" with "A"
\item Replace "á" with "A"
\item Replace "Ç" with "C"
\item Replace "Ç" with "C"
\item Replace "Ê" with "E"
\item Replace "ê" with "E"
\item Replace "È" with "E"
\item Replace "è" with "E"
\item Replace "É" with "E"
\item Replace "é" with "E"
\item Replace "Ï" with "I"
\item Replace "ï" with "I"
\item Replace "Î" with "I"
\item Replace "î" with "I"
\item Replace "Ì" with "I"
\item Replace "ì" with "I"
\item Replace "Í" with "I"
\item Replace "í" with "I"
\item Replace "Ö" with "OE"
\item Replace "ö" with "OE"
\item Replace "Ø" with "O"
\item Replace "ø" with "O"
\item Replace "Ô" with "O"
\item Replace "ô" with "O"
\item Replace "Ò" with "O"
\item Replace "ò" with "O"
\item Replace "Ó" with "O"
\item Replace "ó" with "O"
\item Replace "ß" with "SS"
\item Replace "Ş" with "S"
\item Replace "ş" with "S"
\item Replace "ü" with "UE"
\item Replace "Ü" with "UE"
\item Replace "Ů" with "U"
\item Replace "Û" with "U"
\item Replace "û" with "U"
\item Replace "Ù" with "U"
\item Replace "ù" with "U"}
}
\value{
Returns a character vector with standardized strings.
}



\examples{
strings = c("Päter", " Jürgen", " Roß")
StandardizeString(strings)

}

