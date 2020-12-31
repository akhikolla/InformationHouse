\name{full.ct}
\alias{full.ct}
\alias{full.ct.default}
\alias{full.ct.configTable}
\alias{full.ct.cti}

\title{
Generate all logically possible value configurations of a given set of factors
}
\description{
The function \code{full.ct} generates a \code{configTable} with all logically possible value configurations of the factors defined in the input \code{x}. It is more flexible than \code{\link{allCombs}}. 
\code{x} can be a \code{configTable}, a data frame, an integer, a list specifying the factors' value ranges, or a character vector featuring all admissible factor values.
}

\usage{
full.ct(x, ...)

\method{full.ct}{default}(x, type = c("cs", "mv", "fs"), ...)
\method{full.ct}{configTable}(x, ...)
\method{full.ct}{cti}(x, ...)
}

\arguments{
  \item{x}{A \code{configTable}, a data frame, a matrix, an integer, a list specifying the factors' value ranges,  or a character vector featuring all admissible factor values (see the details and examples below).}
  \item{type}{Character vector specifying the type of \code{x}: \code{"cs"} (crisp-set), \code{"mv"} (multi-value),  or \code{"fs"} (fuzzy-set); passed to \code{\link{configTable}}; only required if \code{x} is a data frame or matrix.}
  \item{\dots}{
  Further arguments passed to methods.
}
}

\details{
\code{full.ct} generates all logically possible value configurations of the factors defined in \code{x}, which can either be a character vector or an integer or a list or a data frame or a matrix. \itemize{
\item{If \code{x} is a character vector, it can be a condition of any of the three types of conditions, \emph{boolean}, \emph{atomic} or \emph{complex} (see \code{\link{condition}}). \code{x} must contain at least one factor. Factor names and admissible values are guessed from the Boolean formulas. If \code{x} contains multi-value factors, only those values are considered admissible that are explicitly contained in \code{x}. Accordingly, in case of multi-value factors, \code{full.ct} should be given the relevant factor definitions by means of a list (see below).
}
\item{If \code{x} is an integer and \code{<=26}, the output will be a logically complete configuration table of type \code{"cs"} with \code{x} factors. The first \code{x} capital letters of the alphabet will be used as the names of the factors.
}
\item{If \code{x} is a list, \code{x} is expected to have named elements each of which provides the factor names with corresponding vectors enumerating their admissible values (i.e. their value ranges). These values must be non-negative integers.
}
\item{If \code{x} is a \code{configTable}, data frame, or matrix, \code{colnames(x)} are interpreted as factor names and the rows as enumerating the admissible values (i.e. as value ranges). If \code{x} is a data frame or a matrix, \code{x} is first converted to a \code{\link{configTable}} (the function \code{configTable} is called with \code{type} as specified in \code{full.ct}), and the \code{configTable} method of \code{full.ct} is then applied to the result. The \code{configTable} method uses all factors and factor values occurring in the \code{configTable}. If \code{x} is of type \code{"fs"}, 0 and 1 are taken as the admissible values.
}
}
In combination with \code{selectCases}, \code{full.ct} is useful for simulating data, which are needed for inverse search trials benchmarking the output of \code{cna}. While \code{full.ct} generates the space of all logically possible configurations of the factors in an analyzed factor set, \code{selectCases} selects those configurations from this space that are compatible with a given data-generating causal structure (i.e. the ground truth), that is, it selects the empirically possible configurations.

The method for class "cti" is for internal use only.
}

\value{
A \code{\link{configTable}} of type \code{"cs"} or \code{"mv"} with the full enumeration of combinations of the factor values.
}

\seealso{
\code{\link{configTable}}, \code{\link{selectCases}}, \code{\link{allCombs}}
}

\examples{
# x is a character vector.
full.ct("A + B*c")
full.ct("A=1*C=3 + B=2*C=1 + A=3*B=1") 
full.ct(c("A + b*C", "a*D"))
full.ct("!A*-(B + c) + F")

# x is a data frame.
full.ct(d.educate)
full.ct(d.jobsecurity, type = "fs")
full.ct(d.pban, type = "mv")

# x is a configTable.
full.ct(csct(d.educate))
full.ct(fsct(d.jobsecurity))
full.ct(mvct(d.pban))

# x is an integer.
full.ct(6)

# x is a list.
full.ct(list(A = 0:1, B = 0:1, C = 0:1))  # cs
full.ct(list(A = 1:2, B = 0:1, C = 1:4))  # mv

# Simulating crisp-set data.
groundTruth.1 <- "(A*b + C*d <-> E)*(E*H + I*k <-> F)"
fullData <- full.ct(groundTruth.1)
idealData <- selectCases(groundTruth.1, fullData) 
# Introduce 20% data fragmentation.
fragData <- idealData[-sample(1:nrow(idealData), nrow(idealData)*0.2), ] 
# Introduce 10% random noise.
(realData <- rbind(ct2df(fullData[sample(1:nrow(fullData), nrow(fragData)*0.1), ]), 
  fragData))

# Simulating multi-value data.
\donttest{groundTruth.2 <- "(JO=3 + TS=1*PE=3 <-> ES=1)*(ES=1*HI=4 + IQ=2*KT=5 <-> FA=1)"
fullData <- full.ct(list(JO=1:3, TS=1:2, PE=1:3, ES=1:2, HI=1:4, IQ=1:5, KT=1:5, FA=1:2))
idealData <- selectCases(groundTruth.2, fullData) 
# Introduce 20% data fragmentation.
fragData <- idealData[-sample(1:nrow(idealData), nrow(idealData)*0.2), ] 
# Introduce 10% random noise.
(realData <- rbind(ct2df(fullData[sample(1:nrow(fullData), nrow(fragData)*0.1), ]), 
  fragData))
}
}
