\name{makeFuzzy}
\alias{makeFuzzy}

\title{
Fuzzifying crisp-set data
}
\description{
The \code{makeFuzzy} function fuzzifies crisp-set data to a customizable degree. 
}
\usage{
makeFuzzy(x, fuzzvalues = c(0, 0.05, 0.1), ...)
}

\arguments{
  \item{x}{Data frame, matrix, or \code{configTable} featuring crisp-set (binary) factors with values 1 and 0 only.}
  \item{fuzzvalues}{Numeric vector of values from the interval [0,1].}
  \item{\dots}{Additional arguments are passed to \code{\link{configTable}}.}
}

\details{In combination with \code{allCombs}, \code{full.ct} and \code{selectCases}, \code{makeFuzzy} is useful for simulating fuzzy-set data, which are needed for inverse search trials benchmarking the output of \code{cna}. \code{makeFuzzy} transforms a data frame or \code{configTable} \code{x} consisting of crisp-set (binary) factors into a fuzzy-set \code{configTable} by adding values selected at random from the argument \code{fuzzvalues} to the 0's and subtracting them from the 1's in \code{x}. \code{fuzzvalues} is a numeric vector of values from the interval [0,1].

\code{selectCases} can be used before and \code{selectCases1} after the fuzzification to select those configurations that are compatible with a given data generating causal structure (see examples below).
}

\value{
A \code{configTable} of type "fs".
}

\seealso{\code{\link{selectCases}}, \code{\link{allCombs}}, \code{\link{full.ct}}, \code{\link{configTable}}, \code{\link{cna}}, \code{\link{ct2df}}, \code{\link{condition}}
}

\examples{
# Fuzzify a crisp-set (binary) 6x3 matrix with default fuzzvalues.
X <- matrix(sample(0:1, 18, replace = TRUE), 6)
makeFuzzy(X)

# ... and with customized fuzzvalues.
makeFuzzy(X, fuzzvalues = 0:5/10)
makeFuzzy(X, fuzzvalues = seq(0, 0.45, 0.01))

# First, generate crisp-set data comprising all configurations of 5 binary factors that 
# are compatible with the causal chain (A*b + a*B <-> C)*(C*d + c*D <-> E) and, 
# second, fuzzify those crisp-set data. 
dat1 <- full.ct(5)
dat2 <- selectCases("(A*b + a*B <-> C)*(C*d + c*D <-> E)", dat1)
(dat3 <- makeFuzzy(dat2, fuzzvalues = seq(0, 0.45, 0.01)))
condition("(A*b + a*B <-> C)*(C*d + c*D <-> E)", dat3)

# First, generate all configurations of 5 dichotomous factors that are compatible with
# the causal chain (A*b + a*B <-> C)*(C*d + c*D <-> E) and, second, fuzzify those 
# crisp-set data.
dat1 <- full.ct(5)
dat2 <- selectCases("(A*b + a*B <-> C)*(C*d + c*D <-> E)", dat1)
(dat3 <- makeFuzzy(dat2, fuzzvalues = seq(0, 0.45, 0.01)))
condition("(A*b + a*B <-> C)*(C*d + c*D <-> E)", dat3)

# Inverse search for the data generating causal structure A*b + a*B + C*D <-> E from
# fuzzy-set data with non-perfect consistency and coverage scores.
dat1 <- full.ct(5)
set.seed(55)
dat2 <- makeFuzzy(dat1, fuzzvalues = 0:4/10)
dat3 <- selectCases1("A*b + a*B + C*D <-> E", con = .8, cov = .8, dat2)
fscna(dat3, ordering = list("E"), strict = TRUE, con = .8, cov = .8)
}
