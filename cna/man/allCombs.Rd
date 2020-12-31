\name{allCombs}
\alias{allCombs}

\title{
Generate all logically possible value configurations of a given set of factors
}

\description{
The function \code{allCombs} generates a data frame of all possible value configurations of \code{length(x)} factors, the first factor having \code{x[1]} values, the second \code{x[2]} values etc. The factors are labeled using capital letters. 
}

\usage{
allCombs(x)
}

\arguments{
  \item{x}{Integer vector with values >0}
}

\details{In combination with \code{selectCases}. \code{makeFuzzy}, and \code{is.submodel}, \code{allCombs} is useful for simulating data, which are needed for inverse search trials benchmarking the output of \code{cna}. In a nutshell, \code{allCombs} generates the space of all logically possible configurations of the factors in an analyzed factor set, \code{selectCases} selects those configurations from this space that are compatible with a given data-generating causal structure (i.e. the ground truth, which can be randomly generated using \code{randomConds}), \code{makeFuzzy} fuzzifies those data, and \code{is.submodel} checks whether the models returned by \code{cna} are true of the ground truth.

The \pkg{cna} package provides another function to the same effect, \code{\link{full.ct}}, which is more flexible than \code{allCombs}.
}

\value{
A data frame.
}

\seealso{\code{\link{selectCases}}, \code{\link{makeFuzzy}}, \code{\link{is.submodel}}, \code{\link{randomConds}}, \code{\link{full.ct}}}

\examples{
# Generate all logically possible configurations of 5 dichotomous factors named "A", "B",
# "C", "D", and "E". 
allCombs(c(2, 2, 2, 2, 2)) - 1
# allCombs(c(2, 2, 2, 2, 2)) generates the value space for values 1 and 2, but as it is
# conventional to use values 0 and 1 for Boolean factors, 1 must be subtracted from
# every value output by allCombs(c(2, 2, 2, 2, 2)) to yield a Boolean data frame.

# Generate all logically possible configurations of 5 multi-value factors named "A", "B",
# "C", "D", and "E", such that A can take on 3 values {1,2,3}, B 4 values {1,2,3,4},
# C 3 values etc.
dat0 <- allCombs(c(3, 4, 3, 5, 3))
head(dat0)
nrow(dat0) # = 3*4*3*5*3

# Generate all configurations of 5 dichotomous factors that are compatible with the causal
# chain (A*b + a*B <-> C)*(C*d + c*D <-> E).
dat1 <- allCombs(c(2, 2, 2, 2, 2)) - 1
(dat2 <- selectCases("(A*b + a*B <-> C)*(C*d + c*D <-> E)", dat1))

# Generate all configurations of 5 multi-value factors that are compatible with the causal
# chain (A=2*B=1 + A=3*B=3 <-> C=1)*(C=1*D=2 + C=4*D=4 <-> E=3).
dat1 <- allCombs(c(3, 3, 4, 4, 3))
dat2 <- selectCases("(A=2*B=1 + A=3*B=3 <-> C=1)*(C=1*D=2 + C=4*D=4 <-> E=3)", dat1, 
                    type = "mv")
nrow(dat1)
nrow(dat2)                    

# Generate all configurations of 5 fuzzy-set factors that are compatible with the causal
# structure A*b + C*D <-> E, such that con = .8 and cov = .8.
dat1 <- allCombs(c(2, 2, 2, 2, 2)) - 1
dat2 <- makeFuzzy(dat1, fuzzvalues = seq(0, 0.45, 0.01))
(dat3 <- selectCases1("A*b + C*D <-> E", con = .8, cov = .8, dat2))

# Inverse search for the data generating causal structure A*b + a*B + C*D <-> E from
# fuzzy-set data with non-perfect consistency and coverage scores.
set.seed(3)
groundTruth <- "A*b + a*B + C*D <-> E"
dat1 <- allCombs(c(2, 2, 2, 2, 2)) - 1
dat2 <- makeFuzzy(dat1, fuzzvalues = 0:4/10)
dat3 <- selectCases1(groundTruth, con = .8, cov = .8, dat2)
ana1 <- fscna(dat3, ordering = list("E"), strict = TRUE, con = .8, cov = .8)
any(is.submodel(asf(ana1)$condition, groundTruth))
}
