
\name{selectCases}
\alias{selectCases}
\alias{selectCases1}

\title{
Select the cases/configurations compatible with a data generating causal structure
}

\description{
\code{selectCases} selects the cases/configurations that are compatible with a Boolean function, in particular (but not exclusively), a data generating causal structure, from a data frame or \code{configTable}.

\code{selectCases1} allows for setting consistency (\code{con}) and coverage (\code{cov}) thresholds. It then selects cases/configurations that are compatible with the data generating structure to degrees \code{con} and \code{cov}.
}

\usage{
selectCases(cond, x = full.ct(cond), type, cutoff = 0.5, 
            rm.dup.factors = FALSE, rm.const.factors = FALSE)
selectCases1(cond, x = full.ct(cond), type, con = 1, cov = 1, 
             rm.dup.factors = FALSE, rm.const.factors = FALSE)
}

\arguments{
  \item{cond}{Character string specifying the Boolean function for which compatible cases are to be selected.}
  \item{x}{Data frame or \code{configTable}; if not specified, \code{\link{full.ct}(cond)} is used.}
  \item{type}{Character vector specifying the type of \code{x}: \code{"cs"} (crisp-set), \code{"mv"} (multi-value),  or \code{"fs"} (fuzzy-set). Defaults to the type of \code{x}, if \code{x} is a \code{configTable} or to \code{"cs"} otherwise.}
  \item{cutoff}{Cutoff value in case of \code{"fs"} data.}
  \item{rm.dup.factors}{Logical; if \code{TRUE}, all but the first of a set of factors with identical value distributions are eliminated.}
  \item{rm.const.factors}{Logical; if \code{TRUE}, constant factors are eliminated.}
  \item{con, cov}{Numeric scalars between 0 and 1 to set the minimum consistency and coverage thresholds.}
}

\details{In combination with \code{allCombs}, \code{full.ct}, \code{randomConds} and \code{makeFuzzy}, \code{selectCases} is useful for simulating data, which are needed for inverse search trials benchmarking the output of the \code{cna} function. 

\code{selectCases} draws those cases/configurations from a data frame or \code{configTable} \code{x} that are compatible with a data generating causal structure (or any other Boolean or set-theoretic function), which is given to \code{selectCases} as a character string \code{cond}. If the argument \code{x} is not specified, configurations are drawn from \code{full.ct(cond)}. \code{cond} can be a condition of any of the three types of conditions, \emph{boolean}, \emph{atomic} or \emph{complex} (see \code{\link{condition}}). To illustrate, if the data generating structure is "A + B <-> C", then a case featuring A=1, B=0, and C=1 is selected by \code{selectCases}, whereas a case featuring A=1, B=0, and C=0 is not (because according to the data generating structure, A=1 must be associated with C=1, which is violated in the latter case). The type of the data frame is specified by the argument \code{type} taking \code{"cs"} (crisp-set), \code{"mv"} (multi-value), and \code{"fs"} (fuzzy-set) as values.

\code{selectCases1} allows for providing consistency (\code{con}) and coverage (\code{cov}) thresholds, such that some cases that are incompatible with \code{cond} are also drawn, as long as \code{con} and \code{cov} remain satisfied. The solution is identified by an algorithm aiming at finding a subset of maximal size meeting the \code{con} and \code{cov} requirements. In contrast to  \code{selectCases}, \code{selectCases1} only accepts a condition of type \emph{atomic} as its \code{cond} argument, i.e. an atomic solution formula. Data drawn by \code{selectCases1} can only be modeled with consistency = \code{con} and coverage = \code{cov}. 

}

\value{
A \code{configTable}.
}

\seealso{\code{\link{allCombs}}, \code{\link{full.ct}}, \code{\link{randomConds}}, \code{\link{makeFuzzy}}, \code{\link{configTable}}, \code{\link{condition}}, \code{\link{cna}}, \code{\link{d.jobsecurity}}}

\examples{
# Generate all configurations of 5 dichotomous factors that are compatible with the causal
# chain (A*b + a*B <-> C) * (C*d + c*D <-> E).
groundTruth.1 <- "(A*b + a*B <-> C) * (C*d + c*D <-> E)"
(dat1 <- selectCases(groundTruth.1))
condition(groundTruth.1, dat1)

# Randomly draw a multi-value ground truth and generate all configurations compatible with it.
dat1 <- allCombs(c(3, 3, 4, 4, 3))
groundTruth.2 <- randomCsf(mvct(dat1), n.asf=2)
(dat2 <- selectCases(groundTruth.2, dat1, type = "mv"))
condition(groundTruth.2, dat2)

# Generate all configurations of 5 fuzzy-set factors compatible with the causal structure
# A*b + C*D <-> E, such that con = .8 and cov = .8.
dat1 <- allCombs(c(2, 2, 2, 2, 2)) - 1
dat2 <- makeFuzzy(dat1, fuzzvalues = seq(0, 0.45, 0.01))
(dat3 <- selectCases1("A*b + C*D <-> E", con = .8, cov = .8, dat2))
condition("A*b + C*D <-> E", dat3)

# Inverse search for the data generating causal structure A*b + a*B + C*D <-> E from
# fuzzy-set data with non-perfect consistency and coverage scores.
dat1 <- allCombs(c(2, 2, 2, 2, 2)) - 1
set.seed(7)
dat2 <- makeFuzzy(dat1, fuzzvalues = 0:4/10)
dat3 <- selectCases1("A*b + a*B + C*D <-> E", con = .8, cov = .8, dat2)
fscna(dat3, ordering = list("E"), strict = TRUE, con = .8, cov = .8)

# Draw cases satisfying specific conditions from real-life fuzzy-set data.
ct.js <- fsct(d.jobsecurity)
selectCases("S -> C", ct.js)  # Cases with higher membership scores in C than in S.
selectCases("S -> C", d.jobsecurity, type = "fs")  # Same.
selectCases("S <-> C", ct.js) # Cases with identical membership scores in C and in S.
selectCases1("S -> C", con = .8, cov = .8, ct.js)  # selectCases1() makes no distinction 
              #  between "->" and "<->".
condition("S -> C", selectCases1("S -> C", con = .8, cov = .8, ct.js))

# selectCases() not only draws cases compatible with Boolean causal models. Any Boolean or 
# set-theoretic function can be given as cond.
selectCases("C > B", allCombs(2:4), type = "mv")
selectCases("C=2 | B!=3", allCombs(2:4), type = "mv")
selectCases("A=1 * !(C=2 + B!=3)", allCombs(2:4), type = "mv")
}