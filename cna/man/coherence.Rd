\name{coherence}
\alias{coherence}

\title{
Calculate the coherence of complex solution formulas
}

\description{
Calculates the coherence measure of complex solution formulas (csf).
}

\usage{
coherence(cond, x, type, tt)
}

\arguments{
  \item{cond}{Character vector specifying an asf or csf.}
  \item{x}{Data frame or \code{configTable}.}
  \item{type}{Character vector specifying the type of \code{x}: \code{"cs"} (crisp-set), \code{"mv"} (multi-value),  or \code{"fs"} (fuzzy-set). Defaults to the type of \code{x}, if \code{x} is a \code{configTable} or to \code{"cs"} otherwise.}
  \item{tt}{Argument \code{tt} is deprecated in \code{coherence()}; use \code{x} instead.}
}

\details{
Coherence is a measure for model fit that is custom-built for complex solution formulas (csf). It measures the degree to which the atomic solution formulas (asf) combined in a csf cohere, i.e. are instantiated together in \code{x} rather than independently of one another. More concretely, coherence is the ratio of the number of cases satisfying all asf contained in a csf to the number of cases satisfying at least one asf in the csf. For example, if the csf contains the three asf asf1, asf2, asf3, coherence amounts to | asf1 * asf2 * asf3 | / | asf1 + asf2 + asf3 |, where |...| expresses the cardinality of the set of cases in \code{x} instantiating the corresponding expression. For asf, \code{coherence} returns 1. For boolean conditions (see \code{\link{condition}}), the coherence measure is not defined and \code{coherence} hence returns \code{NA}. For multiple csf that do not have a factor in common, \code{coherence} returns the minimum of the separate coherence scores. 
}

\value{
Numeric vector of coherence values to which \code{cond} is appended as a "names" attribute. If \code{cond} is a \emph{csf} "asf1*asf2*asf3" composed of asf that do not have a factor in common, the \emph{csf} is rendered with commas in the "names" attribute: "asf1, asf2, asf3".
}

\seealso{\code{\link{cna}}, \code{\link{condition}}, \code{\link{selectCases}}, \code{\link{allCombs}}, \code{\link{full.ct}}, \code{\link{condTbl}}}

\examples{
# Perfect coherence.
dat1 <- selectCases("(A*b <-> C)*(C + D <-> E)")
coherence("(A*b <-> C)*(C + D <-> E)", dat1)
csf(cna(dat1, details = "coherence"))

# Non-perfect coherence.
dat2 <- selectCases("(a*B <-> C)*(C + D <-> E)*(F*g <-> H)")
dat3 <- rbind(ct2df(dat2), c(0,1,0,1,1,1,0,1))
coherence("(a*B <-> C)*(C + D <-> E)*(F*g <-> H)", dat3)
csf(cna(dat3, con = .88, details = "coherence")) 
}

