\name{is.submodel}
\alias{is.submodel}
\alias{identical.model}

\title{
Identify correctness-preserving submodel relations
}
\description{
The function \code{is.submodel} checks for each element of a vector of \code{\link{cna}} solution formulas whether it is a submodel of a specified target model \code{y}. If \code{y} is the true model in an inverse search (i.e. the ground truth), \code{is.submodel} identifies the correct models in the \code{\link{cna}} output (see Baumgartner and Thiem 2017, Baumgartner and Ambuehl 2018).
}

\usage{
is.submodel(x, y, strict = FALSE)
identical.model(x, y)
}

\arguments{
  \item{x}{Character vector of atomic and/or complex solution formulas (asf/csf). Must be of length 1 in \code{identical.model}.}
  \item{y}{Character string of length 1 specifying the target asf or csf.}
  \item{strict}{Logical; if \code{TRUE}, the elements of \code{x} only count as submodels of \code{y} if they are proper parts of \code{y} (i.e. not identical to \code{y}).}
}
\details{

To benchmark the reliability of a method of causal inference it must be tested to what degree the method recovers the true data-generating structure \eqn{\Delta} or proper substructures of \eqn{\Delta} from data of varying quality. Reliability benchmarking is done in so-called \emph{inverse searches}, which reverse the order of causal discovery as normally conducted in scientific practice. An inverse search comprises three steps: (1) a causal structure \eqn{\Delta} is drawn/presupposed (as ground truth), (2) artificial data \eqn{\delta} is simulated from \eqn{\Delta}, possibly featuring various deficiencies (e.g. noise, limited diversity, measurement error etc.), and (3) \eqn{\delta} is processed by the benchmarked method in order to check whether its output meets the tested reliability benchmark (e.g. whether the output is true of or identical to \eqn{\Delta}). 

The main purpose of \code{is.submodel} is to execute step (3) of an inverse search that is tailor-made to test the reliability of \code{\link{cna}} [with \code{\link{randomConds}} and \code{\link{selectCases}} designed for steps (1) and (2), respectively]. A solution formula \code{x} being a submodel of a target formula \code{y} means that all the causal claims entailed by \code{x} are true of \code{y}, which is the case if a causal interpretation of \code{x} entails conjunctive and disjunctive causal relevance relations that are all likewise entailed by a causal interpretation of \code{y}. More specifically, \code{x} is a submodel of \code{y} if, and only if, the following conditions are satisfied: (i) all factor values causally relevant according to \code{x} are also causally relevant according to \code{y}, (ii) all factor values contained in two different disjuncts in \code{x} are also contained in two different disjuncts in \code{y}, (iii) all factor values contained in the same conjunct in \code{x} are also contained in the same conjunct in \code{y}, and (iv) if \code{x} is a csf with more than one asf, (i) to (iii) are satisfied for all asfs in \code{x}. For more details see Baumgartner and Thiem (2020) or Baumgartner and Ambuehl (2018, online appendix).



\code{is.submodel} requires two inputs \code{x} and \code{y}, where \code{x} is a character vector of \code{\link{cna}} solution formulas (asf or csf) and \code{y} is one asf or csf (i.e. a character string of length 1), viz. the target structure or ground truth. The function returns \code{TRUE} for elements of \code{x} that are a submodel of \code{y} according to the definition of submodel-hood given in the previous paragraph. If \code{strict = TRUE}, \code{x} counts as a submodel of \code{y} only if \code{x} is a proper part of \code{y} (i.e. \code{x} is not identical to \code{y}).

The function \code{identical.model} returns \code{TRUE} only if \code{x} (which must be of length 1) and \code{y} are identical. It can be used to test whether \code{y} is completely recovered in an inverse search.
}

\value{Logical vector of the same length as \code{x}.}

\references{
Baumgartner, Michael and Mathias Ambuehl. 2018. \dQuote{Causal Modeling with Multi-Value and Fuzzy-Set Coincidence Analysis.} \emph{Political Science Research and Methods}. doi:10.1017/psrm.2018.45. 

Baumgartner, Michael and Alrik Thiem. 2020. \dQuote{Often Trusted But Never (Properly) Tested: Evaluating Qualitative Comparative Analysis}. \emph{Sociological Methods & Research} 49:279-311.

}

\seealso{
\code{\link{randomConds}}, \code{\link{selectCases}}, \code{\link{cna}}.
}

\examples{
# Binary expressions
# ------------------
trueModel.1 <- "(A*b + a*B <-> C)*(C*d + c*D <-> E)"
candidates.1 <- c("(A + B <-> C)*(C + c*D <-> E)", "(A + B <-> C)", 
                 "(A <->  C)*(C <-> E)", "(C <-> E)")
candidates.2 <- c("(A*B + a*b <-> C)*(C*d + c*D <-> E)", "(A*b*D + a*B <-> C)", 
                 "(A*b + a*B <-> C)*(C*A*D <-> E)", "(D <-> C)", 
                 "(A*b + a*B + E <-> C)*(C*d + c*D <-> E)")

is.submodel(candidates.1, trueModel.1)
is.submodel(candidates.2, trueModel.1)
is.submodel(c(candidates.1, candidates.2), trueModel.1)

is.submodel("C + b*A <-> D", "A*b + C <-> D")
is.submodel("C + b*A <-> D", "A*b + C <-> D", strict = TRUE)
identical.model("C + b*A <-> D", "A*b + C <-> D")

target.1 <- "(A*b + a*B <-> C)*(C*d + c*D <-> E)"
testformula.1 <- "(A*b + a*B <-> C)*(C*d + c*D <-> E)*(A + B <-> C)"
is.submodel(testformula.1, target.1)

# Multi-value expressions
# -----------------------
trueModel.2 <- "(A=1*B=2 + B=3*A=2 <-> C=3)*(C=1 + D=3 <-> E=2)"
is.submodel("(A=1*B=2 + B=3 <-> C=3)*(D=3 <-> E=2)", trueModel.2)
is.submodel("(A=1*B=1 + B=3 <-> C=3)*(D=3 <-> E=2)", trueModel.2)
is.submodel(trueModel.2, trueModel.2)
is.submodel(trueModel.2, trueModel.2, strict = TRUE)

target.2 <- "C=2*D=1*B=3 + A=1 <-> E=5"
testformula.2 <- c("C=2 + D=1 <-> E=5","C=2 + D=1*B=3 <-> E=5","A=1+B=3*D=1*C=2 <-> E=5",
                "C=2 + D=1*B=3 + A=1 <-> E=5","C=2*B=3 + D=1 + B=3 + A=1 <-> E=5")
is.submodel(testformula.2, target.2)
identical.model(testformula.2[3], target.2)
identical.model(testformula.2[1], target.2)
}
