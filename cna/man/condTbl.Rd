\name{condTbl}
\alias{condTbl}
\alias{print.condTbl}
\alias{as.condTbl}
\alias{condTbl}
\alias{msc}
\alias{asf}
\alias{csf}


\title{
Extract conditions and solutions from an object of class \dQuote{cna}
}

\description{
Given a solution object \code{x} produced by \code{\link{cna}}, \code{msc(x)} extracts all minimally sufficient conditions, \code{asf(x)} all atomic solution formulas, and \code{csf(x, n.init)} builds approximately \code{n.init} complex solution formulas. All solution attributes (\code{details}) that are saved in \code{x} are recovered as well. The three functions return a data frame with the additional class attribute \code{condTbl}.

\code{as.condTbl} reshapes the output produced by \code{\link{condition}} in such a way as to make it identical to the output returned by \code{msc}, \code{asf}, and \code{csf}.

\code{condTbl} executes \code{\link{condition}} and returns a concise summary table featuring consistencies and coverages.
}

\usage{
msc(x, details = x$details)
asf(x, details = x$details, warn_details = TRUE)
csf(x, n.init = 1000, details = x$details,
    asfx = asf(x, details, warn_details = FALSE),
    inus.only = x$inus.only, minimalizeCsf = inus.only,
    acyclic.only = x$acyclic.only, cycle.type = x$cycle.type,
    verbose = FALSE)
\method{print}{condTbl}(x, n = 20, digits = 3, quote = FALSE, row.names = TRUE, ...)

condTbl(...)
as.condTbl(x, ...)
}

\arguments{
  \item{x}{Object of class \dQuote{cna}. In \code{as.condTbl}, \code{x} is a list of evaluated conditions as returned by \code{condition}. In \code{print}, \code{x} is an object returned by \code{msc}, \code{asf}, or \code{csf}.}
  \item{details}{Either \code{TRUE/FALSE} or a character vector specifying which solution attributes to print (see \code{\link{cna}}). Note that \code{msc} and \code{asf} can only display attributes that are saved in \code{x}, i.e. those that have been requested in the \code{details} argument within the call of \code{\link{cna}}.}
  \item{warn_details}{Logical; if \code{TRUE}, a warning is issued when some attribute requested in \code{details} is not available in \code{x} (parameter for internal use).}
  \item{n.init}{Integer capping the amount of initial asf combinations. Default at 1000. Serves to control the computational complexity of the csf building process. %After combining the asf to a preliminary set of csf, the number of these csf is reduced to \code{n.init} in order to limit the computational burden.
  }
  \item{asfx}{Object of class \dQuote{condTbl} produced by the \code{asf} function.}
  \item{inus.only}{Logical; if \code{TRUE}, csf are freed of structural redundancies and only csf not featuring partial structural redundancies are retained  (see the fourth stage of \code{\link{cna}}'s search algorithm). Defaults to \code{x$inus.only}.}
  \item{minimalizeCsf}{Logical; if \code{TRUE}, csf are freed of structural redundancies.}
  \item{acyclic.only}{Logical; if \code{TRUE}, csf featuring a cyclic substructure are not returned. \code{FALSE} by default.}
  \item{cycle.type}{Character string specifying what type of cycles to be detected: \code{"factor"} (the default) or \code{"value"} (cf. \code{\link{cyclic}}).}
  \item{verbose}{Logical; if \code{TRUE}, some details on the csf building process are printed. \code{FALSE} by default.}
  \item{n}{Maximal number of msc, asf, or csf to be printed.}
  \item{digits}{Number of digits to print in consistency, coverage, exhaustiveness, faithfulness, and coherence scores.}
  \item{quote, row.names}{As in \code{\link{print.data.frame}}.}
  \item{\dots}{All arguments in \code{condTbl} are passed on to \code{\link{condition}}.}
}

\details{Depending on the processed data, the solutions (models) output by \code{\link{cna}} are often ambiguous, to the effect that many atomic and complex solutions fit the data equally well. To facilitate the inspection of the \code{\link{cna}} output, however, the latter standardly returns only 5 minimally sufficient conditions (msc), 5 atomic solution formulas (asf), and 5 complex solution formulas (csf) for each outcome. \code{msc} can be used to extract \emph{all} msc from an object \code{x} of class \dQuote{cna}, \code{asf} to extract \emph{all} asf, and \code{csf} to build approximately \code{n.init} csf from the asf stored in \code{x}. All solution attributes (\code{details}) that are saved in \code{x} are recovered as well.
The outputs of \code{msc}, \code{asf}, and \code{csf} can be further processed by the \code{\link{condition}} function.

While \code{msc} and \code{asf} merely extract information stored in \code{x}, \code{csf} builds csf from the inventory of asf recovered in the second phase of the \code{\link{cna}} algorithm. That is, the \code{csf} function implements the fourth stage of that algorithm. It proceeds in a stepwise manner as follows.
\enumerate{\item \code{n.init} possible conjunctions featuring one asf of every outcome are built. \item  If \code{inus.only = TRUE} or \code{minimalizeCsf = TRUE}, the solutions resulting from step 1 are freed of structural redundancies (cf. Baumgartner and Falk 2019). 
\item If \code{inus.only = TRUE}, tautologous and contradictory solutions as well as solutions with partial structural redundancies and constant factors are eliminated. 
[If \code{inus.only = FALSE} and \code{minimalizeCsf = TRUE}, only structural redundancies are eliminated, meaning only step 2, but not step 3, is executed.]
\item If \code{acyclic.only = TRUE}, solutions with cyclic substructures are eliminated. \item For those solutions that were modified in the previous steps, consistency and coverage are re-calculated and solutions that no longer reach \code{con} or \code{cov} are eliminated. 
\item The remaining solutions are returned as csf, ordered by complexity and the product of consistency and coverage.}


The argument \code{digits} applies to the \code{print} method. It determines how many digits of consistency, coverage, exhaustiveness, faithfulness, and coherence scores are printed. The default value is 3.

The function \code{as.condTbl} takes a list of objects of class \dQuote{cond} that are returned by the \code{\link{condition}} function as input, and reshapes these objects in such a way as to make them identical to the output returned by \code{msc}, \code{asf}, and \code{csf}.

\code{condTbl(...)} is identical with \code{as.condTbl(condition(...))}.
}

\value{
\code{msc}, \code{asf}, \code{csf}, and \code{as.condTbl} return objects of class \dQuote{condTbl}, a \code{data.frame} which features the following components:
\tabular{rl}{
\code{outcome}: \tab the outcomes\cr
\code{condition}:\tab the relevant conditions or solutions\cr
\code{consistency}:\tab the consistency scores \cr
\code{coverage}: \tab the coverage scores  \cr
\code{complexity}: \tab the complexity scores\cr
\code{inus}: \tab whether the solutions have INUS form\cr
\code{exhaustiveness}: \tab the exhaustiveness scores\cr
\code{faithfulness}: \tab the faithfulness scores\cr
\code{coherence}: \tab the coherence scores\cr
\code{redundant}: \tab whether the csf contain redundant proper parts\cr
\code{cyclic}: \tab whether the csf contain cyclic substructures
}
The latter five measures are optional and will be appended to the table according to the setting of the argument \code{details}.
}


\section{Contributors}{
Falk, Christoph: development, testing
}

\references{
Baumgartner, Michael and Christoph Falk. 2019. \dQuote{Boolean Difference-Making: A Modern Regularity Theory of Causation}. \emph{The British Journal for the Philosophy of Science}. doi:10.1093/bjps/axz047.

Lam, Wai Fung, and Elinor Ostrom. 2010.
\dQuote{Analyzing the Dynamic Complexity of Development Interventions: Lessons
from an Irrigation Experiment in Nepal.}
\emph{Policy Sciences} 43 (2):1-25.
}


\seealso{\code{\link{cna}}, \code{\link{configTable}}, \code{\link{condition}}, \code{\link{minimalizeCsf}}, \code{\link{d.irrigate}}}

\examples{
# Crisp-set data from Lam and Ostrom (2010) on the impact of development interventions
# ------------------------------------------------------------------------------------
# CNA with causal ordering that corresponds to the ordering in Lam & Ostrom (2010); coverage 
# cut-off at 0.9 (consistency cut-off at 1).
cna.irrigate <- cna(d.irrigate, ordering = list(c("A","R","F","L","C"),"W"), cov = .9,
                    maxstep = c(4, 4, 12), details = TRUE)
cna.irrigate

# The previous function call yields a total of 12 complex solution formulas, only
# 5 of which are returned in the default output. 
# Here is how to extract all 12 complex solution formulas along with all 
# solution attributes.
csf(cna.irrigate)
# With only the standard attributes plus exhaustiveness and faithfulness.
csf(cna.irrigate, details = c("e", "f"))

# Extract all atomic solution formulas.
asf(cna.irrigate)

# Extract all minimally sufficient conditions.
msc(cna.irrigate)

# Extract only the conditions (solutions).
csf(cna.irrigate)$condition
asf(cna.irrigate)$condition
msc(cna.irrigate)$condition

# A CNA of d.irrigate without a presupposed ordering is even more ambiguous.
cna2.irrigate <- cna(d.irrigate, cov = .9, maxstep = c(4,4,12), details = TRUE)

# To speed up the construction of complex solution formulas, first extract atomic solutions
# and then pass these asf to csf.
cna2.irrigate.asf <- asf(cna2.irrigate)
\donttest{csf(cna2.irrigate, asfx = cna2.irrigate.asf, details = FALSE)}
# Reduce the initial asf combinations.
csf(cna2.irrigate, asfx = cna2.irrigate.asf, n.init = 50)
# Print the first 20 csf.
csf(cna2.irrigate, asfx = cna2.irrigate.asf, n.init = 50)[1:20, ]
# Also extract exhaustiveness scores.
csf(cna2.irrigate, asfx = cna2.irrigate.asf, n.init = 50,
    details = "e")[1:20, ]

# Print details about the csf building process.
csf(cna.irrigate, verbose = TRUE)

# Return solution attributes with 5 digits.
print(cna2.irrigate.asf, digits = 5)

# Further examples
# ----------------
# An example generating structural redundancies.
target <- "(A*B + C <-> D)*(c + a <-> E)"
dat1 <- selectCases(target)
ana1 <- cna(dat1, maxstep = c(3, 4, 10))
# Run csf with elimination of structural redundancies.
csf(ana1, verbose = TRUE)
# Run csf without elimination of structural redundancies. 
csf(ana1, verbose = TRUE, inus.only = FALSE)

# An example generating partial structural redundancies.
dat2 <- data.frame(A = c(1, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1),
                   B = c(1, 1, 0, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 1, 0, 1),
                   C = c(1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 0),
                   D = c(1, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                   E = c(0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0),
                   F = c(1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0))
ana2 <- cna(dat2, con = .8, cov = .8, maxstep = c(3, 4, 10))
# Run csf without elimination of partial structural redundancies.
csf(ana2, inus.only = FALSE, verbose = TRUE)
# Run csf with elimination of partial structural redundancies.
csf(ana2, verbose = TRUE)
# Only acyclic models.
csf(ana2, verbose = TRUE, acyclic.only = TRUE)

# Feed the outputs of msc, asf, and csf into the condition function to further inspect the 
# properties of minimally sufficient conditions and atomic and complex solution formulas.
condition(msc(ana2)$condition, dat2)
condition(asf(ana2)$condition, dat2)
condition(csf(ana2)$condition, dat2)

# Reshape the output of the condition function in such a way as to make it identical to the
# output returned by msc, asf, and csf.
as.condTbl(condition(msc(ana2)$condition, dat2))
as.condTbl(condition(asf(ana2)$condition, dat2))
as.condTbl(condition(csf(ana2)$condition, dat2))

condTbl(csf(ana2)$condition, dat2) # Same as preceding line
}
