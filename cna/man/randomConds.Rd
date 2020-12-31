\name{randomConds}
\alias{randomConds}
\alias{randomAsf}
\alias{randomCsf}

\title{Generate random solution formulas}

\description{
Based on a set of factors and a corresponding data type---given as a data frame or \code{configTable}---, \code{randomAsf} generates a random atomic solution formula (asf) and \code{randomCsf} a random (acyclic) complex solution formula (csf).
}

\usage{
randomAsf(x, outcome = NULL, compl = NULL, how = c("inus", "minimal"))
randomCsf(x, outcome = NULL, n.asf = NULL, compl = NULL)
}

\arguments{
  \item{x}{Data frame or \code{configTable}; determines the number of factors, their names and their possible values.}
  \item{outcome}{Optional character vector (of length 1 in \code{randomAsf}) specifying the outcome factor(s) in the solution formula; must be a subset of \code{names(x)}.}
  \item{compl}{Integer vector specifying the maximal complexity of the formula (i.e. number of factors in msc; number of msc in asf).}
  \item{how}{Character string, either \code{"inus"} or \code{"minimal"}, specifying whether the generated solution formula is redundancy-free relative to \code{full.ct(x)} or relative to \code{x} (see details).}
  \item{n.asf}{Integer scalar specifying the number of asf in the csf. Is overridden by \cr\code{length(outcome)} if \code{outcome} is not \code{NULL}. Note that \code{n.asf} is limited to \code{ncol(x)-2}.}
}


\details{
\code{randomAsf} and \code{randomCsf} can be used to randomly draw data generating structures (ground truths) in inverse search trials benchmarking the output of \code{\link{cna}}. In the regularity theoretic context in which the CNA method is embedded, a causal structure is a redundancy-free Boolean dependency structure. Hence, \code{randomAsf} and \code{randomCsf} both produce redundancy-free Boolean dependency structures. \code{randomAsf} generates structures with one outcome, i.e. atomic solution formulas (asf), \code{randomCsf} generates structures with multiple outcomes, i.e. complex solution formulas (csf), that are free of cyclic substructures. In a nutshell, \code{randomAsf} proceeds by, first, randomly drawing disjunctive normal forms (DNFs) and by, second, eliminating redundancies from these DNFs. \code{randomCsf} essentially consists in repeated executions of \code{randomAsf}.

The only mandatory argument of \code{randomAsf} and \code{randomCsf} is a data frame or a \code{configTable} \code{x} defining the factors (with their possible values) from which the generated asf and csf shall be drawn. If asf and csf are built from multi-value or fuzzy-set factors, \code{x} must be a \code{configTable}.

The optional argument \code{outcome} determines which factors in \code{x} shall be treated as outcomes. If \code{outcome} is at its default value \code{NULL}, \code{randomAsf} and \code{randomCsf} randomly draw factor(s) from \code{x} to be treated as outcome(s). 

The argument \code{compl} controls the complexity of the generated asf and csf. More specifically, the \emph{initial} complexity of asf and csf (i.e. the number of factors included in msc and the number of msc included in asf prior to redundancy elimination) is drawn from the vector \code{compl}. As this complexity might be reduced in the subsequent process of redundancy elimination, issued asf or csf will often have lower complexity than specified in \code{compl}. The default value of \code{compl} is determined by the number of columns in \code{x}. Assigning unduly high values to \code{compl} results in an error.


\code{randomAsf} has the additional argument \code{how} with the two possible values \code{"inus"} and \code{"minimal"}. \code{how = "inus"} determines that the generated asf is redundancy-free relative to all logically possible configurations of the factors in \code{x}, i.e. relative to \code{full.ct(x)}, whereas in case of \code{how = "minimal"} redundancy-freeness is imposed only relative to all configurations actually contained in \code{x}, i.e. relative to \code{x} itself. Typically \code{"inus"} should be used; the value \code{"minimal"} is relevant mainly in repeated \code{randomAsf} calls from within \code{randomCsf}. Moreover, setting \code{how = "minimal"} will return an error if \code{x} is a \code{\link{configTable}} of type \code{"fs"}.

The argument \code{n.asf} controls the number of asf in the generated csf. Its value is limited to \code{ncol(x)-2} and overridden by \code{length(outcome)}, if \code{outcome} is not \code{NULL}. Analogously to \code{compl}, \code{n.asf} specifies the number of asf prior to redundancy elimination, which, in turn, may further reduce these numbers. That is, \code{n.asf} provides an upper bound for the number of asf in the resulting csf.
}

\value{The randomly generated formula, a character string.}

\seealso{
\code{\link{is.submodel}}, \code{\link{selectCases}}, \code{\link{full.ct}}, \code{\link{configTable}}, \code{\link{cna}}.
}

\examples{# randomAsf
# ---------
# Asf generated from explicitly specified binary factors.
randomAsf(full.ct("H*I*T*R*K"))
randomAsf(full.ct("Johnny*Debby*Aurora*Mars*James*Sonja"))

# Asf generated from a specified number of binary factors.
randomAsf(full.ct(7))

# Asf generated from an existing data frame.
randomAsf(d.educate)

# Specify the outcome.
randomAsf(d.educate, outcome = "G")

# Specify the complexity.
randomAsf(full.ct(7), compl = 2)
randomAsf(full.ct(7), compl = 3:4)

# Redundancy-freeness relative to x instead of full.ct(x).
randomAsf(d.educate, outcome = "G", how = "minimal")

# Asf with multi-value factors (x must be given as a configTable).
randomAsf(mvct(allCombs(c(3,4,3,5,3,4))))

# Asf from fuzzy-set data (x must be given as a configTable).
randomAsf(fsct(d.jobsecurity))
randomAsf(fsct(d.jobsecurity), outcome = "JSR")

# Generate 20 asf.
\donttest{replicate(20, randomAsf(full.ct(7), compl = 2:3))}


# randomCsf
# ---------
# Csf generated from explicitly specified binary factors.
randomCsf(full.ct("H*I*T*R*K*Q*P"))

# Csf generated from a specified number of binary factors.
randomCsf(full.ct(7))

# Specify the outcomes.
randomCsf(d.volatile, outcome = c("RB","SE"))

# Specify the complexity.
randomCsf(d.volatile, outcome = c("RB","SE"), compl = 2)
randomCsf(full.ct(7), compl = 3:4)

# Specify the number of asf.
randomCsf(full.ct(7), n.asf = 3)

# Csf with multi-value factors (x must be given as a configTable).
randomCsf(mvct(allCombs(c(3,4,3,5,3,4))))

# Generate 20 csf.
\donttest{replicate(20, randomCsf(full.ct(7), n.asf = 2, compl = 2:3))


# Inverse searches
# ----------------
# === Ideal Data ===
# Draw the data generating structure. (Every run yields different 
# targets and data.)
target <- randomCsf(full.ct(5), n.asf = 2)
target
# Select the cases compatible with the target.
x <- selectCases(target)
# Run CNA without an ordering.
mycna <- cna(x, rm.dup.factors = FALSE)
# Extract the csf. 
csfs <- csf(mycna)
# Check whether the target is completely returned.
any(unlist(lapply(csfs$condition, identical.model, target)))

# === Data fragmentation (20% missing observations) ===
# Draw the data generating structure. (Every run yields different 
# targets and data.)
target <- randomCsf(full.ct(7), n.asf = 2)
target
# Generate the ideal data.
x <- ct2df(selectCases(target))
# Introduce fragmentation.
x <- x[-sample(1:nrow(x), nrow(x)*0.2), ] 
# Run CNA without an ordering.
mycna <- cna(x, rm.dup.factors = FALSE)
# Extract the csf. 
csfs <- csf(mycna)
# Check whether (a submodel of) the target is returned.
any(is.submodel(csfs$condition, target))

# === Data fragmentation and noise (20% missing observations, noise ratio of 0.05) ===
# Multi-value data.
# Draw the data generating structure. (Every run yields different 
# targets and data.)
fullData <- mvct(allCombs(c(4,4,4,4,4)))
target <- randomCsf(fullData, n.asf=2, compl = 2:3)
target
# Generate the ideal data.
x <- ct2df(selectCases(target, fullData))
# Introduce fragmentation.
x <- x[-sample(1:nrow(x), nrow(x)*0.2), ] 
# Introduce random noise.
x <- rbind(ct2df(fullData[sample(1:nrow(fullData), nrow(x)*0.05), ]), x)  
# Run CNA without an ordering.
mycna <- mvcna(x, con = .75, cov = .75, maxstep = c(3, 3, 12), rm.dup.factors = F)
# Extract the csf. 
csfs <- csf(mycna)
# Check whether no causal fallacy (no false positive) is returned.
if(nrow(csfs)==0) {
      TRUE } else {any(is.submodel(csfs$condition, target))}
}
}