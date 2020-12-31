\name{minimalize}
\alias{minimalize}

\title{
Eliminate logical redundancies from Boolean expressions
}

\description{
\code{minimalize} eliminates logical redundancies from a Boolean expression \code{cond} based on all configurations of the factors in \code{cond} that are possible according to classical Boolean logic. That is, \code{minimalize} performs logical (i.e. not data-driven) redundancy elimination. The output is a set of redundancy-free DNFs that are logically equivalent to \code{cond}.
}

\usage{
minimalize(cond, x = NULL, maxstep = c(4, 4, 12))
}

\arguments{
  \item{cond}{
Character vector specifying Boolean expressions; the acceptable syntax is the same as that of \code{\link{condition}}.
}
  \item{x}{A data frame, a \code{configTable}, or a list determining the possible values for each factor in \code{cond}; \code{x} has no effect for a \code{cond} with only binary factors but is mandatory for a \code{cond} with multi-value factors (see details).
}
  \item{maxstep}{
Maximal complexity of the returned redundancy-free DNFs (see \code{\link{cna}}).
}
}

\details{ 
The regularity theory of causation underlying CNA conceives of causes as parts of redundancy-free Boolean dependency structures. Boolean dependency structures tend to contain a host of redundancies. %Some of these redundancies are relative to the empirical data from which the dependency structures are inferred, others are of principled, logical nature. Logical redundancies are not relative to data; they obtain in all possible data sets. A Boolean expression with logical redundancies, therefore, cannot represent a causal structure for principled reasons. 
Redundancies may obtain relative to an analyzed set of empirical data, which, typically, are fragmented and do not feature all logically possible configurations, or they may obtain for principled logical reasons, that is, relative to all configurations that are possible according to Boolean logic.
Whether a Boolean expression (in disjunctive normal form) contains the latter type of logical redundancies can be checked with the function \code{\link{is.inus}}. 

\code{minimalize} eliminates logical redundancies from \code{cond} and outputs all redundancy-free disjunctive normal forms (DNF) (within some complexity range given by \code{maxstep}) that are logically equivalent with \code{cond}.
If \code{cond} is redundancy-free, no reduction is possible and \code{minimalize} returns \code{cond} itself (possibly as an element of multiple logically equivalent redundancy-free DNFs). If \code{cond} is not redundancy-free, a \code{\link{cna}} with \code{con = 1} and \code{cov = 1} is performed relative to \code{\link{full.ct}(x)} (relative to \code{full.ct(cond)} if \code{x} is \code{NULL}). The output is the set of all redundancy-free DNFs in the complexity range given by \code{maxstep} that are logically equivalent to \code{cond}.

The purpose of the optional argument \code{x} is to determine the space of possible values of the factors in \code{cond}. If all factors in \code{cond} are binary, \code{x} is optional and without influence on the output of \code{minimalize}. If some factors in \code{cond} are multi-value, \code{minimalize} needs to be given the range of these values. \code{x} can be a data frame or \code{configTable} listing all possible value configurations or simply a list of the possible values for each factor in \code{cond} (see examples).

The argument \code{maxstep}, which is identical to the corresponding argument in \code{\link{cna}}, specifies the maximal complexity of the returned DNF. \code{maxstep} expects a vector of three integers \code{c(i, j, k)}  determining that the generated DNFs have maximally \code{j} disjuncts with maximally \code{i} conjuncts each and a total of maximally \code{k} factor values. The default is \code{maxstep = c(4, 4, 12)}.
If the complexity range of the search space given by \code{maxstep} is too low, it may happen that nothing is returned (accompanied by a corresponding warning message). In that case, the \code{maxstep} values need to be increased. 
}

\value{
A list of character vectors of the same length as \code{cond}. Each list element contains one or several redundancy-free disjunctive normal forms (DNFs) that are logically equivalent to \code{cond}.
}

\seealso{
\code{\link{condition}}, \code{\link{is.inus}}, \code{\link{cna}}, \code{\link{full.ct}}.
}

\examples{
# Binary expressions
# ------------------
# DNFs as input.
minimalize(c("A", "A+B", "A + a*B", "A + a", "A*a"))
minimalize(c("F + f*G", "F*G + f*H + G*H", "F*G + f*g + H*F + H*G"))

# Any Boolean expressions (with variable syntax) are admissible inputs.
minimalize(c("!(A*B*C + a*b*c)", "A*!(B*d+E)->F", "-(A+-(E*F))<->H"))

# Proper redundancy elimination may require increasing the maxstep values.
minimalize("!(A*B*C*D*E+a*b*c*d*e)")  
minimalize("!(A*B*C*D*E+a*b*c*d*e)", maxstep = c(3, 5, 15))


# Multi-value expressions
# -----------------------
# In case of expressions with multi-value factors, the relevant range of factor 
# values must be specified by means of x. x can be a list or a configTable:
values <- list(C = 0:3, F = 0:2, V = 0:4)
minimalize(c("C=1 + F=2*V=0", "C=1 + C=0*V=1"), values)
\donttest{minimalize("C=1 + F=2 <-> V=1", values, maxstep=c(3,40,40))}
minimalize(c("C=1 + C=0 * C=2", "C=0 + C=1 + C=2"), mvct(d.pban))


# Eliminating logical redundancies from non-INUS asf inferred from real data
# --------------------------------------------------------------------------
\donttest{fsdata <- fsct(d.jobsecurity)
conds <- asf(cna(fsdata, con = 0.8, cov = 0.8, inus.only = FALSE))$condition
conds <- cna:::lhs(conds)
noninus.conds <- conds[-which(is.inus(conds, fsdata))]
minimalize(noninus.conds)}
}
