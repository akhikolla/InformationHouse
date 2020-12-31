\name{is.inus}
\alias{is.inus}

\title{
Check whether expressions in the syntax of CNA solutions have INUS form
}
\description{
\code{is.inus} checks for each element of a character vector of disjunctive normal forms (DNFs) or expressions in the syntax of CNA solution formulas whether it has INUS form, meaning whether it is free of redundancies in necessary or sufficient conditions, structural redundancies, partial structural redundancies, whether it has constant factors or identical outcomes, and whether it is tautologous or contradictory.


}

\usage{
is.inus(cond, x = NULL, csf.info = FALSE)
}

\arguments{
  \item{cond}{
Character vector of DNFs or expressions in the syntax of CNA solutions (i.e. asf or csf). }
  \item{x}{An optional argument providing a \code{configTable}, a data frame, or a list specifying the factors' value ranges if \code{cond} contains multi-value factors; if \code{x} is not \code{NULL}, \code{is.inus} tests whether \code{cond} has INUS form relative to \code{full.ct(x)}, otherwise relative to \code{full.ct(cond)}.
}
  \item{csf.info}{Logical; if \code{TRUE} and \code{cond} has the syntax of a csf, details about the performed INUS checks are printed. If \code{cond} does not have the syntax of a csf, \code{csf.info} has no effect.
}
}

\details{
A Boolean dependency structure is not interpretable in terms of a deterministic causal structure if it contains at least one of the following (cf. the \dQuote{Examples} section for examples): \enumerate{ \item redundancies in necessary or sufficient conditions, \item structural redundancies,\item partial structural redundancies, \item constant factors, \item  tautologous or contradictory substructures, \item multiple instances of the same outcome. }



The function \code{is.inus} takes a character vector \code{cond} specifying Boolean disjunctive normal forms (DNFs) or expressions in the syntax of CNA solution formulas as input and runs a series of checks on \code{cond}; one for each of the conditions (1) to (6). For instance, whenever a syntactic proper part of \code{cond} is logically equivalent to \code{cond} itself, the surplus in the latter is redundant, meaning that \code{cond} violates condition (1) and is not causally interpretable. To illustrate, \dQuote{A + a*B <-> C} is logically equivalent to \dQuote{A + B <-> C} and, hence, "a" redundant in the first expression, which is not causally interpretable due to a violation of condition (1).
Or the first asf in \dQuote{(a + C <-> D)*(D + G <-> A)} entails that whenever "a" is given, so is "D", while the second asf entails that whenever "D" is given, so is "A". It follows that "a" cannot ever be given, meaning that the factor A takes the constant value 1 and, hence, violates condition (4). As constant factors can neither be causes nor effects, \dQuote{(a + C <-> D)*(D + G <-> A)} is not a well-formed causal structure. 

If an expression passes the \code{is.inus}-check it can be interpreted as a causal structure according to Mackie's (1974) INUS-theory of causation or modern variants thereof (e.g. Grasshoff and May 2001; Baumgartner and Falk 2019). In other words, such an expression has the form of an INUS structure, i.e. it has \emph{INUS form}, for short.

In the function's default call with \code{x = NULL}, the INUS checks are performed relative to \code{full.ct(cond)}; if \code{x} is not \code{NULL}, the checks are performed relative to \code{full.ct(x)}. As \code{full.ct(cond)} and \code{full.ct(x)} coincide in case of binary factors, the argument \code{x} has no effect in the crisp-set and fuzzy-set cases and, hence, does not have to be specified. In case of multi-value factors, however, the argument \code{x} should be specified in order to define the factors' value ranges (see examples below).

If the argument \code{csf.info} is set to its non-default value \code{TRUE} and \code{cond} has the syntax of a csf, the results of the individual checks of conditions (1) to (6) are printed (in that order) to the console.

In its default setting, the \code{cna} function does not output solutions that do not have INUS form. But when \code{cna} is called with \code{inus.only = FALSE}, non-INUS solutions may be returned. The function \code{is.inus} is standardly called from within the \code{cna} function to determine whether its output has INUS form. 

\code{is.inus} also serves an important purpose in the context of benchmark tests. Not any Boolean expression can be interpreted to represent a causal structure; only expressions in INUS form can. That means when simulating data on randomly drawn target structures, it must be ensured that the latter have INUS form. An expression as \dQuote{A + a*B <-> C}, which has a logically equivalent proper part and, hence, does not have INUS form, is not a well-formed causal structure that could be used as a search target in a benchmark test.




}

\value{
Logical vector of the same length as \code{cond}; if \code{cond} is a csf and \code{is.inus} is called with \code{csf.info = TRUE}, an attribute \dQuote{csf.info} is added.
}

\references{
Baumgartner, Michael and Christoph Falk. 2019. \dQuote{Boolean Difference-Making: A Modern Regularity Theory of Causation}. \emph{The British Journal for the Philosophy of Science}. doi:10.1093/bjps/axz047.

Grasshoff, Gerd and Michael May. 2001. \dQuote{Causal Regularities.} In W Spohn, M Ledwig, M Esfeld (eds.), \emph{Current Issues in Causation}, pp. 85-114. Mentis, Paderborn.

Mackie, John L. 1974. \emph{The Cement of the Universe: A Study of Causation}. Oxford: Oxford University Press.
}

\seealso{
\code{\link{condition}}, \code{\link{full.ct}}, \code{\link{redundant}}, \code{\link{minimalize}}, \code{\link{cna}}, \code{\link{minimalizeCsf}}
}


\examples{
# Crisp-set case
# --------------
# Testing disjunctive normal forms. 
is.inus(c("A", "A + B", "A + a*B", "A + a", "A*a", "A*a + B"))
  
# Testing expressions in the syntax of atomic solution formulas.
is.inus(c("A + B <-> C", "A + B <-> c", "A + a*B <-> C", "A*a + B <-> C", "A + a <-> C", 
          "F*G + f*g + H <-> E", "F*G + f*g + H*F + H*G <-> E")) 

# Testing expressions in the syntax of complex solution formulas.
is.inus(c("(A + B <-> C)*(c + E <-> D)", "(A <-> B)*(B <-> C)", "(A <-> B)*(B <-> C)*(C <-> D)", 
          "(A <-> B)*(B <-> a)", "(A*B + c <-> D)*(E + f <-> D)",
          "(A + B <-> C)*(B*c + E <-> D)"))

# A redundancy in necessary or sufficient conditions, i.e. 
# a non-INUS asf in a csf.
is.inus("(A + A*B <-> C)*(B + D <-> E)", csf.info = TRUE)

# A structural redundancy in a csf.
cond1 <- "(e + a*D <-> C)*(C + A*B <-> D)*(a + c <-> E)"
is.inus("(e + a*D <-> C)*(C + A*B <-> D)*(a + c <-> E)", csf.info = TRUE)
# The first asf in cond1 is redundant.
minimalizeCsf(cond1, selectCases(cond1))

# A partial structural redundancy in a csf.
cond2 <- "(A + B*c + c*E <-> D)*(B + C <-> E)"
is.inus(cond2, csf.info = TRUE)
# The second or third disjunct in the first asf of cond2 is redundant.
cna(selectCases(cond2))

# A csf entailing that one factor is constant.
is.inus("(a + C <-> D)*(D + G <-> A)", csf.info = TRUE)

# A contradictory (i.e. logically constant) csf.
is.inus("(A <-> B)*(B <-> a)", csf.info = TRUE)

# A csf with multiple identical outcomes.
is.inus("(A + C <-> B)*(C + E <-> B)", csf.info = TRUE)


# Multi-value case
# ----------------
# In case of multi-value data, is.inus() needs to be given a data set x determining
# the value ranges of the factors in cond.
mvdata <- mvct(setNames(allCombs(c(2, 3, 2, 3)) -1, c("C", "F", "V", "O")))
is.inus("C=1 + F=2*V=0 <-> O=2", mvdata)
# x can also be given to is.inus() as a list.
is.inus("C=1 + F=2*V=0 <-> O=2", list(C=0:1, F=0:2, V=0:1, O=0:2))
# When x is NULL, is.inus() is applied to full.ct("C=1 + F=2*V=0"), which has only
# one single row. That row is then interpreted to be the only possible configuration, 
# in which case C=1 + F=2*V=0 is tautologous and, hence, non-INUS.
is.inus("C=1 + F=2*V=0 <-> O=2") 
        
is.inus("C=1 + C=0*C=2", mvct(d.pban))    # contradictory
is.inus("C=0 + C=1 + C=2", mvct(d.pban))  # tautologous

# A redundancy in necessary or sufficient conditions, i.e. a 
# non-INUS asf in a csf.
fullDat <- full.ct(list(A=1:3, B=1:3, C=1:3, D=1:3, E=1:3))
is.inus("(A=1 + A=1*B=2 <-> C=3)*(B=2 + D=3 <-> E=1)", fullDat, csf.info = TRUE)

# A structural redundancy in a csf.
cond3 <- "(E=2 + C=1*D=3 <-> A=1)*(A=3*E=1 + C=2*D=2 <-> B=3)*(A=1*E=3 + D=2*E=3 <-> C=1)*
          (A=1*C=2 + A=1*C=3 <-> E=2)"
is.inus(cond3, fullDat, csf.info = TRUE)
# The last asf in cond3 is redundant.
minimalizeCsf(cond3, selectCases(cond3, fullDat))

# A partial structural redundancy in a csf.
cond4 <- "(B=2*C=3 + C=2*D=1 + B=2*C=1*D=2*E=1 <-> A=2)*(D=2*E=1 + D=3*E=1 <-> B=1)"
is.inus(cond4, fullDat, csf.info = TRUE)
# The third disjunct in the first asf of cond4 is redundant.
cna(selectCases(cond4, fullDat))

# A csf entailing that one factor is constant. (I.e. D is constantly ~(D=1).)
cond5 <- "(A=1 + B=2 + E=3 <->C=3)*(A=1*C=1 + B=2*C=1 <-> D=1)"
is.inus(cond5, fullDat, csf.info = TRUE)

# A contradictory csf.
is.inus("(A=1 <-> C=1)*(A=1 <-> C=2)*(A=1 <-> C=3)", fullDat, csf.info = TRUE)

# A csf with multiple identical outcomes.
is.inus("(A=1 + B=2 + D=3 <-> C=1)*(A=2 + B=3 + D=2 <-> C=1)", fullDat, csf.info = TRUE)


# Fuzzy-set case 
# --------------
fsdata <- fsct(d.jobsecurity)
conds <- csf(cna(fsdata, con = 0.85, cov = 0.85, inus.only = FALSE))$condition
# Various examples of different types.
is.inus(conds, fsdata, csf.info = TRUE)
is.inus(c("S + s", "S + s*R", "S*s"), fsdata)

# A redundancy in necessary or sufficient conditions, i.e. a 
# non-INUS asf in a csf.
is.inus("(S + s*L <-> JSR)*(R + P <-> V)", fsdata, csf.info = TRUE)

# A structural redundancy in a csf.
is.inus("(s + l*R <-> C)*(C + L*V <-> R)*(l + c <-> S)", fsdata, csf.info = TRUE)

# A partial structural redundancy in a csf.
is.inus("(S + L*c + c*R <-> P)*(L + C <-> R)", fsdata, csf.info = TRUE)

# A csf entailing that one factor is constant.
is.inus("(S + L <-> P)*(L*p <-> JSR)", csf.info = TRUE)

# A contradictory csf.
is.inus("(S <-> JSR)*(JSR <-> s)", fsdata, csf.info = TRUE)

# A csf with multiple identical outcomes.
is.inus("(S*C + V <-> JSR)*(R + P <-> JSR)", fsdata, csf.info = TRUE)
}

