\name{cna}
\alias{cna}
\alias{cscna}
\alias{mvcna}
\alias{fscna}
\alias{print.cna}


\title{Perform Coincidence Analysis}

\description{
The \code{cna} function performs Coincidence Analysis to identify atomic solution formulas (asf) consisting of minimally necessary
disjunctions of minimally sufficient conditions of all outcomes in the data
and combines the recovered asf to complex solution formulas (csf) representing multi-outcome structures, e.g. common-cause and/or
causal chain structures.
}

\usage{
cna(x, type, ordering = NULL, strict = FALSE, con = 1, cov = 1, con.msc = con,
    notcols = NULL, rm.const.factors = TRUE, rm.dup.factors = TRUE,  
    maxstep = c(3, 4, 10), inus.only = only.minimal.msc && only.minimal.asf, 
    only.minimal.msc = TRUE,  only.minimal.asf = TRUE, 
    maxSol = 1e6, suff.only = FALSE, 
    what = if (suff.only) "m" else "ac", cutoff = 0.5, 
    border = c("down", "up", "drop"), details = FALSE, 
    acyclic.only = FALSE, cycle.type = c("factor", "value"))
cscna(...)
mvcna(...)
fscna(...)

\method{print}{cna}(x, what = x$what, digits = 3, nsolutions = 5, 
      details = x$details, show.cases = NULL, inus.only = x$inus.only, 
      acyclic.only = x$acyclic.only, cycle.type = x$cycle.type, 
      verbose = FALSE, ...)
}

\arguments{
  \item{x}{Data frame or \code{configTable} (as output by \code{\link{configTable}}).}
  \item{type}{Character vector specifying the type of \code{x}: \code{"cs"} (crisp-set), \code{"mv"} (multi-value),  or \code{"fs"} (fuzzy-set).}
  \item{ordering}{List of character vectors specifying the causal ordering of
        the factors in \code{x}.}
  \item{strict}{Logical; if \code{TRUE}, factors on the same level of the causal
        ordering are \emph{not} potential causes of each other; if \code{FALSE}, factors on the same level \emph{are} potential causes of each other.}
  \item{con}{Numeric scalar between 0 and 1 to set the minimum consistency threshold every minimally sufficient condition (msc), atomic solution formula (asf), and complex solution formula (csf) must satisfy. (See also the argument \code{con.msc} below).}
  \item{cov}{Numeric scalar between 0 and 1 to set the minimum coverage threshold every asf and csf must satisfy.}
  \item{con.msc}{Numeric scalar between 0 and 1 to set the minimum consistency threshold every msc must satisfy. Overrides \code{con} for msc and, thereby, allows for imposing a consistency threshold on msc that differs from the threshold \code{con} imposes on asf and csf. Defaults to \code{con}.}
  \item{maxstep}{Vector of three integers; the first specifies the maximum number of conjuncts in each disjunct of an asf, the second specifies the maximum number of disjuncts in an asf, the third specifies the maximum \emph{complexity} of an asf. The complexity of an asf is %an integer defined to be the sum of the number of conjuncts in all of its disjuncts, i.e. 
  the total number of exogenous factors in the asf. Default: \code{c(3,4,10)}.}
  \item{inus.only}{Logical; if \code{TRUE}, only disjunctive normal forms that are free of redundancies are retained as asf (see also \code{\link{is.inus}}). \cr Defaults to \code{only.minimal.msc && only.minimal.asf}.}
  \item{only.minimal.msc}{Logical; if \code{TRUE} (the default), only minimal conjunctions are retained as msc. If \code{FALSE}, sufficient conjunctions are not required to be minimal.}
  \item{only.minimal.asf}{Logical; if \code{TRUE} (the default), only minimal disjunctions are retained as asf. If \code{FALSE}, necessary disjunctions are not required to be minimal.}
  \item{maxSol}{Maximum number of asf calculated.}
  \item{suff.only}{Logical; if \code{TRUE}, the function only searches for msc and not for asf and csf.}
  \item{notcols}{Character vector of factors to be negated in \code{x}. If \code{notcols = "all"}, all factors in \code{x} are negated.}
  \item{rm.const.factors, rm.dup.factors}{Logical; if \code{TRUE} (default), factors with constant values are removed and all but the first of a set of duplicated factors are removed. These parameters are passed to \code{\link{configTable}}.}
  \item{what}{Character string specifying what to print; \code{"t"} for the configuration table, \code{"m"} for msc, \code{"a"} for asf, \code{"c"} for csf, and \code{"all"} for all. Defaults to \code{"ac"} if \code{suff.only = FALSE}, and to \code{"m"} otherwise.}
  \item{cutoff}{Minimum membership score required for a factor to count as instantiated in the data and to be integrated in the analysis. Value in the unit interval [0,1]. The default cutoff is 0.5. Only meaningful if \code{type = "fs"}.}
  \item{border}{Character vector specifying whether factors with membership scores equal to \code{cutoff} are rounded up (\code{"up"}), rounded down (\code{"down"}) or dropped from the analysis (\code{"drop"}). Only meaningful if \code{type="fs"}. }
  \item{details}{Either \code{TRUE}/\code{FALSE}, or a character vector with possible elements \code{"exhaustiveness"}, \code{"faithfulness"}, \code{"coherence"}, \code{"redundant"}, \code{"cyclic"}. The strings can also be abbreviated, e.g.  \code{"e"} or \code{"exh"} for \code{"exhaustiveness"}, etc.}
  \item{acyclic.only}{Logical; if \code{TRUE}, csf featuring a cyclic substructure are not returned. \code{FALSE} by default.}
  \item{cycle.type}{Character string specifying what type of cycles to be detected: \code{"factor"} (the default) or \code{"value"} (cf. \code{\link{cyclic}}).}
  \item{verbose}{Logical; if \code{TRUE}, some details on the csf building process are printed. \code{FALSE} by default.}
  \item{digits}{Number of digits to print in consistency, coverage, exhaustiveness, faithfulness, and coherence scores.}
  \item{nsolutions}{Maximum number of msc, asf, and csf to print. Alternatively, \code{nsolutions = "all"} will print all solutions.}
  \item{show.cases}{Logical; if \code{TRUE}, the \code{configTable}'s attribute \dQuote{cases}
        is printed. \cr See \code{\link{print.configTable}}}
  \item{\dots}{
        In \code{cscna}, \code{mvcna}, \code{fscna}: any formal argument of \code{cna} except \code{type}.
        In \code{print.cna}: arguments passed to other \code{print}-methods.}
}

\details{
The \strong{first input \code{x}} of the \code{cna} function is a data frame or a configuration table. To ensure that no misinterpretations of returned asf and csf can occur, users are advised to use only upper case letters as factor (column) names. Column names may contain numbers, but the first sign in a column name must be a letter. Only ASCII signs should be used for column and row names.

The \code{cna} function must be told what \bold{type of data} \code{x} contains, unless \code{x} is a configuration table. In the latter case, the type of \code{x} is already defined. Data that feature factors taking values 1 or 0 only are called \emph{crisp-set}, in which case the \code{type} argument takes its default value \code{"cs"}. If the data contain at least one factor that takes more than two values, e.g. \{1,2,3\}, the data count as \emph{multi-value}, which is indicated by \code{type = "mv"}. Data featuring at least one factor taking real values from the interval [0,1] count as \emph{fuzzy-set}, which is specified by \code{type = "fs"}. (Note that mixing multi-value and fuzzy-set factors in one analysis is not (currently) supported). To abbreviate the specification of the data type using the \code{type} argument, the functions \code{cscna(x, ...)}, \code{mvcna(x, ...)}, and \code{fscna(x, ...)} are available as shorthands for \code{cna(x, type = "cs", ...)}, \code{cna(x, type = "mv", ...)}, and \code{cna(x, type = "fs", ...)}, respectively.

A data frame or configuration table \code{x} with a corresponding type specification is the only mandatory input of the \code{cna} function. In particular, the \code{cna} function does not need an input specifying which factor(s) in \code{x} are endogenous, it tries to infer that from the data. Still, when prior causal knowledge about an investigated process is available, \code{cna} can be prohibited from treating certain factors as potential causes of other factors by means of the argument \code{ordering}. If specified, that argument defines a \bold{causal ordering} for the factors in \code{x}. For example,  \code{ordering = list(c("A",} \code{ "B"), "C")} determines that C is causally located \emph{after} A and B, meaning that C is \emph{not} a potential cause of A and B. In consequence, \code{cna} only checks whether values of A and B can be modeled as causes of values of C; the test for a causal dependency in the other direction is skipped. If the argument \code{ordering} is not specified or if it is given the \code{NULL} value (which is the argument's default value), \code{cna} searches for dependencies between all factors in \code{x}. An \code{ordering} does not need to explicitly mention all factors in \code{x}. If only a subset of the factors are included in the \code{ordering}, the non-included factors are entailed to be causally before the included ones. Hence, \code{ordering = list("C")}, for instance, means that C is causally located after all other factors in \code{x}, meaning that C is the ultimate outcome of the structure under scrutiny. 

The argument \code{strict} determines whether the elements of one level in an ordering can be causally related or not. For example, if \code{ordering = list(c("A", "B"), "C")} and \code{strict = TRUE}, then A and B---which are on the same level of the ordering---are excluded to be causally related and \code{cna} skips corresponding tests. By contrast, if \code{ordering = list(c("A", "B"), "C")} and \code{strict = FALSE}, then \code{cna} also searches for dependencies among A and B. The default is \code{strict} \code{ = FALSE}. If the user knows prior to the analysis that the data contain exactly one endogenous factor E and that the remaining exogenous factors are mutually causally independent, the appropriate function call should feature \code{cna(..., ordering = list("E"), strict = TRUE,...)}.




If no causal ordering is provided, all factor values in \code{x} are treated as potential outcomes; more specifically, in case of \code{"cs"} and \code{"fs"} data, \code{cna} tests for all factors whether their presence (i.e. them taking the value 1) can be modeled as an outcome, and in case of  \code{"mv"} data, \code{cna} tests for all factors whether any of their possible values can be modeled as an outcome. That is done by searching for redundancy-free Boolean functions (in disjunctive normal form) that account for the behavior of an outcome in accordance with \code{cna}'s core model fit parameters of \bold{consistency and coverage} (for details see the \pkg{cna} package vignette or Ragin 2006). First, \code{cna} identifies all minimally sufficient conditions (msc) that meet the threshold given by the consistency threshold \code{con.msc} (resp. \code{con}, if \code{con.msc = con}) for each factor in \code{x}. Then, these msc are disjunctively combined to minimally
necessary conditions that meet the coverage threshold given by \code{cov} such that the whole disjunction meets the solution consistency threshold given by \code{con}. The resulting expressions are the atomic solution formulas (asf) for every factor value that can be modeled as outcome. The default value for \code{con.msc}, \code{con}, and \code{cov} is 1. 



The \code{cna} function builds its models in four stages using a \strong{\emph{bottom-up search algorithm}} (see Baumgartner and Ambuehl 2018). \describe{\item{First stage}{On the basis of the \code{ordering}, the algorithm builds a set of potential outcomes \strong{O} from the  factors in \code{x}.}

\item{Second stage}{The algorithm
checks whether single factor values, e.g. A, b, C, (where "A" stands for "A=1" and "b" for "B=0") or D=3, E=2, etc., (whose membership scores, in case of \code{"fs"} data, meet \code{cutoff} in at least one case) are sufficient for a potential outcome in \strong{O} (where a factor value counts as sufficient iff it meets the threshold given by \code{con.msc}). Next, conjuncts of two factor values, e.g. A*b, A*C, D=3*E=2 etc., (whose membership scores, in case of \code{"fs"} data, meet \code{cutoff} in at least one case) are tested for sufficiency. Then, conjuncts of three factors, and so on. Whenever a conjunction (or a single factor value) is found to be sufficient, all supersets of that conjunction contain redundancies and are, thus, not considered for the further analysis. The result is a set of msc for every outcome in \strong{O}. To recover certain target structures in cases of noisy data, it may be useful to allow \code{cna} to also consider sufficient conditions for further analysis that are not minimal. This can be accomplished by setting \code{only.minimal.msc} to \code{FALSE}. A concrete example illustrating the utility of \code{only.minimal.msc} is provided in the \dQuote{Examples} section below. (The ordinary user is advised not to change the default value of this argument.)} 

\item{Third stage}{Minimally necessary disjunctions are built for each outcome in \strong{O} by first testing whether single msc are necessary, then disjunctions of two msc, then of three, etc. (where a disjunction of msc counts as necessary iff it meets the threshold given by \code{cov}). Whenever a disjunction of msc (or a single msc) is found to be necessary, all supersets of that disjunction contain redundancies and are, thus, excluded from the further analysis. Finally, all and only those disjunctions of msc that meet both \code{cov} and \code{con} are issued as redundancy-free \strong{atomic solution formulas} (asf). To recover certain target structures in cases of noisy data, it may be useful to allow \code{cna} to also consider necessary conditions for further analysis that are not minimal. This can be accomplished by setting \code{only.minimal.asf} to \code{FALSE}, in which case \emph{all} disjunctions of msc reaching the con and cov thresholds will be returned. (The ordinary user is advised not to change the default value of this argument.)  %A concrete example illustrating the purpose of \code{only.minimal.msc} is provided in the example section below.

As the combinatorial search space for asf is potentially too large to be exhaustively scanned in reasonable time, the argument \bold{\code{maxstep}} allows for setting an upper bound for the complexity of the generated asf. \code{maxstep} takes a vector of three integers \code{c(i, j, k)} as input, entailing that the generated asf have maximally \code{j} disjuncts with maximally \code{i} conjuncts each and a total of maximally \code{k} factor values (\code{k} is the maximal complexity). The default is \code{maxstep = c(3, 4, 10)}.

Note that when the data feature noise due to uncontrolled background influences the default \code{con} and \code{cov} thresholds of 1 will often not yield any asf. In such cases, \code{con} and \code{cov} may be set to suitable values in the interval [0.7, 1]. \code{con} and \code{cov} should neither be set too high, in order to avoid overfitting, nor too low, in order to avoid underfitting. The \strong{overfitting danger} is severe in causal modeling with CNA (and configurational causal modeling more generally). For a discussion of this problem see Parkkinen and Baumgartner (2020), who also introduce a procedure for robustness assessment that explores all threshold settings in a given interval---in an attempt to reduce both the over- and underfitting. }

\item{Fourth stage}{If \code{cna} finds asf, it builds \strong{complex solution formulas} (csf) from those asf. This is done in a stepwise manner as follows. First, all logically possible conjunctions featuring one asf of every outcome are built. Second, if \code{inus.only = TRUE}, the solutions resulting from step 1 are freed of structural redundancies (cf. Baumgartner and Falk 2019), and  tautologous and contradictory solutions as well as solutions with partial structural redundancies and constant factors are eliminated (cf. \code{\link{is.inus}}). Third, if \code{acyclic.only = TRUE}, solutions with cyclic substructures are eliminated. Fourth, for those solutions that were modified in the previous steps, consistency and coverage are re-calculated and solutions that no longer reach \code{con} or \code{cov} are eliminated. The remaining solutions are returned as csf. (See also \code{\link{csf}}.)}
}


The \bold{default output} of \code{cna} lists asf and csf, ordered by complexity and the product of consistency and coverage. It provides the consistency and coverage scores of each solution, a complexity score, which corresponds to the number of exogenous factor values in a solution, and a column \dQuote{\code{inus}} indicating whether a solution has INUS form, meaning whether it is redundancy-free as required by the \emph{INUS-theory} of causation (Mackie 1974, ch. 3; Baumgartner and Falk 2019). If \code{inus.only = TRUE}, all solutions automatically have INUS form, but if \code{only.minimal.msc} or 
\code{only.minimal.asf} are set to \code{FALSE}, non-INUS solutions may also be returned. 

Apart from the standard solution attributes, \code{cna} can calculate a number of  \bold{further solution attributes}: \code{exhaustiveness}, \code{faithfulness}, \code{coherence}, \code{redundant}, and \code{cyclic} all of which are recovered by setting \code{details} to its non-default value \code{TRUE} or to a character vector specifying the attributes to be calculated. 
These attributes require explication (see also the package vignette):
\itemize{


\item \code{exhaustiveness} and \code{faithfulness} are two measures of model fit that quantify the degree of correspondence between the configurations that are, in principle, compatible with a solution and the configurations contained in the data from which that solution is derived. %Roughly, exhaustiveness is high when \emph{all} or most configurations \emph{compatible} with a solution are in the data, whereas faithfulness is high when \emph{no} or only few configurations that are \emph{incompatible} with a solution are in the data. More specifically,
\itemize{\item \code{exhaustiveness} amounts to the ratio of the number of configurations in the data that are compatible with a solution to the number of configurations in total that are compatible with a solution.}\itemize{ \item \code{faithfulness} amounts to the ratio of the number of configurations in the data that are compatible with a solution to the total number of configurations in the data.} High exhaustiveness and faithfulness means that the configurations in the data are all and only the configurations that are compatible with the solution. Low exhaustiveness and/or faithfulness means that the data do not contain all configurations compatible with the solution and/or the data contain many configurations not compatible with the solution. In general, solutions with higher exhaustiveness and faithfulness scores are preferable over solutions with lower scores because they are better supported by the evidence in the data.

\item  \code{coherence} measures the degree to which the asf combined in a csf cohere, i.e. are instantiated together in the data rather than independently of one another. For more details see \code{\link{coherence}}. 

\item \code{redundant} determines whether a csf contains structurally redundant proper parts. A csf with \code{redundant = TRUE} should not be causally interpreted. If \code{inus.only = TRUE}, all csf are free of structural redundancies.  For more details see \code{\link{redundant}}. 

\item \code{cyclic} determines whether a csf contains a cyclic substructure. For more details see \code{\link{cyclic}}.
}



 

The argument \code{notcols} is used to calculate asf and csf
for \bold{negative outcomes} in data of \code{type} \code{"cs"} and \code{"fs"} (in \code{"mv"} data \code{notcols} has no meaningful interpretation and, correspondingly, issues an error message). If \code{notcols = "all"}, all factors in \code{x} are negated,
i.e. their membership scores i are replaced by 1-i. If \code{notcols} is given a character vector 
of factors in \code{x}, only the factors in that vector are negated. For example, \code{notcols = c("A", "B")}
determines that only factors A and B are negated. The default is no negations, i.e. \code{notcols = NULL}.

\code{suff.only} is applicable whenever a complete \code{cna} analysis cannot be performed for reasons of computational complexity. In such a case, \code{suff.only = TRUE} forces \code{cna} to stop the analysis after the identification of msc, which will normally yield results even in cases when a complete analysis does not terminate. In that manner, it is possible to shed at least some light on the dependencies among the factors in \code{x}, in spite of an incomputable solution space.

\code{rm.const.factors} and \code{rm.dup.factors} are used to determine the handling of \bold{constant factors}, i.e. factors with constant values in all cases (rows) in \code{x}, and of \bold{duplicated factors}, i.e. factors that take identical value distributions in all cases in \code{x}. If \code{rm.const.factors = TRUE}, which is the default value, constant factors are removed from the data prior to the analysis, and if \code{rm.dup.factors = TRUE} (the default) all but the first of a set of duplicated factors are removed. From the perspective of configurational causal modeling, factors with constant values in all cases can neither be modeled as causes nor as outcomes; therefore, they can be removed prior to the analysis. Factors that take identical values in all cases cannot be distinguished configurationally, meaning they are one and the same factor as far as configurational causal modeling is concerned. Therefore, only one factor of a set of duplicated factors is standardly retained by \code{cna}.

The argument \code{what} can be specified both for the \code{cna} and the \code{print()}
function. It regulates what items of the output of \code{cna} are printed. If
\code{what} is given the value \dQuote{\code{t}}, the configuration table is printed; if
it is given an \dQuote{\code{m}}, the msc are printed; if it is given an \dQuote{\code{a}},
the asf are printed; if it is given a \dQuote{\code{c}}, the csf are printed.
\code{what = "all"} or \code{what = "tmac"} determine that all output items are
printed. Note that \code{what} has no effect on the computations that are performed when executing \code{cna}; it only determines how the result is printed.
The default output of \code{cna} is \code{what = "ac"}. It first returns the implemented ordering. Second, the top 5 asf and, third, the top 5 csf are reported, along with an indication of how many solutions in total exist. To print all msc, asf, and csf, the corresponding functions in \code{\link{condTbl}} should be used. %If csf are the same as asf, this is indicated by "Same as asf". 
In case of \code{suff.only = TRUE}, \code{what} defaults to \code{"m"}. msc are printed with an attribute \code{minimal} specifying whether a sufficient condition is minimal as required by the INUS-theory of causation. If \code{inus.only = TRUE}, all msc are minimal by default. 

\code{cna} only includes factor configurations in the analysis that are actually instantiated in the data. The argument \code{cutoff} determines the minimum membership score required for a factor or a combination of factors to count as instantiated. It takes values in the unit interval [0,1] with a default of 0.5. \code{border} specifies whether configurations with membership scores equal to \code{cutoff} are rounded up (\code{border = "up"}), rounded down (\code{border = "down"}), which is the default, or dropped from the analysis (\code{border = "drop"}).

The arguments \code{digits}, \code{nsolutions}, and \code{show.cases} apply to the \strong{\code{print()} method}, which takes an object of class \dQuote{cna} as first input. \code{digits} determines how many digits of consistency, coverage, coherence, exhaustiveness, and faithfulness scores
are printed, while \code{nsolutions} fixes the number of conditions and solutions
to print. \code{nsolutions} applies separately to minimally sufficient conditions,
atomic solution formulas, and complex solution formulas. \code{nsolutions = "all"} recovers all minimally sufficient conditions, atomic and complex solution formulas. \code{show.cases} is applicable if the \code{what} argument is given the value \dQuote{\code{t}}. In that case, \code{show.cases = TRUE} yields a configuration table featuring a \dQuote{cases} column, which assigns cases to configurations.

The option \dQuote{spaces} controls how the conditions are rendered. The current setting is queried by typing \code{getOption("spaces")}. The option specifies characters that will be printed with a space before and after them. The default is \code{c("<->","->","+")}. A more compact output is obtained with \code{option(spaces = NULL)}. 
}

\value{
\code{cna} returns an object of class \dQuote{cna}, which amounts to a list with the following elements:

\tabular{rl}{
\code{call}: \tab the executed function call\cr
\code{x}:\tab the processed data frame or configuration table\cr
\code{ordering}:\tab the implemented ordering\cr
\code{configTable}: \tab the object of class \dQuote{configTable}, as input to \code{cna}\cr
\code{configTable_out}: \tab the object of class \dQuote{configTable}, after modification according to \code{notcols}\cr
\code{solution}: \tab the solution object, which itself is composed of lists exhibiting msc, asf, and csf for\cr\tab all factors in \code{x}\cr
\code{what}:\tab the values given to the \code{what} argument\cr
\code{details}:\tab the calculated solution attributes\cr
\code{...}:\tab plus additional list elements reporting the values given to the parameters \code{con}, \code{cov},\cr
\tab\code{con.msc}, \code{inus.only}, \code{acyclic.only}, and \code{cycle.type}. 
  }
}

\note{In the first example described below (in \emph{Examples}), the two resulting complex solution formulas represent a common cause structure and a causal chain, respectively. The common cause structure is graphically depicted in figure (a) below, the causal chain in figure (b).

\if{html}{\figure{structures3.png}{Causal Structures}}
\if{latex}{\figure{structures3.png}{options: width=13.5cm}}
}


\section{Contributors}{
Epple, Ruedi: development, testing\cr
Thiem, Alrik: testing
}



\references{
Basurto, Xavier. 2013. \dQuote{Linking Multi-Level Governance to Local Common-Pool 
Resource Theory using Fuzzy-Set Qualitative Comparative Analysis: Insights from 
Twenty Years of Biodiversity Conservation in Costa Rica.} \emph{Global Environmental Change} 23(3):573-87.

%Baumgartner, Michael. 2008. \dQuote{Regularity Theories Reassessed.}
%\emph{Philosophia} 36:327-354.

Baumgartner, Michael. 2009a. \dQuote{Inferring Causal Complexity.}
\emph{Sociological Methods & Research} 38(1):71-101.

Baumgartner, Michael and Mathias Ambuehl. 2018. \dQuote{Causal Modeling with Multi-Value and Fuzzy-Set Coincidence Analysis}. \emph{Political Science Research and Methods.} doi:10.1017/psrm.2018.45.

Baumgartner, Michael and Christoph Falk. 2019. \dQuote{Boolean Difference-Making: A Modern Regularity Theory of Causation}. \emph{The British Journal for the Philosophy of Science}. doi:10.1093/bjps/axz047.

Hartmann, Christof, and Joerg Kemmerzell. 2010. \dQuote{Understanding Variations 
in Party Bans in Africa.} \emph{Democratization} 17(4):642-65.
doi: 10.1080/13510347.2010.491189.

Krook, Mona Lena. 2010.
\dQuote{Women's Representation in Parliament: A Qualitative Comparative Analysis.}
\emph{Political Studies} 58(5):886-908.

Mackie, John L. 1974. \emph{The Cement of the Universe: A Study of Causation.} Oxford: Oxford University Press.

Parkkinen, Veli-Pekka and Michael Baumgartner. 2020. \dQuote{Robustness and Model Selection in Configurational Causal Modeling}, forthcoming in \emph{Sociological Methods & Research}. \cr url: https://people.uib.no/mba110/docs/rccm.pdf

Ragin, Charles C. 2006. \dQuote{Set Relations in Social Research: Evaluating Their Consistency and Coverage}. \emph{Political Analysis} 14(3):291-310.

Wollebaek, Dag. 2010.
\dQuote{Volatility and Growth in Populations of Rural Associations.}
\emph{Rural Sociology} 75:144-166.
}

\seealso{\code{\link{configTable}}, \code{\link{condition}}, \code{\link{cyclic}}, \code{\link{condTbl}}, \code{\link{selectCases}}, \code{\link{makeFuzzy}}, \code{\link{some}},  \code{\link{coherence}}, \code{\link{minimalizeCsf}}, \code{\link{randomConds}}, \code{\link{is.submodel}},  \code{\link{is.inus}},  \code{\link{redundant}}, \code{\link{full.ct}}, \code{\link{d.educate}}, \code{\link{d.women}}, \code{\link{d.pban}}, \code{\link{d.autonomy}}}

\examples{
# Ideal crisp-set data from Baumgartner (2009a) on education levels in western democracies
# ----------------------------------------------------------------------------------------
# Exhaustive CNA without constraints on the search space; print atomic and complex 
# solution formulas (default output).
cna.educate <- cna(d.educate)
cna.educate
# The two resulting complex solution formulas represent a common cause structure 
# and a causal chain, respectively. The common cause structure is graphically depicted 
# in (Note, figure (a)), the causal chain in (Note, figure (b)).

# Print only complex solution formulas.
print(cna.educate, what = "c")

# Print only atomic solution formulas.
print(cna.educate, what = "a")

# Print only minimally sufficient conditions.
print(cna.educate, what = "m")

# Print only the configuration table.
print(cna.educate, what = "t")

# CNA with negations of the factors E and L.
cna(d.educate, notcols = c("E","L"))

# CNA with negations of all factors.
cna(d.educate, notcols = "all")

# Print msc, asf, and csf with all solution attributes.
cna(d.educate, what = "mac", details = TRUE)

# Add only the non-standard solution attributes "exhaustiveness" and  "faithfulness".
cna(d.educate, details = c("e", "f"))

# Print solutions without spaces before and after "+".
options(spaces = c("<->", "->" ))
cna(d.educate, details = c("e", "f"))

# Print solutions with spaces before and after "*".
options(spaces = c("<->", "->", "*" ))
cna(d.educate, details = c("e", "f"))

# Restore the default of the option "spaces".
options(spaces = c("<->", "->", "+"))


# Crisp-set data from Krook (2010) on representation of women in western-democratic parliaments
# ---------------------------------------------------------------------------------------------
# This example shows that CNA can infer which factors are causes and which ones
# are effects from the data. Without being told which factor is the outcome, 
# CNA reproduces the original QCA of Krook (2010).
\donttest{ana1 <- cna(d.women, details = c("e", "f"))
ana1}

# The two resulting asf only reach an exhaustiveness score of 0.438, meaning that
# not all configurations that are compatible with the asf are contained in the data
# "d.women". Here is how to extract the configurations that are compatible with 
# the first asf but are not contained in "d.women":
\donttest{library(dplyr)
setdiff(ct2df(selectCases(asf(ana1)$condition[1], full.ct(d.women))),
        d.women)
}

# Highly ambiguous crisp-set data from Wollebaek (2010) on very high volatility of 
# grassroots associations in Norway
# --------------------------------------------------------------------------------
# csCNA with ordering from Wollebaek (2010) [Beware: due to massive ambiguities, this analysis
# will take about 20 seconds to compute.]
\donttest{cna(d.volatile, ordering = list("VO2"), maxstep = c(6, 6, 16))}
              
# Using suff.only, CNA can be forced to abandon the analysis after minimization of sufficient 
# conditions. [This analysis terminates quickly.]
cna(d.volatile, ordering = list("VO2"), maxstep = c(6, 6, 16), suff.only = TRUE)

# Similarly, by using the default maxstep, CNA can be forced to only search for asf and csf
# with reduced complexity.
\donttest{cna(d.volatile, ordering = list("VO2"))}


# Multi-value data from Hartmann & Kemmerzell (2010) on party bans in Africa
# ---------------------------------------------------------------------------
# mvCNA with causal ordering that corresponds to the ordering in Hartmann & Kemmerzell 
# (2010); coverage cutoff at 0.95 (consistency cutoff at 1), maxstep at c(6, 6, 10).
cna.pban <- mvcna(d.pban, ordering = list(c("C","F","T","V"),"PB"), cov = .95,
                  maxstep = c(6, 6, 10), what = "all")
cna.pban

# The previous function call yields a total of 14 asf and csf, only 5 of which are 
# printed in the default output. Here is how to extract all 14 asf and csf.
asf(cna.pban)
csf(cna.pban)

# [Note that all of these 14 causal models reach better consistency and 
# coverage scores than the one model Hartmann & Kemmerzell (2010) present in their paper, 
# which they generated using the TOSMANA software, version 1.3: 
# T=0 + T=1 + C=2 + T=1*V=0 + T=2*V=0 <-> PB=1]
condTbl("T=0 + T=1 + C=2 + T=1*V=0 + T=2*V=0 <-> PB = 1", mvct(d.pban))

# Extract all minimally sufficient conditions.
msc(cna.pban)

# Alternatively, all msc, asf, and csf can be recovered by means of the nsolutions
# argument of the print function.
print(cna.pban, nsolutions = "all")

# Print the configuration table with the "cases" column.
print(cna.pban, what = "t", show.cases = TRUE)

# Build solution formulas with maximally 4 disjuncts.
\donttest{mvcna(d.pban, ordering = list(c("C","F","T","V"),"PB"), cov = .95, maxstep = c(4, 4, 10))

# Only print 2 digits of consistency and coverage scores.
print(cna.pban, digits = 2)

# Build all but print only two msc for each factor and two asf and csf.
print(mvcna(d.pban, ordering = list(c("C","F","T","V"),"PB"), cov = .95,
      maxstep = c(6, 6, 10), what = "all"), nsolutions = 2)

# Lowering the consistency instead of the coverage threshold yields further models with
# excellent fit scores; print only asf.
mvcna(d.pban, ordering = list(c("C","F","T","V"),"PB"), con = .93, what = "a",
      maxstep = c(6, 6, 10))

# Importing an ordering from prior causal knowledge is unnecessary for d.pban. PB  
# is the only factor in those data that could possibly be an outcome.
mvcna(d.pban, cov = .95, maxstep = c(6, 6, 10))
}

# Fuzzy-set data from Basurto (2013) on autonomy of biodiversity institutions in Costa Rica
# ---------------------------------------------------------------------------------------
# Basurto investigates two outcomes: emergence of local autonomy and endurance thereof. The 
# data for the first outcome are contained in rows 1-14 of d.autonomy, the data for the second
# outcome in rows 15-30. For each outcome, the author distinguishes between local ("EM",  
# "SP", "CO"),  national ("CI", "PO") and international ("RE", "CN", "DE") conditions. Here,   
# we first apply fsCNA to replicate the analysis for the local conditions of the endurance of 
# local autonomy.
dat1 <- d.autonomy[15:30, c("AU","EM","SP","CO")]
fscna(dat1, ordering = list("AU"), strict = TRUE, con = .9, cov = .9)

# The fsCNA model has significantly better consistency (and equal coverage) scores than the 
# model presented by Basurto (p. 580): SP*EM + CO <-> AU, which he generated using the 
# fs/QCA software.
fscond("SP*EM + CO <-> AU", dat1) # both EM and CO are redundant to account for AU

# If we allow for dependencies among the conditions by setting strict = FALSE, CNA reveals 
# that SP is a common cause of both AU and EM:
fscna(dat1, ordering = list("AU"), strict = FALSE, con = .9, cov = .9)

# Here is the analysis for the international conditions of autonomy endurance, which
# yields the same model as the one presented by Basurto (plus one model Basurto does not mention):
dat2 <- d.autonomy[15:30, c("AU","RE", "CN", "DE")]
fscna(dat2, ordering = list("AU"), con = .9, con.msc = .85, cov = .85)

# But there are other models (here printed with all solution attributes)
# that fare equally well.
fscna(dat2, ordering = list("AU"), con = .85, cov = .9, details = TRUE)

# Finally, here is an analysis of the whole data set, showing that across the whole period 
# 1986-2006, the best causal model of local autonomy (AU) renders that outcome dependent
# only on local direct spending (SP):
\donttest{fscna(d.autonomy, ordering = list("AU"), strict = TRUE, con = .85, cov = .9, 
      maxstep = c(5, 5, 11), details = TRUE)}

# Also build non-INUS solutions.
\donttest{asf(fscna(d.autonomy, ordering = list("AU"), strict = TRUE, con = .85, cov = .9, 
      maxstep = c(5, 5, 11), details = TRUE, inus.only = FALSE))}


# Highly ambiguous artificial data to illustrate exhaustiveness and acyclic.only
# ------------------------------------------------------------------------------
mycond <- "(D + C*f <-> A)*(C*d + c*D <-> B)*(B*d + D*f <-> C)*(c*B + B*f <-> E)"
dat1 <- selectCases(mycond)
\donttest{ana1 <- cna(dat1, details = c("e","cy"))
# There exist almost 2M csf. This is how to build the first 1076 of them, with 
# additional messages about the csf building process:
first.csf <- csf(ana1, verbose = TRUE)
# Most of these csf are compatible with more configurations than are contained in 
# dat1. Only 105 csf in first.csf are perfectly exhaustive (i.e. all compatible 
# configurations are contained in dat1):
subset(first.csf, exhaustiveness == 1)

# 1020 of the csf in first.csf contain cyclic substructures.
subset(first.csf, cyclic == TRUE)

# Here's how to only build acyclic csf:
ana2 <- cna(dat1, details = c("e","cy"), acyclic.only = TRUE)
csf(ana2, verbose = TRUE)
}

# Inverse search trials to assess the correctness of CNA
# ------------------------------------------------------
# 1. Ideal mv data, i.e. perfect consistencies and coverages, without data fragmentation.
\donttest{# Define the target and generate data on the target.
target <- "(A=1*B=2 + A=4*B=3 <-> C=1)*(C=4*D=1 + C=2*D=4 <-> E=4)"
dat1 <- allCombs(c(4, 4, 4, 4, 4)) 
dat2 <- selectCases(target, dat1, type = "mv")
# Analyze the simulated data with CNA.
test1 <- mvcna(dat2, maxstep = c(3,2,9))
# Check whether a correctness-preserving submodel of the target is among the 
# returned solutions. 
is.submodel(csf(test1)$condition, target)

# Same test as above with data fragmentation, i.e. with non-ideal data:
# only 100 of 472 observable configurations are actually
# observed. [Repeated runs will generate different data.]
dat3 <- some(dat2, n = 100, replace = TRUE)
test2 <- mvcna(dat3, maxstep = c(3,2,9))
is.submodel(csf(test2)$condition, target)

# 2. Fs data with imperfect consistencies (con = 0.8) and coverages (cov = 0.8); 
# about 150 cases (depending on the seed). Randomly generated target asf. 
# [Repeated runs will generate different targets and data. In some runs, no solutions
# are found.]
target <- randomAsf(full.ct(5), compl = c(2,3))
outcome <- as.vector(cna:::rhs(target))
# Simulate the data with con =  cov = 0.8.
dat1 <- allCombs(c(2, 2, 2, 2, 2)) - 1
dat2 <- some(configTable(dat1), n = 200, replace = TRUE)
dat3 <- makeFuzzy(ct2df(dat2), fuzzvalues = seq(0, 0.45, 0.01))
dat4 <- selectCases1(target, con = .8, cov = .8, type = "fs", dat3)
# Analyze the simulated data with CNA.
test3 <- fscna(dat4, ordering = list(outcome), strict = TRUE, con = .8, cov = .8)
# Check whether a correctness-preserving submodel of the target is among the 
# returned solutions. 
is.submodel(asf(test3)$condition, target)

# Same test as above with data fragmentation: only 80 of about 150 possible
# cases are actually observed. [Repeated runs will generate different data.]
dat5 <- some(dat4, n = 80, replace = TRUE)
test4 <- fscna(dat5, ordering = list(outcome), strict = TRUE, con = .8, cov = .8)
is.submodel(asf(test4)$condition, target)
}

# Illustration of only.minimal.msc = FALSE
# ----------------------------------------
# Simulate noisy data on the causal structure "a*B*d + A*c*D <-> E"
set.seed(1324557857)
mydata <- allCombs(rep(2, 5)) - 1
dat1 <- makeFuzzy(mydata, fuzzvalues = seq(0, 0.5, 0.01))
dat1 <- ct2df(selectCases1("a*B*d + A*c*D <-> E", con = .8, cov = .8, dat1))

# In dat1, "a*B*d + A*c*D <-> E" has the following con and cov scores:
as.condTbl(fscond("a*B*d + A*c*D <-> E", dat1))

# The standard algorithm of CNA will, however, not find this structure with
# con = cov = 0.8 because one of the disjuncts (a*B*d) does not meet the con
# threshold:
as.condTbl(fscond(c("a*B*d <-> E", "A*c*D <-> E"), dat1))
fscna(dat1, ordering=list("E"), strict = TRUE, con = .8, cov = .8)

# With the argument con.msc we can lower the con threshold for msc, but this does not
# recover "a*B*d + A*c*D <-> E" either:
cna2 <- fscna(dat1, ordering=list("E"), strict = TRUE, con = .8, cov = .8, con.msc = .78)
cna2
msc(cna2)

# The reason is that "A*c -> E" and "c*D -> E" now also meet the con.msc threshold and,
# therefore,  "A*c*D -> E" is not contained in the msc---because of violated minimality.
# In a situation like this, lifting the minimality requirement via 
# only.minimal.msc = FALSE allows CNA to find the intended target:
fscna(dat1, ordering=list("E"), strict=TRUE, con = .8, cov = .8, con.msc = .78,
      only.minimal.msc = FALSE)
}

