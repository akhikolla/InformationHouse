\name{cna-package}
\alias{cna-package}


\docType{package}

\title{
cna: A Package for Causal Modeling with Coincidence Analysis
}

\description{
\dfn{Coincidence Analysis} (CNA) is a configurational comparative method of causal data analysis that was first introduced for crisp-set (i.e. binary) data in Baumgartner (2009a, 2009b, 2013) and generalized for multi-value and fuzzy-set data in Baumgartner and Ambuehl (2018). The \pkg{cna} package reflects and implements the method's latest stage of development. 

CNA infers causal structures as defined by modern variants of Mackie's (1974) INUS-theory of causation%, or modern variants thereof 
(e.g. Grasshoff and May 2001; Baumgartner and Falk 2019) from empirical data. INUS structures have two characteristic features: \emph{conjunctivity}---complex causes only become operative when all of their components are co-instantiated (each of which, in isolation, is ineffective)---and \emph{disjunctivity}---effects can be brought about along alternative causal routes such that, when suppressing one route, the effect may still be produced via another one. 
Causal structures featuring conjunctivity or disjunctivity often do not exhibit linear dependencies between pairs of exogenous and endogenous factors, which gives rise to various problems when it comes to discovering them in data. Most importantly, there may not be any dependencies between an individual component X of a conjunctive cause and the corresponding effect Y as long as not all components are properly co-instantiated and alternative causal pathways to Y suppressed. It follows that X cannot be identified as an INUS cause of Y by searching for suitable pairwise dependencies between X and Y but only by embedding X in a complex Boolean structure over many factors and fitting that structure as a whole to the data. But the space of Boolean functions over even a handful of factors is vast. So, a method for INUS-discovery must find ways to purposefully and efficiently navigate in that vast space of possibilities.



Methods of data analysis custom-built for the handling of INUS structures have been developed in various disciplines---independently of one another. In epidemiology, for example, Rothman's (1976) sufficient-component cause model with its various further developments (e.g. VanderWeele and Robins 2009) is designed for that purpose. In biostatistics, the method targeting conjunctivity and disjunctivity is called \emph{logic regression} (LR; Ruczinski et al. 2003), which is implemented in the R package  \pkg{\href{https://cran.r-project.org/package=LogicReg}{LogicReg}}. In the social sciences, the corresponding methods are called \emph{configurational comparative methods}, with Qualitative Comparative Analysis (QCA; Ragin 2008) as its best known representative, which is implemented in the R packages \pkg{\href{https://cran.r-project.org/package=QCApro}{QCApro}} and \pkg{\href{https://cran.r-project.org/package=QCA}{QCA}}.

CNA has commonalities and differences with both QCA and LR. Like LR, CNA exclusively searches for redundancy-free Boolean functions of causally modeled outcomes, whereas QCA also considers functions with redundant elements (in so-called intermediate and conservative solutions). Like QCA, CNA's core parameters of model fit are consistency and coverage (Ragin 2006), while LR draws on standard regression analytic fit criteria (e.g. the residual sum of squares). Unlike QCA, CNA does not generate causal models from the top down by first building maximal Boolean dependency structures and then gradually eliminating redundant elements (using e.g. Quine-McCluskey optimization); rather, like LR, CNA builds causal models from the bottom up by gradually combining single factor values to complex dependency structures until the requested thresholds of model fit are met, such that the resulting models are automatically redundancy-free. As a direct consequence, contrary to QCA, CNA can handle data fragmentation (limited diversity) without resorting to counterfactual reasoning. Contrary to LR, CNA does not rely on a search heuristic but exhaustively scans the space of data-fitting Boolean functions within user-defined bounds for model complexity. Finally, unlike both QCA and LR, CNA allows for multiple effects (outcomes/endogenous factors), meaning it can analyze common-cause and causal chain structures---it not only groups causes conjunctively and disjunctively but also sequentially. 

The new functionalities provided by version 3.0 of the \pkg{cna} package mainly concern the construction of models with multiple effects, so-called complex solution formulas (\emph{csf}). In its default setting, the \code{csf()} function now automatically eliminates structural redundancies (Baumgartner and Falk 2019) and ensures that all issued \emph{csf} have INUS form. Various other functions have been adapted accordingly, for instance, \code{cna()} or \code{is.inus()}. The stand-alone \code{minimalizeCsf()} function is now dispensable---but kept in the package for backwards compatibility. Apart from these added functionalities, the main change in version 3.0 is terminological. In order to process data, the \code{cna()} function first produces a compact representation of the data, which in previous package versions was called a ``truth table''. That label, however, led to confusion with the QCA terminology where a very different type of object is also called a ``truth table''. To prevent that confusion, the compact data representation is now called a ``configuration table'' in the \pkg{cna} package. In consequence, a number of functions received new names, e.g. \code{truthTab()} became \code{configTable()}, \code{full.tt()} became \code{full.ct()}, or \code{tt2df()} is now called \code{ct2df()}. For backwards compatibility, all old names are kept as aliases in the package.    
The package vignette, which presents the theoretical background of CNA and introduces to causal modeling with \pkg{cna}, has been updated accordingly. 



}

\details{
    \tabular{ll}{
Package: \tab cna\cr
Type:    \tab Package\cr
Version: \tab 3.0.1\cr
Date:    \tab 2020-11-06\cr
License: \tab GPL (>= 2)\cr
    }
}
\author{
\bold{Authors}:\cr
Mathias Ambuehl\cr
\email{mathias.ambuehl@consultag.ch}

Michael Baumgartner\cr
Department of Philosophy\cr
University of Bergen\cr
\email{michael.baumgartner@uib.no}


\bold{Maintainer}:\cr
Mathias Ambuehl
}

\references{
Baumgartner, Michael. 2009a. \dQuote{Inferring Causal Complexity.}
\emph{Sociological Methods & Research} 38(1):71-101.

Baumgartner, Michael. 2009b. \dQuote{Uncovering Deterministic Causal Structures:
A Boolean Approach.} \emph{Synthese} 170(1):71-96.

Baumgartner, Michael. 2013. \dQuote{Detecting Causal Chains in Small-n Data.}
\emph{Field Methods} 25 (1):3-24.

Baumgartner, Michael and Mathias Ambuehl. 2018. \dQuote{Causal Modeling with Multi-Value and Fuzzy-Set Coincidence Analysis.} \emph{Political Science Research and Methods}. doi:10.1017/psrm.2018.45. 

Baumgartner, Michael and Christoph Falk. 2019. \dQuote{Boolean Difference-Making: A Modern Regularity Theory of Causation}. \emph{The British Journal for the Philosophy of Science}. doi:10.1093/bjps/axz047.

Baumgartner, Michael and Alrik Thiem. 2015. \dQuote{Identifying Complex Causal Dependencies in Configurational Data with Coincidence Analysis}, \emph{The R Journal} 7:176-184.

Grasshoff, G. and M. May. 2001. \dQuote{Causal Regularities}. In W. Spohn, M. Ledwig, and
M. Esfeld (Eds.), \emph{Current Issues in Causation}, pp. 85-114. Paderborn: Mentis.

Mackie, John L. 1974. \emph{The Cement of the Universe: A Study of Causation}. Oxford: Oxford University Press.

Ragin, Charles C. 1987. \emph{The Comparative Method}. Berkeley: University of California Press.

Ragin, Charles C. 2006. \dQuote{Set Relations in Social Research: Evaluating Their Consistency and Coverage}. \emph{Political Analysis} 14(3):291-310.

Ragin, Charles C. 2008. \emph{Redesigning Social Inquiry: Fuzzy Sets and Beyond}. Chicago:
University of Chicago Press.

Ruczinski, I., C. Kooperberg, and M. LeBlanc. 2003. \dQuote{Logic Regression}. \emph{Journal of Computational and Graphical Statistics} 12:475-511.

VanderWeele, T. J. and J.M. Robins. 2009. \dQuote{Minimal Sufficient Causation and Directed
Acyclic Graphs}. \emph{Ann. Statist.} 37:1437-1465.
}

\keyword{package}
