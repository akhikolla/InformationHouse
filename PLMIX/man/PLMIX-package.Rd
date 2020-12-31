\name{PLMIX-package}
\alias{PLMIX-package}
\alias{PLMIX}
\docType{package}
\title{
Bayesian Analysis of Finite Mixtures of Plackett-Luce Models for Partial Rankings/Orderings
}
\description{
The \pkg{PLMIX} package for R provides functions to fit and analyze finite mixtures of Plackett-Luce models for partial top rankings/orderings within the Bayesian framework. It provides MAP point estimates via EM algorithm and posterior MCMC simulations via Gibbs Sampling. It also fits MLE as a special case of the noninformative Bayesian analysis with vague priors. 

In addition to inferential techniques, the package assists other fundamental phases of a model-based analysis for partial rankings/orderings, by including functions for data manipulation, simulation, descriptive summary, model selection and goodness-of-fit evaluation.

Specific S3 classes and methods are also supplied to enhance the usability and foster exchange with other packages. Finally, to address the issue of computationally demanding procedures typical in
ranking data analysis, \pkg{PLMIX} takes advantage of a hybrid code linking the R environment with the C++ programming language.
}
\details{
The Plackett-Luce model is one of the most popular and frequently applied parametric distributions to analyze partial top rankings/orderings of a finite set of items. The present package allows to account for unobserved sample heterogeneity of partially ranked data with a model-based analysis relying on Bayesian finite mixtures of Plackett-Luce models. The package provides a suite of functions that covers the fundamental phases of a model-based analysis: 

\strong{Ranking data manipulation}
  \describe{
      \item{\code{\link{binary_group_ind}}}{Binary group membership matrix from the mixture component labels.}
      \item{\code{\link{freq_to_unit}}}{From the frequency distribution to the dataset of individual orderings/rankings.}
      \item{\code{\link{make_complete}}}{Random completion of partial orderings/rankings data.}
      \item{\code{\link{make_partial}}}{Censoring of complete orderings/rankings data.}
      \item{\code{\link{rank_ord_switch}}}{From rankings to orderings and vice-versa.}
      \item{\code{\link{unit_to_freq}}}{From the dataset of individual orderings/rankings to the frequency distribution.}
    }
\strong{Ranking data simulation}
  \describe{
   \item{\code{\link{rPLMIX}}}{Random sample from a finite mixture of Plackett-Luce models.}
    }
\strong{Ranking data description}
  \describe{
    \item{\code{\link{paired_comparisons}}}{Paired comparison frequencies.}
    \item{\code{\link{rank_summaries}}}{Summary statistics of partial ranking/ordering data.}
  }    
\strong{Model estimation}
  \describe{
   \item{\code{\link{gibbsPLMIX}}}{Bayesian analysis with MCMC posterior simulation via Gibbs sampling.}
   \item{\code{\link{label_switchPLMIX}}}{Label switching adjustment of the Gibbs sampling simulations.}
   \item{\code{\link{likPLMIX}}}{Likelihood evaluation for a mixture of Plackett-Luce models.}
   \item{\code{\link{loglikPLMIX}}}{Log-likelihood evaluation for a mixture of Plackett-Luce models.}
   \item{\code{\link{mapPLMIX}}}{MAP estimation via EM algorithm.}
   \item{\code{\link{mapPLMIX_multistart}}}{MAP estimation via EM algorithm with multiple starting values.}
    }
\strong{Class coercion and membership}
  \describe{
   \item{\code{\link{as.top_ordering}}}{Coercion into top-ordering datasets.}
   \item{\code{\link{gsPLMIX_to_mcmc}}}{From the Gibbs sampling simulation to an MCMC class object.}
   \item{\code{\link{is.top_ordering}}}{Test for the consistency of input data with a top-ordering dataset.}
    }
\strong{S3 class methods}
  \describe{
   \item{\code{\link{plot.gsPLMIX}}}{Plot of the Gibbs sampling simulations.}
   \item{\code{\link{plot.mpPLMIX}}}{Plot of the MAP estimates.}
   \item{\code{\link{print.gsPLMIX}}}{Print of the Gibbs sampling simulations.}
   \item{\code{\link{print.mpPLMIX}}}{Print of the MAP estimation algorithm.}
   \item{\code{\link{summary.gsPLMIX}}}{Summary of the Gibbs sampling procedure.}
   \item{\code{\link{summary.mpPLMIX}}}{Summary of the MAP estimation.} 
    }
\strong{Model selection}
  \describe{
   \item{\code{\link{bicPLMIX}}}{BIC value for the MLE of a mixture of Plackett-Luce models.}
   \item{\code{\link{selectPLMIX}}}{Bayesian model selection criteria.}
    }
\strong{Model assessment}
  \describe{
   \item{\code{\link{ppcheckPLMIX}}}{Posterior predictive diagnostics.}
   \item{\code{\link{ppcheckPLMIX_cond}}}{Posterior predictive diagnostics conditionally on the number of ranked items.}
    }
\strong{Datasets}
  \describe{
   \item{\code{\link{d_apa}}}{American Psychological Association Data (partial orderings).}
   \item{\code{\link{d_carconf}}}{Car Configurator Data (partial orderings).}
   \item{\code{\link{d_dublinwest}}}{Dublin West Data (partial orderings).}
   \item{\code{\link{d_gaming}}}{Gaming Platforms Data (complete orderings).}
   \item{\code{\link{d_german}}}{German Sample Data (complete orderings).}
   \item{\code{\link{d_nascar}}}{NASCAR Data (partial orderings).}
   \item{\code{\link{d_occup}}}{Occupation Data (complete orderings).}
   \item{\code{\link{d_rice}}}{Rice Voting Data (partial orderings).}
    }
Data have to be supplied as an object of class \code{matrix}, where missing positions/items are denoted with zero entries and Rank = 1 indicates the most-liked alternative. For a more efficient implementation of the methods, partial sequences with a single missing entry should be preliminarily filled in, as they correspond to complete rankings/orderings. In the present setting, ties are not allowed. Some quantities frequently recalled in the manual are the following: 
\describe{
    \item{\eqn{N}}{Sample size.}
    \item{\eqn{K}}{Number of possible items.}
    \item{\eqn{G}}{Number of mixture components.}
    \item{\eqn{L}}{Size of the final posterior MCMC sample (after burn-in phase).}
     }
}
\author{
Cristina Mollica and Luca Tardella

Maintainer: Cristina Mollica <cristina.mollica@uniroma1.it>
}
\references{
Mollica, C. and Tardella, L. (2017). Bayesian Plackett-Luce mixture models for partially ranked data. \emph{Psychometrika}, \bold{82}(2), pages 442--458, ISSN: 0033-3123, \url{http://dx.doi.org/10.1007/s11336-016-9530-0}.

Mollica, C. and Tardella, L. (2014). Epitope profiling via mixture modeling for ranked data. \emph{Statistics in Medicine}, \bold{33}(21), pages 3738--3758, ISSN: 0277-6715, \url{http://onlinelibrary.wiley.com/doi/10.1002/sim.6224/full}.
}
\concept{Ranking data}
\concept{Plackett-Luce model}
\concept{Mixture models}
\concept{Bayesian inference}
\concept{MAP estimation}
\concept{MCMC approximation}

