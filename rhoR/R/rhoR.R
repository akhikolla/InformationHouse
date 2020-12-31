#' @title rhoR: A package for computing rho
#' @name rhoR
#' 
#' @description 
#' Rho is used to test the generalization of inter rater reliability (IRR)
#' statistics, in this case Cohen's Kappa.
#' 
#' Rho is a Monte Carlo rejective method of interrater reliability statistics, implemented here for Cohen's Kappa. Rho constructs a collection of data sets in which kappa is below a specified threshold, and computes the empirical distribution on kappa based on the specified sampling procedure. Rho returns the percent of the empirical distribution greater than or equal to an observed kappa. As a result, Rho quantifies the type 1 error in generalizing from an observed test set to a true value of agreement between two raters.
#' 
#' Rho starts with an observed kappa value, calculated on a subset of a \code{\link{codeSet}}, known as an observed \code{\link[=getTestSet]{testSet}}, and a \emph{kappa threshold} which indicates what is considered significant agreement between raters.
#' 
#' It then generates a collection of fully-coded, simulated 
#' \code{\link[=getTestSet]{codeSets}} (ScS), further described in 
#' \code{\link{createSimulatedCodeSet}}, all of which have a kappa value below the kappa 
#' threshold and similar properties as the original \code{\link{codeSet}}.
#' 
#' Then, kappa is calculated on a \code{\link[=getTestSet]{testSet}} sampled from each of the ScSs in the
#' collection to create a null hypothesis distribution. These \code{\link[=getTestSet]{testSets}} mirror the observed \code{\link[=getTestSet]{testSet}} in their size and sampling method.  How these \code{\link[=getTestSet]{testSets}} are sampled is futher described in \code{\link[=getTestSet]{testSet}}.
#' 
#' The null hypothesis is that the observed \code{\link[=getTestSet]{testSet}}, was sampled from a data set, which, if both raters were to code in its entirety, would result in a level of agreement below the kappa threshold.
#' 
#' For example, using an alpha level of 0.05, if the observed kappa is greater than 95 percent of the kappas in the null hypothesis distribution, the null hypothesis is rejected. Then one can conclude that the two raters would have acceptable agreement had they coded the entire data set.
#' 
#' @section rho: 
#' Use \code{\link{rho}} \code{\link{rhoK}} \code{\link{rhoSet}} \cr \code{\link{rhoCT}}
#' 
#' @section kappa: 
#' Use \code{\link{kappa}} \code{\link{kappaSet}} \cr \code{\link{kappaCT}}
#' 
#' @section rhoMin: 
#' Use \code{\link{rhoMin}}
#'   
#' @import stats
#' @import utils
#' @importFrom methods is
#' @importFrom Rcpp sourceCpp
#' 
#' @useDynLib rhoR
NULL

# @title additional attribtues
# @description Additional attribtues appending to rating.sets
additional_attributes <- c("agreement", "baserate", "kappa")
