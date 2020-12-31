#' ptsuite: Package for Pareto Parameter Estimation
#'
#' The ptsuite package provides functions to estimate parameters of Type I
#' Pareto data. The details of the methods used and the equations implemented
#' are given in the vignette. The package contains the following functions:
#'
#' @section Estimator Functions:
#' \itemize{
#'     \item Maximum Likelihood Estimator
#'     \{\code{\link{alpha_mle}}\}
#'     \item Weighted Least Squares Estimator
#'     \{\code{\link{alpha_wls}}\}
#'     \item Hill's Estimator
#'     \{\code{\link{alpha_hills}}\}
#'     \item Method of Moments Estimator
#'     \{\code{\link{alpha_moment}}\}
#'     \item Method of Percentiles Estimator
#'     \{\code{\link{alpha_percentile}}\}
#'     \item Method of Modified Percentiles Estimator
#'     \{\code{\link{alpha_modified_percentile}}\}
#'     \item Method of Geometric Percentiles
#'     \{\code{\link{alpha_geometric_percentile}}\}
#'     \item Least Squares Estimator
#'     \{\code{\link{alpha_ls}}\}
#' }
#'
#' @section Other Functions:
#' \itemize{
#'     \item Generate Pareto Data
#'     \{\code{\link{generate_pareto}}\}
#'     \item Estimates from all estimators
#'     \{\code{\link{generate_all_estimates}}\}
#'     \item Q-Q Plot to test for Pareto Distribution
#'     \{\code{\link{pareto_qq_test}}\}
#'     \item Pareto Test
#'     \{\code{\link{pareto_test}}\}
#' }
#' 
#' @section References:
#' Newman MEJ (2005). "Power Laws, Pareto Distributions And Zipf's Law." 
#' Contemporary Physics, 46, 323-351.
#' 
#' Nair J, Wierman A, Zwart B (2019). "The Fundamentals Of Heavy Tails: 
#' Properties, Emergence, And Identification." 
#' http://users.cms.caltech.edu/ adamw/heavytails.html.
#' 
#' Pokorna M (2016). Estimation and Application of the Tail Index. 
#' Bachelor's thesis, Charles University in Prague, Faculty of Social 
#' Sciences, Institute of Economic Studies.
#' 
#' Hill B (1975). "A Simple General Approach To Inference About The Tail 
#' Of A Distribution."The Annals of Statistics, 3(5), 1163-1174.
#' 
#' Rytgaard M (1990). "Estimation In The Pareto Distribution." ASTIN
#' Bulletin: The Journal Of The IAA, 20(2), 201-216.
#' 
#' Brazauskas V, Serfling R (2000). "Robust and Efficient Estimation Of 
#' The Tail Index Of A Single-Parameter Pareto Distribution." North 
#' American Actuarial Journal, 4, 12-27.
#'
#' Bhatti SH, Hussain S, Ahmad T, Aslam M, Aftab M, Raza MA (2018). "Efficient
#' estimation of Pareto model: Some modified percentile estimators." 
#' PLoS ONE, 13(5), 1-15.
#' 
#' Zaher HM, El-Sheik AA, El-Magd NATA (2014). "Estimation of Pareto 
#' Parameters Using a Fuzzy Least-Squares Method and Other Known Techniques 
#' with a Comparison." British Journal of Mathematics & Computer Science, 
#' 4(14), 2067-2088.
#' 
#' Gulati S, Shapiro S (2008). "Goodness-of-Fit Tests for Pareto 
#' Distribution." In F Vonta (ed.), Statistical Models and Methods for 
#' Biomedical and Technical Systems, chapter 19, pp. 259-274. Birkhauser 
#' Basel. ISBN 978-0-8176-4619-6. doi:10.1007/978-0-8176-4619-6.
#' 
#' @docType package
#' @name ptsuite-package
#' @useDynLib ptsuite
#' @import stats
#' @import graphics
#' @import utils
#' @importFrom Rcpp sourceCpp
NULL
