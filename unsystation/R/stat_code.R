#' A second-order stationarity of time series based on unsystematic sub-samples
#'
#' The function implements a stationarity test procedure, where the main statistic is obtained from measuring the difference in the second-order structure over pairs of randomly drawn intervals. Maximising the main statistics after AR Sieve bootstrap-based variance stabilisation, the test statistic is obtained which is reported along with the corresponding pair of intervals and the test outcome.
#' @param x input time series
#' @param M number of randomly drawn intervals
#' @param sig.lev significance level between \code{0} and \code{1}
#' @param max.scale  number of wavelet scales used for wavelet periodogram computation; \code{max.scale = NULL} activates the default choice (\code{max.scale = round(log(log(length(x), 2), 2))})
#' @param m minimum length of a random interval; \code{m = NULL} activates the default choice (\code{m = round(sqrt(length(x)))})
#' @param B bootstrap sample size
#' @param eps a parameter used for random interval generation, see the supplementary document of Cho (2016)
#' @param use.all if \code{use.all=TRUE}, all \code{M*M} pairs of random intervals are considered in test statistic computation; if \code{use.all=FALSE}, only \code{10*M} pairs are used; regardless, the whole \code{M*M} pairs are considered in test criterion generation
#' @param do.parallel number of copies of R running in parallel, if \code{do.parallel = 0}, \%do\% operator is used, see also \link{foreach}
#' @return 
#' \item{intervals}{a pair of intervals corresponding to the test statistic, exhibiting the most distinct second-order behaviour}
#' \item{test.stat}{test statistic}
#' \item{test.criterion}{test criterion}
#' \item{test.res}{if \code{test.res=TRUE}, the null hypothesis of stationarity is rejected at the given significance level}
#' @references
#' H. Cho (2016) A second-order stationarity of time series based on unsystematic sub-samples. Stat, vol. 5, 262-277.
#' @examples
#' \dontrun{
#' x <- rnorm(200)
#' unsys.station.test(x, M=1000)
#' }
#' @import Rcpp foreach doParallel iterators
#' @importFrom stats ar qnorm mad qt median
#' @useDynLib unsystation, .registration = TRUE
#' @export
unsys.station.test <- function(x, M=2000, sig.lev=.05, max.scale=NULL, m=NULL, B=200, eps=5, use.all=FALSE, do.parallel=0){

	T <- length(x)
	if(is.null(max.scale)) max.scale <- round(log(log(T, 2), 2))
	if(is.null(m)) m <- round(sqrt(T))
	if(do.parallel > 0){ cl <- parallel::makeCluster(do.parallel); doParallel::registerDoParallel(cl) }
	`%mydo%` <- ifelse(do.parallel > 0, `%dopar%`, `%do%`)
	
	top.cand0 <- bottom.cand0 <- NULL
	y.mat <- matrix(0, ncol=max.scale, nrow=T)
	for(k in 1:max.scale){
	  y.mat[, k] <- y <- func_coef(x, -k)^2
	  ref <- sort(y, decreasing=TRUE, index.return=TRUE)
	  top.cand0 <- c(top.cand0, setdiff(ref$ix[1:(eps + 2*2^k)], c(1:(2^k), (T-2^k+1):T))[1:eps])
	  bottom.cand0 <- c(bottom.cand0, setdiff(ref$ix[T:(T-eps-2*2^k+1)], c(1:(2^k), (T-2^k+1):T))[1:eps])
	}
	top.cand <- base::sample(top.cand0, M, replace=TRUE); bottom.cand <- base::sample(bottom.cand0, M, replace=TRUE)

	fr <- funcRes(y.mat, M, m, rep(1/T, T), top.cand, bottom.cand, apply(y.mat, 2, mean))

	ind <- which((!duplicated(fr$res[, 5])) & fr$res[, 5] > 0)
	if(use.all) ref <- ind else{
	  ref <- ind[sort(fr$res[ind, 4+max.scale+1], decreasing=TRUE, index.return=TRUE)$ix[1:(10*M)]]
	}
	se.mat <- fr$res[ref, 1:4]
	I <- length(ind); R <- length(ref)

	arx <- stats::ar(x, order.max=log(T), method='yw')
	coef <- arx$ar
	ep <- arx$resid[!is.na(arx$resid)];
	sig <- stats::mad(ep)
	ep <- ep[abs(ep-stats::median(ep)) < sig*stats::qt(1-.005, 10)]
	ep <- ep-mean(ep)
	if(length(coef)==0) boot.x <- matrix(base::sample(ep, B*T, replace=TRUE), ncol=B) else boot.x <- funcSimX(coef, matrix(base::sample(ep, (T+length(coef))*B, replace=TRUE), ncol=B))

	b <- 0
	null.stat <- foreach::foreach(b=iterators::iter(1:B), .combine=rbind, .packages=c('Rcpp', 'RcppArmadillo', 'unsystation')) %mydo% {
		bx <- boot.x[, b]
		by.mat <- matrix(0, ncol=max.scale, nrow=T)
		for(k in 1:max.scale){
			by.mat[, k] <- func_coef(bx, -k)^2
		}
		tmp <- funcResVar(by.mat, se.mat, apply(by.mat, 2, mean))
		c(tmp)
	}
	stat <- abs(fr$res[ref, 4+1:max.scale])/funcApplyVar(null.stat, max.scale, R)

	k <- which.max(apply(stat, 1, max))
	intervals <- fr$res[ref[k], 1:4]
	test.stat <- max(stat[k,])
	test.criterion <- stats::qnorm(1-sig.lev/2/I/max.scale)
	test.res <- test.stat > test.criterion

	if(do.parallel > 0) parallel::stopCluster(cl)

	return(list(intervals=intervals, test.stat=test.stat, test.criterion=test.criterion, test.res=test.res))

}
