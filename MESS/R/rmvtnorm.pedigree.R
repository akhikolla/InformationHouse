#' Simulate residual multivariate Gaussian data from a polygenic model
#'
#' Simulates residual multivariate Gaussian response data from a pedigree where
#' the additive genetic, dominance genetic, and shared environmental effects
#' are taken into account.
#'
#' The three parameters should have a sum: h2+c2+d2 that is less than 1. The
#' total variance is set to 1, and the mean is zero.
#'
#' @param n numeric. The number of simulations to generate
#' @param pedigree a \code{pedigree} object
#' @param h2 numeric. The heritability
#' @param c2 numeric. The environmentability
#' @param d2 numeric. The dominance deviance effect
#' @return Returns a matrix with the simulated values with n columns (one for
#' each simulation) and each row matches the corresponding individual from the
#' pedigree
#' @author Claus Ekstrom \email{claus@@rprimer.dk}
#' @seealso \code{pedigree}, \code{kinship},
#' @keywords datagen
#' @examples
#'
#' library(kinship2)
#' library(mvtnorm)
#' mydata <- data.frame(id=1:5,
#'                      dadid=c(NA, NA, 1, 1, 1),
#'                      momid=c(NA, NA, 2, 2, 2),
#'                      sex=c("male", "female", "male", "male", "male"),
#'                      famid=c(1,1,1,1,1))
#' relation <- data.frame(id1=c(3), id2=c(4), famid=c(1), code=c(1))
#' ped <- pedigree(id=mydata$id, dadid=mydata$dadid, momid=mydata$momid,
#'                 sex=mydata$sex, relation=relation)
#' rmvtnorm.pedigree(2, ped, h2=.25)
#'
#' @export rmvtnorm.pedigree
rmvtnorm.pedigree <- function(n=1, pedigree, h2=0, c2=0, d2=0) {

# simulate.pedigree

  if (! "pedigree" %in% class(pedigree))
    stop("requires a single pedigree as input")

  if ((sum(h2, c2, d2)>=1) | min(h2, c2, d2) <0)
    stop("input parameters wrong. Must me non-negative and sum to less than 1")

  # Twins are automatically handled by the kinship function from kinship2 if coded correctly in the
  # pedigree

  if (d2!=0)
    stop("cannot handle d2 != 0 yet")

  # Could be programmed easily for non-inbred pedigrees using the fact that Delta7 = phi_km phi_ln + phi_kn phi_lm
  # (see Lange p. 80).

  pedsize <- length(pedigree$id)

  # VCmat per block
  vcmat <- (1-h2-c2)*diag(pedsize) + h2*2*kinship2::kinship(pedigree) + c2
  t(mvtnorm::rmvnorm(n, mean=rep(0, pedsize), sigma=vcmat))
}



#' Simulate residual multivariate t-distributed data from a polygenic model
#'
#' Simulates residual multivariate t-distributed response data from a pedigree
#' where the additive genetic, dominance genetic, and shared environmental
#' effects are taken into account.
#'
#' The three parameters should have a sum: h2+c2+d2 that is less than 1. The
#' total variance is set to 1, and the mean is zero.
#'
#' @param n numeric. The number of simulations to generate
#' @param pedigree a \code{pedigree} object
#' @param h2 numeric. The heritability
#' @param c2 numeric. The environmentability
#' @param d2 numeric. The dominance deviance effect
#' @param df numeric. The degrees of freedom for the t distribution
#' @return Returns a matrix with the simulated values with n columns (one for
#' each simulation) and each row matches the corresponding individual from the
#' pedigree
#' @author Claus Ekstrom \email{claus@@rprimer.dk}
#' @seealso \code{pedigree}, \code{kinship},
#' @keywords datagen
#' @examples
#'
#' library(kinship2)
#' library(mvtnorm)
#' mydata <- data.frame(id=1:5,
#'                      dadid=c(NA, NA, 1, 1, 1),
#'                      momid=c(NA, NA, 2, 2, 2),
#'                      sex=c("male", "female", "male", "male", "male"),
#'                      famid=c(1,1,1,1,1))
#' relation <- data.frame(id1=c(3), id2=c(4), famid=c(1), code=c(1))
#' ped <- pedigree(id=mydata$id, dadid=mydata$dadid, momid=mydata$momid,
#'                 sex=mydata$sex, relation=relation)
#' rmvt.pedigree(2, ped, h2=.25, df=4)
#'
#' @export rmvt.pedigree
rmvt.pedigree <- function(n=1, pedigree, h2=0, c2=0, d2=0, df=1) {

# simulate.pedigree

  if (! "pedigree" %in% class(pedigree))
    stop("requires a single pedigree as input")

  if ((sum(h2, c2, d2)>=1) | min(h2, c2, d2) <0)
    stop("input parameters wrong. Must me non-negative and sum to less than 1")

  # Twins are automatically handled by the kinship function from kinship2 if coded correctly in the
  # pedigree

  if (d2!=0)
    stop("cannot handle d2 != 0 yet")

  # Could be programmed easily for non-inbred pedigrees using the fact that Delta7 = phi_km phi_ln + phi_kn phi_lm
  # (see Lange p. 80).

  pedsize <- length(pedigree$id)

  # VCmat per block
  vcmat <- (1-h2-c2-d2)*diag(pedsize) + h2*2*kinship2::kinship(pedigree) + c2
  t(mvtnorm::rmvt(n, df=df, sigma=vcmat))
}


#rmvnorm.matlist <- function(n, mean=0,
#                            sigma, varlist, sigmaI=1) {

#  if (is.matrix(varlist))
#    varlist <- list(varlist)
#  if (!is.list(varlist))
#    stop("varlist must be a list")

#  if (length(sigma) != length(varlist))
#    stop("lengths of sigma and varlist must match")

  # Check that all matrices have thesame dimensions?
#  size <- nrow(varlist[[1]])

#  vcmat <- sigmaI*diag(size)

#  sapply(1:length(sigma), function(i) {vcmat <- vcmat + sigma[i]*varlist[[i]]} )

#  for (i in 1:length(sigma)) {
#    lapply(varlist)
#  }
#  as.vector(rmvnorm(1, mean=rep(0, pedsize), sigma=vcmat))


#}
