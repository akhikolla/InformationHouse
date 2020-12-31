####==========================================================================####
## The R code for fitting the disease models with K nearest-neighbor regression.##
##                                                                              ##
####==========================================================================####
#' @title K nearest-neighbor (KNN) regression
#'
#' @description  \code{rhoKNN} uses the KNN approach to estimate the probabilities of the disease status in case of three categories.
#'
#' @param X  a numeric design matrix.
#' @param Dvec   a n * 3  binary matrix with three columns, corresponding to the three classes of the disease status. In row i, 1 in column j indicates that the i-th subject belongs to class j, with j = 1, 2, 3. A row of \code{NA} values indicates a non-verified subject.
#' @param V  a binary vector containing the verification status (1 verified, 0 not verified).
#' @param K  an integer value/vector, which indicates the number of nearest neighbors. It should be less than the number of the verification subjects.
#' @param type  a distance measure.
#' @param trace switch for tracing estimation process. Default \code{FALSE}.
#'
#' @details \code{type} should be selected as one of \code{"eucli"}, \code{"manha"}, \code{"canber"}, \code{"lagran"}, \code{"mahala"} corresponding to Euclidean, Manhattan, Canberra, Lagrange and Mahalanobis distance. In practice, the selection of a suitable distance is typically dictated by features of the data and possible subjective evaluations. For example, if the covariates are heterogeneous with respect to their variances (which is particularly true when the variables are measured on heterogeneous scales), the choice of the Mahalanobis distance may be a good choice.
#'
#' For the number of nearest neighbors, a small value of \code{K}, within the range 1-3, may be a good choice. In general, the choice of \code{K} may depend on the dimension of the feature space, and propose to use cross--validation to find \code{K} in case of high--dimensional covariate. See \code{\link{CVknn}}.
#'
#'
#' @return \code{rhoKNN} returns a list containing the following components:
#'  \item{values}{estimates of the probabilities.}
#'  \item{X}{a design model matrix.}
#'  \item{K}{the number of nearest neighbors.}
#'  \item{type}{the chosen distance.}
#'
#' @references
#' To Duc, K., Chiogna, M., Adimari, G. (2016): Nonparametric Estimation of ROC Surfaces Under Verification Bias. \url{https://arxiv.org/abs/1604.04656v1}. Submitted.
#'
#' @examples
#' data(EOC)
#' XX <- cbind(EOC$CA125, EOC$CA153, EOC$Age)
#' Dna <- preDATA(EOC$D, EOC$CA125)
#' Dvec.na <- Dna$Dvec
#'
#' ## Euclidean distance, K = 1
#' out.ecul.1nn <- rhoKNN(XX, Dvec.na, EOC$V, K = 1, type = "eucli")
#'
#' ## Manhattan distance, K = 1
#' out.manh.1nn <- rhoKNN(XX, Dvec.na, EOC$V, K = 1, type = "manha")
#'
#' ## Canberra distance, K = 3
#' out.canb.1nn <- rhoKNN(XX, Dvec.na, EOC$V, K = 3, type = "canber")
#'
#' ## Lagrange distance, K = 3
#' out.lagr.1nn <- rhoKNN(XX, Dvec.na, EOC$V, K = 3, type = "lagran")
#'
#' ## Mahalanobis distance, K = c(1,3)
#' out.maha.13nn <- rhoKNN(XX, Dvec.na, EOC$V, K = c(1,3), type = "mahala")
#'
#' @export
rhoKNN <- function(X, Dvec, V, K, type = c("eucli", "manha", "canber", "lagran",
                                            "mahala"), trace = FALSE){
  if(!inherits(X, "matrix")) stop("\"XX\" not a matrix \n")
  if(!inherits(Dvec, "matrix") | ncol(Dvec) != 3 | !all(is.element(na.omit(Dvec), c(0,1)))) stop("variable \"Dvec\" must be a binary matrix with 3 columns")
  if(nrow(X) != nrow(Dvec)) stop(gettextf("arguments imply differing number of observation: %d", nrow(X)), gettextf(", %d", nrow(Dvec)), domain = NA)
  if(missing(V)) stop("object \"V\" is missing \n")
  if(missing(K)) stop("object \"K\" is missing \n")
  nk <- length(K)
  if(nk == 0) stop("what is the value of K!? \n")
  is.wholenumber <- function(x) all(x == round(x)) & all(x > 0)
  if(!is.wholenumber(K)) stop("\"K\" not a positive integer number or vector \n")
  n <- nrow(X)
  S.inv <- solve(cov(X))
  type <- match.arg(type)
  dista <- switch(type,
                  eucli = function(xi, xx, Sx.inv){
                    ss <- colSums((xi - t(xx))^2)
                    return(sqrt(ss))
                  },
                  manha = function(xi, xx, Sx.inv){
                    ss <- colSums(abs(xi - t(xx)))
                    return(ss)
                  },
                  canber = function(xi, xx, Sx.inv){
                    ss <- colSums( abs(xi - t(xx))/( abs(xi) + t(abs(xx)) ) )
                    return(ss)
                  },
                  lagran = function(xi, xx, Sx.inv){
                    ss <- apply(abs(xi - t(xx)), 2, max)
                    return(ss)
                  },
                  mahala = function(xi, xx, Sx.inv){
                    ss1 <- (xi - t(xx))
                    ss <- diag(t(ss1)%*%Sx.inv%*%ss1)
                    return(sqrt(ss))
                  }
  )
  if(trace){
    cat("Fitting the disease model by using K nearest-neighbor imputation.\n")
    cat("K:", K, "\n")
    cat("Distance measure:", switch(type,
                                    eucli = "Euclidean",
                                    manha = "Manhattan",
                                    canber = "Canberra",
                                    lagran = "Lagrange",
                                    mahala = "Mahalanobis"), "\n")
    cat("\n")
  }
  Dvec.veri <- Dvec[V == 1, ]
  X.veri <- X[V == 1, ]
  id.ms <- which(V == 0)
  id.veri <- which(V == 1)
  if(nk == 1){
    res <- matrix(0, nrow = n, ncol = 3)
    for(i in id.ms){
      dist.tem <- dista(X[i,], X.veri, S.inv)
      id.tem <- order(dist.tem)
      y.temp <- Dvec.veri[id.tem, ]
      res[i,] <- c(mean(y.temp[1:K,1]), mean(y.temp[1:K,2]), mean(y.temp[1:K,3]))
    }
    for(i in id.veri){
      dist.tem <- dista(X[i,], X.veri, S.inv)
      id.tem <- order(dist.tem)[-1]
      y.temp <- Dvec.veri[id.tem, ]
      res[i,] <- c(mean(y.temp[1:K,1]), mean(y.temp[1:K,2]), mean(y.temp[1:K,3]))
    }
    fit <- list(values = res, X = X, K = K, type = type)
    class(fit) <- "prob_dise_knn"
  }
  else{
    res <- list()
    for(j in 1:nk){
      res[[j]] <- matrix(0, nrow = n, ncol = 3)
    }
    for(i in id.ms){
      dist.tem <- dista(X[i,], X.veri, S.inv)
      id.tem <- order(dist.tem)
      y.temp <- Dvec.veri[id.tem, ]
      for(j in 1:nk){
        res[[j]][i,] <- c(mean(y.temp[1:K[j],1]), mean(y.temp[1:K[j],2]),
                          mean(y.temp[1:K[j],3]))
      }
    }
    for(i in id.veri){
      dist.tem <- dista(X[i,], X.veri, S.inv)
      id.tem <- order(dist.tem)[-1]
      y.temp <- Dvec.veri[id.tem, ]
      for(j in 1:nk){
        res[[j]][i,] <- c(mean(y.temp[1:K[j],1]), mean(y.temp[1:K[j],2]),
                          mean(y.temp[1:K[j],3]))
      }
    }
    fit <- list()
    for(j in 1:nk){
      fit[[j]] <- list(values = res[[j]], X = X, K = K[j], type = type)
      class(fit[[j]]) <- "prob_dise_knn"
    }
  }
  invisible(fit)
}
