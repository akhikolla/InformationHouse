#' Permutation-based variable importance metric for high dimensional datasets appropriate for time to event outcomes,
#' in the presence of imperfect self-reports or laboratory-based diagnostic tests.
#'
#' @description Let N and P denote the number of subjects and number of variables in the dataset, respectively.
#'              Let N** denote the total number of visits, summed over all subjects in the study
#'              [i.e. N** denotes the number of diagnostic test results available for all subjects in the study].
#'               This algorithm builds a user-defined number of survival trees, using bootstrapped datasets.
#'               Using the out of bag (OOB) data in each tree, a permutation-based measure of
#'               variable importance for each of the P variables is obtained.
#'
#' @param data name of the data frame that includes the variables subject, testtimes, result
#' @param subject vector of subject IDs of length N**x1.
#' @param testtimes vector of visit or test times of length N**x1.
#' @param result vector of binary diagnostic test results (0 = negative for event of interest; 1 = positive
#'        for event of interest) of length N**x1.
#' @param sensitivity the sensitivity of the diagnostic test.
#' @param specificity the specificity of the diagnostic test.
#' @param Xmat a N x P matrix of covariates.
#' @param root.size the minimum number of subjects in a terminal node.
#' @param ns number of covariate selected at each node to split the tree.
#' @param pval P-value threshold of the Likelihood Ratio Test.

#' @return a vector of the ensembled variable importance for modified random survival forest (icRSF).
#' @export
#' @examples
#' data(Xmat)
#' data(pheno)
#' tree <- treebuilder(data=pheno, subject=ID, testtimes=time, result=result, sensitivity=1,
#'                    specificity=1, Xmat=Xmat, root.size=30, ns=sqrt(ncol(Xmat)), pval=1)
#'
#' @useDynLib icRSF
#' @importFrom Rcpp evalCpp
#' @importFrom stats optim rbinom rexp
#' @import parallel icensmis
#'
#'






treebuilder <- function(data, subject, testtimes, result, sensitivity, specificity, Xmat, root.size, ns, pval=1){

  id <- eval(substitute(subject), data, parent.frame())
  time <- eval(substitute(testtimes), data, parent.frame())
  result <- eval(substitute(result), data, parent.frame())
  ord <- order(id, time)
  if (is.unsorted(ord)) {
    id <- id[ord]
    time <- time[ord]
    result <- result[ord]
    data <- data[ord, ]
  }
  utime <- sort(unique(time))

  stopifnot(is.numeric(sensitivity), is.numeric(specificity),
            is.numeric(root.size),
            is.numeric(pval), is.numeric(ns))
  stopifnot(length(sensitivity) == 1, sensitivity >= 0, sensitivity <= 1,
            length(specificity) == 1, specificity >= 0, specificity <= 1,
            length(root.size) == 1, root.size >= 1, root.size <= nrow(Xmat),

            length(ns) ==1, ns >=1, ns<=ncol(Xmat),
            all(time > 0), all(is.finite(time)))
  if (!all(result %in% c(0, 1))) stop("result must be 0 or 1")

  if (any(tapply(time, id, anyDuplicated)))
    stop("existing duplicated visit times for some subjects")
  #############################################################################
  # Compute D matrix
  #############################################################################

  timen0 <- (time != 0)
  Dmat <- dmat(id[timen0], time[timen0], result[timen0], sensitivity,
               specificity, 1)
  row.names(Dmat) <- unique(id)

  hdidx <- 1
  tree <- 1
  thres <- stats::qchisq(1-pval, df=1)

  getchisq <- function(Dmat, x) {
    # J <- dim(Dmat)[2]-1
    # parmi <- c(0, log(-log((J:1)/(J+1))))
    # parmD <- c(parmi[1], diff(c(0, parmi[-1])))
    # parmD0 <- parmD[-1]

    ## log-likelihood function
    if(length(unique(x))==1){
      chisq <- 0
      parm <- rep(0, J+1)
    }else{
      q <- try(optim(parmD, loglikCD, gradlikCD, lower = c(rep(-Inf, 2), rep(1e-8, J-1)), Dmat = Dmat, x = x, method="L-BFGS-B"))
      q0 <- try(optim(parmD0, loglikCD0, gradlikCD0, lower = c(-Inf, rep(1e-8, J-1)), Dmat = Dmat, method="L-BFGS-B"))
      if (class(q)=="try-error" | class(q0)=="try-error") return(rep(0, J+2))
      chisq <- 2*(q0$value - q$value)
      parm <- q$par
    }
    return(c(chisq, parm))
  }
  #############################################################################
  # Compute D matrix
  #############################################################################
  #
  #   timen0 <- (time != 0)
  #   Dm <- icensmis:::dmat(id[timen0], time[timen0], result[timen0], sensitivity,
  #                         specificity, 1)
  simdata <- cbind(Dmat, Xmat)
  J <- ncol(Dmat) - 1
  nsub <- nrow(Dmat)
  ncov <- ncol(Xmat)

  parmi <- c(0, log(-log((J:1)/(J+1))))
  parmD <- c(parmi[1], diff(c(0, parmi[-1])))
  parmD0 <- parmD[-1]

  bootid <- sample(1:nrow(simdata),nrow(simdata), replace = T)
  bootdata <- simdata[bootid, ]
  OOBid <- setdiff(1:nsub, unique(bootid))
  OOB <- simdata[setdiff(1:nsub, unique(bootid)), ]
  nOOB <- nrow(OOB)
  #  OOB <- simdata[setdiff(1:nsub, unique(bootid)), ]
  #  nOOB <- nrow(OOB)
  ################# TREE BUILDING ####################

  newtree <- function(firstval, firstc, inc, parm, firstnobs) {
    m <- matrix(rep(NA, inc*(J+6)), nrow = inc, ncol = J + 6)
    m[1, 3] <- firstval
    m[1, 4] <- firstc
    m[1, 5:(J+5)] <- parm
    m[1, J+6] <- firstnobs
    return(list(mat=m, nxt=2, inc=inc))
  }

  ## function to insert value to the tree
  insert <- function(hdidx, dir, tr, newval, newc, parm, nobs) {
    newidx <- tr$nxt
    if (tr$nxt == nrow(tr$mat) + 1) {
      tr$mat <- rbind(tr$mat, matrix(rep(NA, tr$inc*(J+6)), nrow=tr$inc, ncol=J+6))
    }
    tr$mat[newidx, 3] <- newval
    tr$mat[newidx, 4] <- newc
    tr$mat[newidx, 5:(J+5)] <- parm
    tr$mat[newidx, (J+6)] <- nobs

    tr$mat[hdidx, dir] <- newidx
    tr$nxt <- tr$nxt + 1
    return(tr)
  }

  training <- function(bootdata, hdidx, dir, tree, root.size, ns, pval) {
    ## stopping rule
    if (length(unique(row.names(bootdata))) > root.size) {
      Dmat <- bootdata[, 1:(J+1)]
      Xmat <- bootdata[, -(1:(J+1))]
      #selid <- 1:10
      # set.seed(1)
      selid <- sample(1:ncov, ns, replace=FALSE)
      #ind <- apply(Xmat[,selid], 2, function(x) length(unique(x))>1)
      #selid <- selid[ind]
      X <- Xmat[, selid]
      bsplit <- apply(X, 2, function(t) splitpointC(Dmat, t, getchisq))
      if(length(selid) > 0){
        id <- which.max(bsplit[2, ])
        # print(max(bsplit[2, ]))
        if(max(bsplit[2, ])>= thres){
          bvar <- selid[id]
          bcp <- bsplit[1, id] ##change here if change to continuous
          bootdata.L <- bootdata[X[,id] <= bcp, ]
          bootdata.R <- bootdata[X[, id] > bcp, ]
          parm <- bsplit[-(1:2), id]
          nobs <- nrow(bootdata)
          if (nrow(bootdata) == nsub) {
            tree <- newtree(bvar, bcp, 3, parm, nobs)
          } else {
            tree <- insert(hdidx, dir, tree, bvar, bcp, parm, nobs)
          }
          a <- tree$nxt - 1
          tree <- training(bootdata.L, a, 1, tree, root.size, ns, pval)
          tree <- training(bootdata.R, a, 2, tree, root.size, ns, pval)
        }
      }
    }
    return(tree)
  }

  tr <- training(bootdata, hdidx, dir, tree, root.size, ns, pval)
  tr<- tr$mat
  tr <- tr[!is.na(tr[, 3]), ]
  return(list(tree=tr, OOBid=OOBid))
}


