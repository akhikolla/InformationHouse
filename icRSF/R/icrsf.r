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
#' @param root.size minimum number of subjects in a terminal node.
#' @param ntree number of survival trees.
#' @param ns number of covariate selected at each node to split the tree.
#' @param node For parallel computation, specify the number of nodes.
#' @param pval P-value threshold of the Likelihood Ratio Test.

#' @return a vector of the ensembled variable importance for modified random survival forest (icRSF).
#' @export
#' @examples
#' library(parallel)
#' data(Xmat)
#' data(pheno)
#' vimp <- icrsf(data=pheno, subject=ID, testtimes=time, result=result, sensitivity=1,
#'              specificity=1, Xmat=Xmat, root.size=30, ntree=1, ns=sqrt(ncol(Xmat)), node=1, pval=1)
#'
#' @useDynLib icRSF
#' @importFrom Rcpp evalCpp
#' @importFrom stats optim rbinom rexp
#' @import parallel icensmis
#'
#'
icrsf <- function(data, subject, testtimes, result, sensitivity, specificity, Xmat,
                  root.size, ntree, ns, node, pval=1){

  # get parameter function
  getparm <- function(id, Xmat, treemat, hdidx = 1) {
    stopifnot(is.numeric(id), length(id) == 1)
    vec <- treemat[hdidx, ]
    obs <- Xmat[id, ]
    dir <- (obs[vec[3]] >= vec[4]) + 1
    if (is.na(vec[dir])) {
      fparm <- vec[-(1:5)] ### changes here
      fparm[1] <- fparm[1] + (dir-1)*vec[5]*obs[vec[3]] ### changes here
      return(list(fparm=fparm, id=hdidx)) ### changes here
    } else {
      getparm(id, Xmat, treemat, vec[dir])
    }
  }

  # calculate varimp function
  varimp1 <- function(Dmat, Xmat, treemat, OOBid){
    J <- ncol(Dmat) - 1

    nOOB <- length(OOBid)
    ncov <- ncol(Xmat)

    newpred <- function(inc){
      m <- rep(NA, inc)
      m[1] <- 1
      #    m[1, -1] <- parm
      return(list(mat=m, nxt=2, inc=inc))
    }
    insertpred <- function(pred, newval){
      newidx <- pred$nxt
      if (pred$nxt == length(pred) + 1) {
        pred$mat <- c(pred$mat, rep(NA, pred$inc))
      }
      pred$mat[newidx] <- newval
      #  pred$mat[newidx, -1] <- parm
      pred$nxt <- pred$nxt + 1
      return(pred)
    }
    getparm <- function(x, treemat, hdidx = 1) {
      vec <- treemat[hdidx, ]
      dir <- (x[vec[3]] >= vec[4]) + 1
      if (is.na(vec[dir])) {
        fparm <- vec[6:(J+5)]
        fparm[1] <- fparm[1] + (dir-1)*vec[5]*x[vec[3]]
        return(list(fparm=fparm, id=hdidx)) ### Also store the id of
      } else {
        getparm(x, treemat, vec[dir])
      }
    }
    getpred <- function(x, treemat, hdidx = 1, pmat) {
      vec <- treemat[hdidx, ]
      dir <- (x[vec[3]] >= vec[4]) + 1
      if(is.na(hdidx) == FALSE){
        if(hdidx == 1){
          pmat <- newpred(3)
        } else
        {
          pmat <- insertpred(pmat, hdidx)
        }
        hdidx <- vec[dir]
        pmat <- getpred(x, treemat, hdidx, pmat)
      }
      return(pmat)
    }
    #############################################################################
    # Retreive terminal parameters and compute log likelihood
    #############################################################################
    OOBloglik <- function(Dmat, Xmat, OOBid, varid, treemat){
      simdata <- cbind(Dmat, Xmat)
      OOB <- simdata[OOBid, ]
      OOBDmat <- OOB[, 1:(J+1)]
      OOBDes <- OOB[, -(1:(J+1))]
      nOOB <- nrow(OOB)
      if (varid > 0) OOBDes[, varid] <- sample(OOBDes[, varid], nOOB, replace=FALSE)
      p <- apply(OOBDes, 1, function(x) getparm(x, treemat))
      parms <- t(sapply(p, function(x) cumsum(x$fparm)))
      # LIKELIHOOD
      l <- rowSums(OOBDmat*cbind(1, exp(-exp(parms))))
      freq <- unlist(sapply(1:nOOB, function(x) getpred(OOBDes[x, ], treemat, 1, 1)[[1]]))
      freq <- freq[!is.na(freq)]
      l1 <- ifelse(l>0, log(l), NA)
      return(l1)
    }
    cov <- unique(treemat[, 3])
    cov <- cov[!is.na(cov)]
    varimp <- rep(0, ncov)
    lik.org <- OOBloglik(Dmat, Xmat, OOBid, 0, treemat)

    ## This is to get the loglik after permuation of each of the selected covariates in the tree
    permutelik <- lapply(cov, function(x) OOBloglik(Dmat, Xmat, OOBid, x, treemat))
    varimp[cov] <- sapply(permutelik, function(x) sum(lik.org - x, na.rm=T))
    return(varimp)
  }


  #############################################################################
  # Function to build trees
  #############################################################################
  # treebuilder(Dmat, Xmat, 10, 10, 0.05)

  treebuilder <- function(Dmat, Xmat, root.size, ns, pval=1){
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

  ## Main body
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

  #############################################################################
  # Check variables and input parameters check
  #############################################################################
  stopifnot(is.numeric(sensitivity), is.numeric(specificity),
            is.numeric(root.size),  is.numeric(node),
            is.numeric(pval), is.numeric(ns))
  stopifnot(length(sensitivity) == 1, sensitivity >= 0, sensitivity <= 1,
            length(specificity) == 1, specificity >= 0, specificity <= 1,
            length(root.size) == 1, root.size >= 1, root.size <= nrow(Xmat),
            length(node) ==1, node >=1,
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

  # Dmat <- icensmis:::dmat(id[timen0], time[timen0], result[timen0], sensitivity,
  #                        specificity, 1)

  trees <- mclapply(1:ntree, function(x) treebuilder(Dmat, Xmat, root.size, ns, pval),
                    mc.cores=node)
  treemats <- mclapply(trees, function(x) x[[1]], mc.cores=node)
  OOBids <- mclapply(trees, function(x) x[[2]], mc.cores=node)

  vimp <- mclapply(1:ntree, function(x) varimp1(Dmat, Xmat, treemats[[x]], OOBid=OOBids[[x]]))

  vimp.adj <- mclapply(vimp, function(x) ifelse(x>0, x, 0))
  vimp1 <- do.call("rbind", vimp.adj)
  vimp2 <- colSums(vimp1, na.rm=T)
  return(ensemble.vimp=vimp2)
}

# vimp <- icrsf(data=pheno, subject=ID, testtimes=time, result=result, sensitivity=1,
#            specificity=1, Xmat=Xmat, root.size=5, ntree=1, ns=sqrt(ncol(Xmat)), node=1, pval=0.1)


