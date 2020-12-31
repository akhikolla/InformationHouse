#' @title Fit pk searcher efficiency models.
#'
#' @description Searcher efficiency is modeled as a function of the number of
#'  times a carcass has been missed in previous searches and any number of
#'  covariates. Format and usage parallel that of common \code{R} functions
#'  \code{lm}, \code{glm}, and \code{gam}. However, the input data
#'  (\code{data}) is structured differently to accommodate the
#'  multiple-search searcher efficiency trials (see Details), and model
#'  formulas may be entered for both \code{p} (akin to an intercept) and
#'  \code{k} (akin to a slope).
#'
#' @details
#'   The probability of finding a carcass that is present at the time of
#'   search is \code{p} on the first search after carcass arrival and is
#'   assumed to decrease by a factor of \code{k} each time the carcass is
#'   missed in searches. Both \code{p} and \code{k} may depend on covariates
#'   such as ground cover, season, species, etc., and a separate model format
#'   (\code{formula_p} and \code{formula_k}) may be entered for each. The
#'   models are entered as they would be in the familiar \code{lm} or
#'   \code{glm} functions in R. For example, \code{p} might vary with
#'   \code{A} and \code{B}, while \code{k} varies
#'   only with \code{A}. A user might then enter \code{p ~ A + B}
#'   for \code{formula_p} and \code{k ~ A} for
#'   \code{formula_k}. Other R conventions for defining formulas may also be
#'   used, with \code{A:B} for the interaction between covariates
#'   A and B and \code{A * B} as short-hand for \code{A + B + A:B}.
#'
#'   Search trial \code{data} must be entered in a data frame with data in
#'   each row giving the fate of a single carcass in the field trials. There
#'   must be a column for each search occassion, with 0, 1, or NA depending on
#'   whether the carcass was missed, found, or not available (typically
#'   because it was found and removed on a previous search, had been earlier
#'   removed by  scavengers, or was not searched for) on the given search
#'   occasion. Additional columns with values for categorical covariates
#'   (e.g., visibility = E, M, or D) may also be included.
#'
#'  When all trial carcasses are either found on the first search or
#'  are missed on the first search after carcass placement, pkm effects a
#'  necessary adjustment to the for accuracy; otherwise, the model would not be
#'  able to determine the uncertainty and would substantially over-estimate the
#'  variance of the parameter estimates, giving \eqn{\hat{p}} essentially equal
#'  to 0 or 1 with approximately equal probability. The adjustment is to fit the
#'  model on an adjusted data set with duplicated copies of the original data
#'  (\code{2n} observations) but with one carcass having the opposite fate of the
#'  others. For example, in field trials with very high searcher efficiency and
#'  \code{n = 10} carcasses, all of which are found in the first search after
#'  carcass placement, the original data set would have a carcass observation
#'  column consisting of 1s (\code{rep(1, 10)}). The adjusted data set would
#'  have an observation column consisting of \code{2n - 1} 1s and one 0. In this
#'  case, the point estimate of \code{p} is \code{1/(2n)} with distribution that
#'  closely resembling the Bayesian posterior distributions of \code{p} with a
#'  uniform or a Jeffreys prior. The adjustment is applied on a cellwise basis
#'  in full cell models (e.g., 1, A, B, A * B). In the additive model with two
#'  predictors (A + B), the adjustment is made only when a full level of
#'  covariate A or B is all 0s or 1s.
#'
#' @param formula_p Formula for p; an object of class \code{\link{formula}}
#'   (or one that can be coerced to that class): a symbolic description of
#'   the model to be fitted. Details of model specification are given under
#'   "Details".
#'
#' @param formula_k Formula for k; an object of class \code{\link{formula}}
#'   (or one that can be coerced to that class): a symbolic description of the
#'   model to be fitted. Details of model specification are given under
#'   "Details".
#'
#' @param data Data frame with results from searcher efficiency trials and any
#' covariates included in \code{formula_p} or \code{formula_k} (required).
#'
#' @param obsCol Vector of names of columns in \code{data} where results
#'  for each search occasion are stored (optional). If \code{obsCol} is not
#'  provided, \code{pkm} uses as \code{obsCol} all columns with names that
#'  begin with an \code{"s"} or \code{"S"} and end with a number, e.g., "s1",
#'  "s2", "s3", etc. This option is included as a convenience for the user,
#'  but care must be taken that other data are not stored in columns with
#'  names matching that pattern. Alternatively, \code{obsCol} may be
#'  entered as a vector of names, like \code{c("s1", "s2", "s3")},
#'  \code{paste0("s", 1:3)}, or \code{c("initialSearch", "anotherSearch",
#'  "lastSearch")}. The columns must be in chronological order, that is, it is
#'  assumed that the first column is for the first search after carcass arrival,
#'  the second column is for the second search, etc.
#'
#' @param kFixed Parameter for user-specified \code{k} value (optional). If a
#'   value is provided, \code{formula_k} is ignored and the model is fit under
#'   the assumption that the \code{k} parameter is fixed and known to be
#'   \code{kFixed} \eqn{\in [0, 1]}. If a \code{sizeCol} is provided, \code{kFixed}
#'   may either be \code{NULL}, a single number in [0, 1], or a vector with
#'   \code{kFixed} values for two or more of the carcass size classes. For
#'   example, if there are three sizes (\code{S}, \code{M}, and \code{L}),
#'   \code{kFixed} could be \code{c(S = 0.3, M = 0.8, L = 1.0)} to assign fixed
#'   \code{k} values to each size. To fit \code{k} for size \code{S} and to assign
#'   values of 0.8 and 1.0 to sizes \code{M} and \code{L}, resp., use
#'   \code{kFixed = c(S = 0.3, M = 0.8, L = 1.0)}. If there are more than one size
#'   classes and \code{kFixed} is a scalar, then all size classes are assigned the
#'   same \code{kFixed} value (unless \code{kFixed} is named, e.g.,
#'   \code{kFixed = c(S = 0.5)}, in which case only the named size is assigned the
#'   \code{kFixed}).
#'
#' @param allCombos logical. If \code{allCombos = FALSE}, then the single model
#'  expressed by \code{formula_p} and \code{formula_k} is fit using a call to
#'  \code{pkm0}. If \code{allCombos = TRUE}, a full set of \code{\link{pkm}}
#'  submodels derived from combinations of the given covariates for \code{p}
#'  and \code{k} is fit. For example, submodels of \code{formula_p = p ~ A * B}
#'  would be \code{p ~ A * B}, \code{p ~ A + B}, \code{p ~ A}, \code{p ~ B},
#'  and \code{p ~ 1}. Models for each pairing of a \code{p} submodel with a
#' \code{k} submodel are fit via \code{pkmSet}, which fits each model
#'  combination using successive calls to \code{pkm0}, which fits a
#'  single model.
#'
#' @param sizeCol character string. The name of the column in \code{data} that
#'  gives the carcass class of the carcasses in the field trials. If
#'  \code{sizeCol = NULL}, then models are not segregated by size. If a
#'  \code{sizeCol} is provided, then separate models are fit for the \code{data}
#'  subsetted by \code{sizeCol}.
#'
#' @param CL numeric value in (0, 1). confidence level
#'
#' @param kInit numeric value in (0, 1). Initial value used for numerical
#'  optimization of \code{k}. Default is \code{kInit = 0.7}. It is rarely
#' (if ever) necessary to use an alternative initial value.
#'
#' @param quiet Logical indicator of whether or not to print messsages
#'
#' @param ... additional arguments passed to subfunctions
#'
#' @return an object of an object of class \code{pkm}, \code{pkmSet},
#'  \code{pkmSize}, or \code{pkmSetSize}.
#' \describe{
#'  \item{\code{pkm0()}}{returns a \code{pkm} object, which is a description
#'    of a single, fitted pk model. Due to the large number and complexity of
#'    components of a\code{pkm} model, only a subset of them is printed
#'    automatically; the rest can be viewed/accessed via the \code{$} operator
#'    if desired. These are described in detail in the '\code{pkm} Components'
#'    section.}
#'  \item{\code{pkmSet()}}{returns a list of \code{pkm} objects, one for each
#'    of the submodels, as described with parameter \code{allCombos = TRUE}.}
#'  \item{\code{pkmSize()}}{returns a list of \code{pkmSet} objects (one for
#'    each 'size') if \code{allCombos = T}, or a list of \code{pkm} objects (one
#'    for each 'size') if \code{allCombos = T}}
#'  \item{\code{pkm}}{returns a \code{pkm}, \code{pkmSet}, \code{pkmSize}, or
#'    \code{pkmSetSize} object:
#'     \itemize{
#'        \item \code{pkm} object if \code{allCombos = FALSE, sizeCol = NULL}
#'        \item \code{pkmSet} object if \code{allCombos = TRUE, sizeCol = NULL}
#'        \item \code{pkmSize} object if \code{allCombos = FALSE, sizeCol != NULL}
#'        \item \code{pkmSetSize} object if \code{allCombos = TRUE, sizeCol != NULL}
#'     }
#'  }
#' }
#' @section \code{pkm} Components:
#'
#' The following components of a \code{pkm} object are displayed automatically:
#'
#' \describe{
#'  \item{\code{call}}{the function call to fit the model}
#'  \item{\code{formula_p}}{the model formula for the \code{p} parameter}
#'   \item{\code{formula_k}}{the model formula for the \code{k} parameter}
#'  \item{\code{predictors}}{list of covariates of \code{p} and/or \code{k}}
#'  \item{\code{AICc}}{the AIC value as corrected for small sample size}
#'  \item{\code{convergence}}{convergence status of the numerical optimization
#'   to find the maximum likelihood estimates of \code{p} and \code{k}. A
#'    value of \code{0} indicates that the model was fit successfully. For
#'    help in deciphering other values, see \code{\link{optim}}.}
#'  \item{\code{cell_pk}}{summary statistics for estimated cellwise estimates
#'    of \code{p} and \code{k}, including the number of carcasses in each cell,
#'    medians and upper & lower bounds on CIs for each parameter, indexed by
#'    cell (or combination of covariate levels).}
#' }
#'
#' The following components are not printed automatically but can be accessed
#' via the \code{$} operator:
#' \describe{
#'   \item{\code{data}}{the data used to fit the model}
#'   \item{\code{data0}}{\code{$data} with NA rows removed}
#'   \item{\code{betahat_p, betahat_k}}{parameter estimates for the terms in the
#'    regression model for for \code{p} and \code{k} (logit scale). If \code{k}
#'    is fixed or not provided, then \code{betahat_k} is not calculated.}
#'   \item{\code{varbeta}}{the variance-covariance matrix of the estimators
#'     for \code{c(betahat_p, betahat_k)}.}
#'  \item{\code{cellMM_p, cellMM_k}}{cellwise model (design) matrices for
#'    covariate structures of \code{p_formula} and \code{k_formula}}
#'   \item{\code{levels_p, levels_k}}{all levels of each covariate of \code{p}
#'     and \code{k}}
#'   \item{\code{nbeta_p, nbeta_k}}{number of parameters to fit the \code{p}
#'     and \code{k} models}
#'   \item{\code{cells}}{cell structure of the pk-model, i.e., combinations of
#'     all levels for each covariate of \code{p} and \code{k}. For example, if
#'     \code{covar1} has levels \code{"a"}, \code{"b"}, and \code{"c"}, and
#'     \code{covar2} has levels \code{"X"} and \code{"Y"}, then the cells
#'     would consist of \code{a.X}, \code{a.Y}, \code{b.X}, \code{b.Y},
#'     \code{c.X}, and \code{c.Y}.}
#'   \item{\code{ncell}}{total number of cells}
#'  \item{\code{predictors_k, predictors_p}}{covariates of \code{p} and \code{k}}
#'  \item{\code{observations}}{observations used to fit the model}
#'  \item{\code{kFixed}}{the input \code{kFixed}}
#'  \item{\code{AIC}}{the
#'    \href{https://en.wikipedia.org/wiki/Akaike_information_criterion}{AIC}
#'    value for the fitted model}
#'  \item{\code{carcCells}}{the cell to which each carcass belongs}
#'  \item{\code{CL}}{the input \code{CL}}
#'  \item{\code{loglik}}{the log-liklihood for the maximum likelihood estimate}
#'  \item{\code{pOnly}}{a logical value telling whether \code{k} is included in
#'    the model. \code{pOnly = TRUE} if and only if \code{length(obsCol) == 1)}
#'    and \code{kFixed = NULL}}.
#'  \item{\code{data_adj}}{\code{data0} as adjusted for the 2n fix to accommodate
#'    scenarios in which all trial carcasses are either found or all are not
#'    found on the first search occasion (uncommon)}
#'  \item{\code{fixBadCells}}{vector giving the names of cells adjusted for the
#'    2n fix}
#' }
#'
#' @section Advanced:
#'  \code{pkmSize} may also be used to fit a single model for each carcass class if
#'  \code{allCombos = FALSE}. To do so, \code{formula_p} and \code{formula_k}
#'  must be a named list of formulas with names matching the sizes listed in
#'  \code{unique(data[, sizeCol])}. The return value is then a list of
#'  \code{pkm} objects, one for each size.
#'
#' @seealso \code{\link{rpk}}, \code{\link{qpk}}, \code{\link{aicc}},
#'  \code{\link{plot.pkm}}
#'
#' @examples
#'  head(data(wind_RP))
#'  mod1 <- pkm(formula_p = p ~ Season, formula_k = k ~ 1, data = wind_RP$SE)
#'  class(mod1)
#'  mod2 <- pkm(formula_p = p ~ Season, formula_k = k ~ 1, data = wind_RP$SE,
#'    allCombos = TRUE)
#'  class(mod2)
#'  names(mod2)
#'  class(mod2[[1]])
#'  mod3 <- pkm(formula_p = p ~ Season, formula_k = k ~ 1, data = wind_RP$SE,
#'    allCombos = TRUE, sizeCol = "Size")
#'  class(mod3)
#'  names(mod3)
#'  class(mod3[[1]])
#'  class(mod3[[1]][[1]])
#'
#' @export
#'
pkm <- function(formula_p, formula_k = NULL, data, obsCol = NULL, kFixed = NULL,
    allCombos = FALSE, sizeCol = NULL, CL = 0.90, kInit = 0.7, quiet = FALSE,
    ...){
  if (!is.null(kFixed) && !is.numeric(kFixed))
    stop("kFixed must be NULL or numeric")
  if (any(na.omit(kFixed) < 0 | na.omit(kFixed) > 1)){
    badk <- names(which(na.omit(kFixed) < 0 | na.omit(kFixed) > 1))
    stop("invalid k for ", paste0(badk, collapse = ", "))
  }
  if (is.null(allCombos) || is.na(allCombos) || !is.logical(allCombos)){
    stop("allCombos must be TRUE or FALSE")
  }
  if (is.null(sizeCol) || is.na(sizeCol)){
    if (!allCombos){ # single model
      out <- pkm0(formula_p = formula_p, formula_k = formula_k, data = data,
        obsCol = obsCol, kFixed = kFixed, CL = CL, kInit = kInit, quiet = quiet)
    } else { # allCombos of p and k subformulas
      out <- pkmSet(formula_p = formula_p, formula_k = formula_k, data = data,
        obsCol = obsCol, kFixed = kFixed, kInit = kInit, CL = CL, quiet = quiet)
    }
  } else { # specified formula for p and k, split by carcass class
    out <- pkmSize(formula_p = formula_p, formula_k = formula_k, data = data,
      obsCol = obsCol, kFixed = kFixed, sizeCol = sizeCol, allCombos = allCombos,
      CL = CL, kInit = kInit, quiet = quiet)
   }
  return(out)
}

#' @rdname pkm
#' @export
pkm0 <- function(formula_p, formula_k = NULL, data, obsCol = NULL,
    kFixed = NULL, kInit = 0.7, CL = 0.90, quiet = FALSE){
  i <- sapply(data, is.factor)
  data[i] <- lapply(data[i], as.character)
  if (!is.null(kFixed) && is.na(kFixed)) kFixed <- NULL
  if (!is.null(kFixed) && !is.numeric(kFixed))
    stop("kFixed must be NULL or numeric")
  if (!is.null(kFixed) && (kFixed < 0 | kFixed > 1)){
    stop("invalid fixed value for k")
  }
  if(any(! obsCol %in% colnames(data))){
    stop("Observation column provided not in data.")
  }
  if (length(obsCol) == 0){
    obsCol <- grep("^[sS].*[0-9]$", names(data), value = TRUE)
    nobsCol <- length(obsCol)
    if (nobsCol == 0){
      stop("No obsCol provided and no appropriate column names found.")
    }
  }
  predCheck <- c(all.vars(formula_p[[3]]), all.vars(formula_k[[3]]))
  if (any(!(predCheck %in% colnames(data)))){
    stop("User-supplied formula includes predictor(s) not found in data.")
  }
  if (length(kFixed) >= 1){
    if (!is.numeric(kFixed[1])){
      stop("User-supplied kFixed must be numeric (or NULL)")
    }
    if (length(kFixed) > 1){
      kFixed <- kFixed[1]
      if (!quiet){
        message("Vector-valued kFixed. Only the first element will be used.")
      }
    }
    if (kFixed < 0 || kFixed > 1){
      stop("User-supplied kFixed is outside the supported range [0, 1].")
    }
    if (length(formula_k) > 0 & quiet == FALSE){
      message("Formula and fixed value provided for k, fixed value used.")
      formula_k <- NULL
    }
  }
  pOnly <- FALSE
  if (is.null(kFixed)){
    if (length(obsCol) == 1){
      if (!is.null(formula_k) && is.language(formula_k) && quiet == FALSE){
        message("Only one search occasion per carcass. k not estimated.")
      }
      pOnly <- TRUE
      formula_k <- NULL
      kFixed <- 1
    } else {
      # flag to indicate no estimation of k
      if(is.null(formula_k) || !is.language(formula_k)){
        pOnly <- TRUE
        obsCol <- obsCol[1] # use data from first search only
        message("p is estimated from data in first search occasion only.")
        formula_k <- NULL
        kFixed <- 1
      }
    }
  }
  nsearch <- length(obsCol)
  if (ncol(data) == nsearch) data$id <- 1:nrow(data)
  obsData <- as.matrix(data[ , obsCol], ncol = nsearch)

  # replace all non-zero/non-one data with NA:
  if (!is.numeric(obsData)){
    obsData[!is.na(obsData) & !(obsData %in% as.character(0:1))] <- NA
    obsData <- matrix(as.numeric(obsData), ncol = nsearch)
  } else {
    obsData[!(obsData %in% 0:1)] <- NA
  }

  # remove rows that are all NAs
  onlyNA <- (rowSums(is.na(obsData)) == nsearch)
  obsData <- as.matrix(obsData[!onlyNA, ], ncol = nsearch)
  data00 <- data[!onlyNA, ]
  data0 <- data00 # may be modified later to adjust for bad cells
  ncarc <- nrow(obsData)

  preds_p <- all.vars(formula_p[[3]])
  formulaRHS_p <- formula(delete.response(terms(formula_p)))
  levels_p <- .getXlevels(terms(formulaRHS_p), data0)

  preds_k <- character(0)
  if (is.language(formula_k)){
    preds_k <- all.vars(formula_k[[3]])
    formulaRHS_k <- formula(delete.response(terms(formula_k)))
    levels_k <- .getXlevels(terms(formulaRHS_k), data0)
  }
  if (length(kFixed) == 1){
    preds_k <- character(0)
    formulaRHS_k <- formula(~1)
    formula_k <- c(fixedk = kFixed)
    levels_k <- .getXlevels(terms(formulaRHS_k), data0)
  }

  preds <- unique(c(preds_p, preds_k))
  for (pri in preds)
    if (is.numeric(data0[ , pri])) data0[ , pri] <- paste0("_", data0[ , pri])
  cells <- combinePreds(preds, data0)
  ncell <- nrow(cells)
  cellNames <- cells$CellNames
  cellMM_p <- model.matrix(formulaRHS_p, cells)
  cellMM_k <- model.matrix(formulaRHS_k, cells)
  cellMM <- cbind(cellMM_p, cellMM_k)

  if (length(preds) == 0){
    carcCells0 <- rep("all", ncarc)
  } else if (length(preds) == 1){
    carcCells0 <- data0[ , preds]
  } else if (length(preds) > 1){
    carcCells0 <- do.call(paste, c(data0[ , preds], sep = '.'))
  }
  pInitCellMean <- tapply(data0[ , obsCol[1]], INDEX = carcCells0, FUN = mean, na.rm = TRUE)
### prepSE.r stops at this point
  fixBadCells <- NULL
  if (NCOL(cellMM_p) == prod(unlist(lapply(levels_p, length))) & # full cell model
      any(pInitCellMean %in% 0:1)){# ...with bad cells
    # employ the 2n fix for bad cells:
    if (length(preds_p) == 0){ # no predictors
      if (all(data0[ , obsCol[1]] == 1) || all(data0[ , obsCol[1]] == 0)){
        data0 <- rbind(data0, data0) # use 2n
        data0[1, obsCol[1]] <- 1 - data0[1, obsCol[1]]
        if (length(obsCol) > 1 && data0[1, obsCol[1]] == 1){
          data0[1, obsCol[-1]]  <- NA
        }
        fixBadCells <- "all"
      }
    } else {
      fixBadCells <- names(pInitCellMean)[pInitCellMean %in% 0:1]
      badCells <- cells[cells$CellNames %in% fixBadCells, ]
      badCells$CellNames <- NULL
      for (ci in 1:NROW(badCells)){
        cellind <- which(matrixStats::colProds( # factor levels match cell
          t(data0[ , colnames(badCells)]) ==
          as.character(badCells[ci, ])) == 1)
        data0 <- rbind(data0[cellind, ], data0) # the 2n fix
        data0[1, obsCol[1]] <- 1 - data0[1, obsCol[1]]
        if (data0[1, obsCol[1]] == 1 & length(obsCol) > 1){
          data0[1, obsCol[-1]] <- NA
        }
      }
    }
  } else if (length(preds_p) == 2 &  # two predictors but not full cell
        NCOL(cellMM_p) < prod(unlist(lapply(levels_p, length)))){# "+" model
    for (predi in names(levels_p)){
      for (li in levels_p[[predi]]){
        cellind  <- which(data0[ , predi] == li)
        if (all(data0[cellind, obsCol[1]] == 0) |
            all(data0[cellind, obsCol[1]] == 1))
          stop("Initial search has all 0s or 1s for ", preds_p[1], " = ", li,
             " in additive model.")
      }
    }
  }

  obsData <- as.matrix(data0[ , obsCol])
  colnames(obsData) <- obsCol
  misses <- matrixStats::rowCounts(obsData, value = 0, na.rm =T)
  foundOn <- matrixStats::rowMaxs(
      obsData * matrixStats::rowCumsums(1 * !is.na(obsData)), na.rm = T)

  dataMM_p <- model.matrix(formulaRHS_p, data0)
  dataMM_k <- model.matrix(formulaRHS_k, data0)
  dataMM <- t(cbind(dataMM_p, dataMM_k))

  nbeta_k <- ncol(dataMM_k)
  nbeta_p <- ncol(dataMM_p)
  nbeta <- nbeta_p + nbeta_k
  if (length(preds) == 0){
    carcCells <- rep("all", dim(data0)[1])
  } else if (length(preds) == 1){
    carcCells <- data0[ , preds]
  } else if (length(preds) > 1){
    carcCells <- do.call(paste, c(data0[ , preds], sep = '.'))
  }
  cellByCarc <- match(carcCells, cellNames)
  pInitCellMean <- tapply(data0[ , obsCol[1]], INDEX = carcCells, FUN = mean, na.rm = TRUE)

  pInit <- as.vector(pInitCellMean[match(carcCells, names(pInitCellMean))])
  pInit[which(pInit < 0.1)] <- 0.1
  pInit[which(pInit > 0.9)] <- 0.9

  cellMatrix_p <- solve(t(dataMM_p) %*% dataMM_p)
  cellImpact_p <- t(dataMM_p) %*% logit(pInit)
  betaInit_p <- cellMatrix_p %*% cellImpact_p
  betaInit_k <- logit(rep(kInit, nbeta_k))
  betaInit <- c(betaInit_p, betaInit_k)

  if (length(kFixed) == 1){
    betaInit <- betaInit[-length(betaInit)]
  }
  for (ki in 1:3){ # three tries at different initial values for k to correct if var < 0
    MLE <- tryCatch(
             optim(par = betaInit, fn = pkLogLik,
               hessian = TRUE, cellByCarc = cellByCarc, misses = misses,
               maxmisses = max(misses), foundOn = foundOn, cellMM = cellMM,
               nbeta_p = nbeta_p, kFixed = kFixed, method = "BFGS"),
             error = function(x) {NA}
           )
    if (length(MLE) == 1 && is.na(MLE)) stop("Failed optimization.")
    varbeta <- tryCatch(solve(MLE$hessian), error = function(x) {NA})
    if (is.na(varbeta)[1]) stop("Unable to estimate variance.")
    if (sum(diag(varbeta) < 0) > 0) {
      if (ki == 1) {
        betaInit_k <- logit(rep(0.05, nbeta_k))
        betaInit <- c(betaInit_p, betaInit_k)
      } else if (ki == 2) {
        betaInit_k <- logit(rep(0.95, nbeta_k))
        betaInit <- c(betaInit_p, betaInit_k)
      } else {
        stop("Unable to estimate k")
      }
    } else {
      break
    }
  }
  convergence <- MLE$convergence
  betahat <- MLE$par
  if (is.null(fixBadCells)){
    llik <- -MLE$value
  } else {
    obsData <- as.matrix(data00[ , obsCol])
    colnames(obsData) <- obsCol
    misses <- matrixStats::rowCounts(obsData, value = 0, na.rm =T)
    foundOn <- matrixStats::rowMaxs(
        obsData * matrixStats::rowCumsums(1 * !is.na(obsData)), na.rm = T)

    if (length(preds) == 0){
      carcCells <- rep("all", dim(data00)[1])
    } else if (length(preds) == 1){
      carcCells <- data00[ , preds]
    } else if (length(preds) > 1){
      carcCells <- do.call(paste, c(data00[ , preds], sep = '.'))
    }
    cellByCarc <- match(carcCells, cellNames)
    llik <- -pkLogLik(misses = misses, foundOn = foundOn, beta = betahat,
      nbeta_p = nbeta_p, cellByCarc = cellByCarc, maxmisses = max(misses),
      cellMM = cellMM, kFixed = kFixed)
    ncarc <- nrow(obsData)
  }
  nparam <- length(betahat)
  AIC <- 2*nparam - 2*llik
  AICcOffset <- (2 * nparam * (nparam + 1)) / (nrow(obsData) - nparam - 1)
  AICc <- round(AIC + AICcOffset, 2)
  AIC <- round(AIC, 2)

  betahat_p <- betahat[1:nbeta_p]
  names(betahat_p) <- colnames(dataMM_p)
  betahat_k <- NULL
  if (length(kFixed) == 0){
    betahat_k <- betahat[(nbeta_p + 1):(nbeta)]
    names(betahat_k) <- colnames(dataMM_k)
  }

  varbeta_p <- varbeta[1:nbeta_p, 1:nbeta_p]
  cellMean_p <- cellMM_p %*% betahat_p
  cellVar_p <- cellMM_p %*% varbeta_p %*% t(cellMM_p)
  cellSD_p <- suppressWarnings(sqrt(diag(cellVar_p)))

  if (is.null(kFixed) || is.na(kFixed)){
    which_k <- (nbeta_p + 1):(nbeta)
    varbeta_k <- varbeta[which_k, which_k]
    cellMean_k <- cellMM_k %*% betahat_k
    cellVar_k <- cellMM_k %*% varbeta_k %*% t(cellMM_k)
    cellSD_k <- suppressWarnings(sqrt(diag(cellVar_k)))
  } else {
    cellMean_k <- rep(logit(kFixed), ncell)
    cellSD_k <- rep(0, ncell)
  }

  probs <- list(0.5, (1 - CL) / 2, 1 - (1 - CL) / 2)
  cellTable_p <- lapply(probs, qnorm, mean = cellMean_p, sd = cellSD_p)
  cellTable_p <- matrix(unlist(cellTable_p), ncol = 3)
  cellTable_p <- round(alogit(cellTable_p), 6)
  colnames(cellTable_p) <- c("p_median", "p_lwr", "p_upr")
  cell_n <- as.numeric(table(carcCells0)[cellNames])
  names(cell_n) <- NULL
  if (!pOnly){
    cellTable_k <- lapply(probs, qnorm, mean = cellMean_k, sd = cellSD_k)
    cellTable_k <- matrix(unlist(cellTable_k), nrow = ncell, ncol = 3)
    cellTable_k <- round(alogit(cellTable_k), 6)
    colnames(cellTable_k) <- c("k_median", "k_lwr", "k_upr")
    cellTable <- data.frame(
      cell = cellNames,
      n = cell_n,
      cellTable_p,
      cellTable_k)
  } else {
    cellTable <- data.frame(cell = cellNames, n = cell_n, cellTable_p)
  }

  output <- list()
  output$call <- match.call()
  output$data <- data
  output$data0 <- data00
  output$formula_p <- formula_p
  if (!pOnly) output$formula_k <- formula_k
  output$predictors <- preds
  output$predictors_p <- preds_p
  if (!pOnly) output$predictors_k <- preds_k
  output$AIC <- AIC
  output$AICc <- AICc
  output$convergence <- convergence
  output$varbeta <- varbeta
  output$cellMM_p <- cellMM_p
  if (!pOnly) output$cellMM_k <- cellMM_k
  output$nbeta_p <- nbeta_p
  if (!pOnly) output$nbeta_k <- nbeta_k
  output$betahat_p <- betahat_p
  if (!pOnly) output$betahat_k <- betahat_k
  output$levels_p <- levels_p
  if (!pOnly) output$levels_k <- levels_k
  output$cells <- cells
  output$ncell <- ncell
  output$cell_pk <- cellTable
  output$CL <- CL
  output$observations <- obsData
  if (!pOnly) output$kFixed <- kFixed
  output$carcCells <- carcCells
  output$loglik <- llik
  output$pOnly <- pOnly
  if (is.null(fixBadCells)) data_adj <- NULL
  output$data_adj <- data0
  output$fixBadCells <- fixBadCells
  class(output) <- c("pkm", "list")
  attr(output, "hidden") <- c("data", "data0", "predictors_p", "predictors_k",
    "kFixed", "betahat_p", "betahat_k", "cellMM_p", "cellMM_k", "nbeta_p",
    "nbeta_k", "varbeta", "levels_p", "levels_k", "carcCells", "AIC", "cells",
    "ncell", "observations", "loglik", "pOnly", "data_adj", "fixBadCells")
  return(output)
}
#' @title Print a \code{\link{pkm}} model object
#'
#' @description Print a \code{\link{pkm}} model object
#'
#' @param x a \code{\link{pkm}} model object
#'
#' @param ... to be passed down
#'
#' @export
#'
print.pkm <- function(x, ...){
  hid <- attr(x, "hidden")
  notHid <- !names(x) %in% hid
  print(x[notHid])
}
 
#' @title Calculate the negative log-likelihood of a searcher efficiency model
#' 
#' @description The function used to calculate the negative-loglikelihood of
#'   a given searcher efficiency model (\code{\link{pkm}}) with a given data
#'   set
#'
#' @param misses Number of searches when carcass was present but
#'  not found.
#'
#' @param foundOn Search on which carcass was found.
#'
#' @param beta Parameters to be optimized.
#'
#' @param nbeta_p Number of parameters associated with p.
#'
#' @param cellByCarc Which cell each observation belongs to.
#'
#' @param maxmisses Maximum possible number of misses for a carcass.
#'
#' @param cellMM Combined pk model matrix.
#'
#' @param kFixed Value of k if fixed. 
#'
#' @return Negative log likelihood of the observations, given the parameters.
#'
#' @export 
#'
pkLogLik <- function(misses, foundOn, beta, nbeta_p, cellByCarc, maxmisses, 
                     cellMM, kFixed = NULL){
  if (!is.null(kFixed) && is.na(kFixed)) kFixed <- NULL
  if (length(kFixed) == 1){
    beta <- c(beta, logit(kFixed))
  }

  ncell <- nrow(cellMM)
  nbeta <- length(beta)
  which_p <- 1:nbeta_p
  which_k <- (nbeta_p + 1):nbeta

  beta_p <- beta[which_p]
  beta_k <- beta[which_k]
  Beta <- matrix(0, nrow = nbeta, ncol = 2)
  Beta[which_p, 1] <- beta[which_p]
  Beta[which_k, 2] <- beta[which_k]

  pk <- alogit(cellMM %*% Beta)
  p <- pk[ , 1]
  k <- pk[ , 2]

  powk <- matrix(k, nrow = ncell, ncol = maxmisses + 1)
  powk[ , 1] <- 1
  powk <- matrixStats::rowCumprods(powk)

  pmiss <- matrix(1 - (p * powk[ , 1:(maxmisses + 1)]), nrow = ncell)
  pmiss <- matrixStats::rowCumprods(pmiss)
  pfind <- matrixStats::rowDiffs(1 - pmiss)
  pfind_si <- cbind(pk[ , 1], pfind)

  notFoundCell <- cellByCarc[foundOn == 0]
  notFoundMisses <- misses[foundOn == 0]
  notFoundCellMisses <- cbind(notFoundCell, notFoundMisses)
  foundCell <- cellByCarc[foundOn > 0]
  foundFoundOn <- foundOn[foundOn > 0]
  foundCellFoundOn <- cbind(foundCell, foundFoundOn)

  ll_miss <- sum(log(pmiss[notFoundCellMisses]))
  ll_found <- sum(log(pfind_si[foundCellFoundOn]))
  nll_total <- -(ll_miss + ll_found)
 
  return(nll_total)
}

#' @rdname pkm
#' @export
#'
pkmSet <- function(formula_p, formula_k = NULL, data, obsCol = NULL, 
                   kFixed = NULL, kInit = 0.7, CL = 0.90, quiet = FALSE){
  if (!is.null(kFixed) && is.na(kFixed)) kFixed <- NULL
  if (length(kFixed) > 0){
    if (length(kFixed) > 1){
      warning(
        "More than one fixed kFixed value provided by user. ",
        "Only the first element will be used."
      )
      kFixed <- kFixed[1]
    }
    if (is.numeric(kFixed) && !is.na(kFixed) && length(formula_k) > 0){ 
      if(quiet == FALSE){
        message("Formula and fixed value provided for k, fixed value used.")
      }
      formula_k <- NULL
    }
  }
  if (sum(obsCol %in% colnames(data)) == 1 & length(formula_k) > 0){
    if (quiet == FALSE){
      message("Only one search occasion for each carcass; k not estimated.")
    }
    formula_k <- NULL
  }

  unfixk <- FALSE
  if (length(formula_k) == 0){
    if (length(kFixed) == 0){
      kFixed <- 0.5
      unfixk <- TRUE
    }
  }

  # create the set of models to explore, based on the input parameters
  terms_p <- attr(terms(formula_p), "term.labels")
  if (length(formula_k) == 0){
    terms_k <- NULL
  } else {
    terms_k <- attr(terms(formula_k), "term.labels")
  }
  nterms_p <- length(terms_p)
  nterms_k <- length(terms_k)
  nformula_p <- 2^(nterms_p)
  nformula_k <- 2^(nterms_k)

  dropComplex_p <- rep(1:nterms_p, choose(nterms_p, 1:nterms_p))
  dropWhich_p <- numeric(0)
  if (nterms_p > 0){
    for (termi in 1:nterms_p){
      specificDrop <- seq(1, choose(nterms_p, (1:nterms_p)[termi]))
      dropWhich_p <- c(dropWhich_p, specificDrop)
    }
  }
  optionFormula_p <- vector("list", nformula_p)
  optionFormula_p[[1]] <- formula_p
  keepFormula_p <- rep(TRUE, nformula_p)
  if (nformula_p > 1){
    for (formi in 2:nformula_p){
      termDropComplex <- combn(terms_p, dropComplex_p[formi - 1])
      termDropSpec <- termDropComplex[ , dropWhich_p[formi - 1]]
      termDrop <- paste(termDropSpec, collapse = " - ")
      formulaUpdate <- paste(format(~.), "-", termDrop)
      updatedFormula <- update.formula(formula_p, formulaUpdate)
      optionFormula_p[[formi]] <- updatedFormula
      keepFormula_p[formi] <- checkComponents(updatedFormula)
    }
    nkeepFormula_p <- sum(keepFormula_p)
    whichKeepFormula_p <- which(keepFormula_p == TRUE)
    keptFormula_p <- vector("list", nkeepFormula_p)
    for (kepti in 1:nkeepFormula_p){
      keptFormula_p[[kepti]] <- optionFormula_p[[whichKeepFormula_p[kepti]]]
    }
  } else {
    keptFormula_p <- optionFormula_p
  }
  
  dropComplex_k <- rep(1:nterms_k, choose(nterms_k, 1:nterms_k))
  dropWhich_k <- numeric(0)
  if (nterms_k > 0){
    for (termi in 1:nterms_k){
      specificDrop <- seq(1, choose(nterms_k, (1:nterms_k)[termi]))
      dropWhich_k <- c(dropWhich_k, specificDrop)
    }
  }
  optionFormula_k <- vector("list", nformula_k)
  optionFormula_k[[1]] <- formula_k
  keepFormula_k <- rep(TRUE, nformula_k)
  if (nformula_k > 1){
    for (formi in 2:nformula_k){
      termDropComplex <- combn(terms_k, dropComplex_k[formi - 1])
      termDropSpec <- termDropComplex[ , dropWhich_k[formi - 1]]
      termDrop <- paste(termDropSpec, collapse = " - ")
      formulaUpdate <- paste(format(~.), "-", termDrop)
      updatedFormula <- update.formula(formula_k, formulaUpdate)
      optionFormula_k[[formi]] <- updatedFormula
      keepFormula_k[formi] <- checkComponents(updatedFormula)
    }
    nkeepFormula_k <- sum(keepFormula_k)
    whichKeepFormula_k <- which(keepFormula_k == TRUE)
    keptFormula_k <- vector("list", nkeepFormula_k)
    for (kepti in 1:nkeepFormula_k){
      keptFormula_k[[kepti]] <- optionFormula_k[[whichKeepFormula_k[kepti]]]
    }
  } else {
    keptFormula_k <- optionFormula_k
  }
  if (length(kFixed) == 1){
    keptFormula_k <- NA
  }

  expandedKeptFormulae <- expand.grid(keptFormula_p, keptFormula_k)
  keptFormula_p <- expandedKeptFormulae[ , 1]
  keptFormula_k <- expandedKeptFormulae[ , 2]
  if (length(kFixed) == 1){
    keptFormula_k <- NULL
  }
  if (unfixk == TRUE){
    kFixed <- NULL
  }
  nmod <- nrow(expandedKeptFormulae) 
  output <- vector("list", nmod)
  for (modi in 1:nmod){
    formi_p <- keptFormula_p[modi][[1]]
    formi_k <- keptFormula_k[modi][[1]]
    pkm_i <- tryCatch(
               pkm(formula_p = formi_p, formula_k = formi_k, data = data,
                 obsCol = obsCol, kFixed = kFixed, CL = CL, kInit = kInit,
                 quiet = quiet),
               error = function(x) {
                       paste("Failed model fit: ", geterrmessage(), sep = "")
               }
             )
    name_p <- paste(format(formi_p), collapse = "")
    name_p <- gsub("    ", "", name_p)
    name_k <- paste(format(formi_k), collapse = "")
    name_k <- gsub("    ", "", name_k)
    if (length(kFixed) == 1){
      name_k <- paste("k fixed at ", kFixed, sep = "")
    }
    modName <- paste(name_p, "; ", name_k, sep = "")
    modName <- gsub("NULL", "k not estimated", modName)
    output[[modi]] <- pkm_i
    names(output)[modi] <- modName
  }
  class(output) <- c("pkmSet", "list")
  return(output)
} 

#' @rdname pkm
#' @export
#'
pkmSize <- function(formula_p, formula_k = NULL, data, kFixed = NULL,
    obsCol = NULL, sizeCol = NULL, allCombos = FALSE, kInit = 0.7,
    CL = 0.90, quiet = FALSE){

  if (length(sizeCol) == 0 || is.na(sizeCol)){
    pkfunc <- ifelse(allCombos, "pkmSet", "pkm0")
    out <- list()
    out[["all"]] <- get(pkfunc)(formula_p = formula_p, formula_k = formula_k,
      data = data, obsCol = obsCol, kFixed = kFixed, kInit = kInit,
      CL = CL, quiet = quiet
    )
    class(out) <- c(ifelse(allCombos, "pkmSetSize", "pkmSize"), "list")
    return(out)
  }

  if (!(sizeCol %in% colnames(data))){
    stop("sizeCol not in data set.")
  }

  sizeclasses <- sort(unique(as.character(data[ , sizeCol])))

  if (all(is.na(kFixed))){
    kFixed <- NULL
  }
  if (length(kFixed) >= 1){
    kFixed <- kFixed[!is.na(kFixed)]
    if (length(kFixed) == 1 && is.null(names(kFixed))){
    # if exactly one kFixed is provided and no name is given,
    # all sizes take on the same kFixed value
      kFixed <- rep(kFixed, length(sizeclasses))
      names(kFixed) <- sizeclasses
      if (quiet == FALSE){
        message("One unnamed kFixed value provided by user. ",
          "All classes are assumed to have kFixed = ", kFixed[1], ". ",
          "To specify specific classes to apply kFixed values to, ",
          "class names must be provided in kFixed vector."
        )
      }
    } else {
      if (is.null(names(kFixed)) || !all(names(kFixed) %in% sizeclasses)){
        stop("kFixed names must be names of carcass classes.")
      }
    }
  }
  if (!allCombos){
    if ("list" %in% intersect(class(formula_p), class(formula_k))){
      # then fit the specific models for each formula and corresponding size
      if (!setequal(names(formula_p), names(formula_k)) ||
          !setequal(names(formula_p), unique(data[ , sizeCol]))){
        stop("p and k formula names must match carcass classes")
      }
      formlist <- TRUE
    } else {
      formlist <- FALSE
    }
  }
  out <- list()
  for (sci in sizeclasses){
    if (allCombos){
      out[[sci]] <- pkmSet(formula_p = formula_p, formula_k = formula_k,
        data = data[data[, sizeCol] == sci, ], obsCol = obsCol,
        kFixed = kFixed[sci], CL = CL, kInit = kInit, quiet = quiet)
    } else {
      if (formlist){
        out[[sci]] <- pkm0(
          formula_p = formula_p[[sci]], formula_k = formula_k[[sci]],
          data = data[data[, sizeCol] == sci, ], obsCol = obsCol,
          kFixed = kFixed[sci], CL = CL, kInit = kInit, quiet = quiet
        )
      } else {
        out[[sci]] <- pkm0(
          formula_p = formula_p, formula_k = formula_k,
          data = data[data[, sizeCol] == sci, ], obsCol = obsCol,
          kFixed = kFixed[sci], CL = CL, kInit = kInit, quiet = quiet
        )
      }
    }
  }
  class(out) <- c(ifelse(allCombos, "pkmSetSize", "pkmSize"), "list")
  return(out)
} 

#' @title Create the AICc tables for a set of searcher efficiency models
#' 
#' @description Generates model comparison tables based on AICc values for
#'   a set of pk models generated by \code{\link{pkmSet}}
#' 
#' @param x Set of searcher efficiency models fit to the same
#'   observations
#' 
#' @param ... further arguments passed to or from other methods
#'
#' @param quiet Logical indicating if messages should be printed
#' 
#' @param app Logical indicating if the table should have the app model names
#'
#' @return AICc table
#' 
#' @examples
#'   data(wind_RP)
#'   mod <- pkmSet(formula_p = p ~ Season, formula_k = k ~ Season, data = wind_RP$SE)
#'  aicc(mod)
#'
#' @export 
#'
aicc.pkmSet <- function(x, ... , quiet = FALSE, app = FALSE){
  pkmset <- x
  nmod <- length(pkmset)
  formulas <- names(pkmset)
  formulas_p <- rep(NA, nmod)
  formulas_k <- rep(NA, nmod)
  AICc <- rep(NA, nmod)
  deltaAICc <- rep(NA, nmod)

  if (nmod == 1){
    splitFormulas <- strsplit(formulas, "; ")[[1]]
    formulas_p <- splitFormulas[1]
    formulas_k <- splitFormulas[2]
    AICc <- tryCatch(pkmset[[1]]$AICc, error = function(x) {1e7})
    deltaAICc <- 0
    AICcOrder <- 1
  } else {
    for (modi in 1:nmod){
      splitFormulas_i <- strsplit(formulas[modi], "; ")[[1]]
      formulas_p[modi] <- splitFormulas_i[1]
      formulas_k[modi] <- splitFormulas_i[2]
      AICc[modi] <- tryCatch(pkmset[[modi]]$AICc, error = function(x) {1e7})
    }
    AICcOrder <- order(AICc)
    deltaAICc <- round(AICc - min(AICc), 2)
    which_fails <- which(AICc == 1e7)
    AICc[which_fails] <- NA
    deltaAICc[which_fails] <- NA
  }

  if (app){
    formulas_p <- gsub("~ 1", "~ constant", formulas_p)
    formulas_k <- gsub("~ 1", "~ constant", formulas_k)
  }

  output <- data.frame(formulas_p, formulas_k, AICc, deltaAICc)
  output <- output[AICcOrder, ]
  colnames(output) <- c("p Formula", "k Formula", "AICc", "\u0394AICc")
  whichAICcNA <- which(is.na(output$AICc))
  whichAICcMax <- which(output$AICc == 1e7)
  if (length(whichAICcNA) > 0 & quiet == FALSE){
    message("Models with incorrect specification were removed from output.")
    output <- output[-whichAICcNA, ]
  }
  if (length(whichAICcMax) > 0 & quiet == FALSE){
    message("Models that failed during fit were removed from output.")
    output <- output[-whichAICcMax, ]
  }
  class(output) <- c("corpus_frame", "data.frame")
  return(output)  # AIC
}

#' @title extract AICc value from pkm object
#'
#' @description extract AICc value from pkm object
#'
#' @param x object of class \code{pkm}
#'
#' @param ... further arguments passed to or from other methods
#'
#' @return Data frame with the formulas for p and k and the AICc of the model
#'
#' @export
#'
aicc.pkm <- function(x,...){
  return(
    data.frame(cbind(
      "formula_p" = deparse(x$formula_p),
      "formula_k" = deparse(x$formula_k),
      "AICc" = x$AICc
    ))
  )
}

#' @title Create the AICc tables for a list of sets of searcher efficiency models
#'
#' @description Generates model comparison tables based on AICc values for
#'   a set of pk models generated by \code{\link{pkm}} with
#'   \code{allCombos = TRUE} and a non-\code{NULL} \code{sizeCol}.
#'
#' @param x List of set of searcher efficiency models fit to the same
#'   observations
#'
#' @param ... further arguments passed to or from other methods
#'
#' @return AICc table
#'
#' @examples
#'   data(wind_RP)
#'   mod <- pkmSet(formula_p = p ~ Season, formula_k = k ~ Season, data = wind_RP$SE)
#'  aicc(mod)
#'
#' @export
#'
aicc.pkmSetSize <- function(x, ... ){
  return(lapply(x, FUN = function(y){
    class(y) <- c("pkmSet", "list")
    aicc(y)
  }))
}

#' @title Create the AICc tables for a list of sets of searcher efficiency models
#'
#' @description Generates model comparison tables based on AICc values for
#'   a set of pk models generated by \code{\link{pkm}} with
#'   \code{allCombos = FALSE} and a non-\code{NULL} \code{sizeCol}.
#'
#' @param x List of set of searcher efficiency models fit to the same
#'   observations
#'
#' @param ... further arguments passed to or from other methods
#'
#' @return AICc table
#'
#' @examples
#'   data(wind_RP)
#'   mod <- pkmSet(formula_p = p ~ Season, formula_k = k ~ Season, data = wind_RP$SE)
#'  aicc(mod)
#'
#' @export
#'
aicc.pkmSize <- function(x, ... ){
  return(lapply(x, FUN = function(y){
    class(y) <- c("pkm", "list")
    aicc(y)
  }))
}

#' @title Simulate parameters from a fitted pk model
#'
#' @description Simulate parameters from a \code{\link{pkm}} model object
#'
#' @param n the number of simulation draws
#'
#' @param model A \code{\link{pkm}} object (which is returned from 
#'   \code{pkm()})
#'
#' @return list of pairs of matrices of \code{n} simulated \code{p} and
#'  \code{k} for cells defined by the \code{model} object.
#'
#' @seealso \code{\link{rpk}}, \code{\link{pkm}}
#'
#' @examples
#'   data(wind_RP)
#'   mod <- pkm(formula_p = p ~ 1, formula_k = k ~ Season, data = wind_RP$SE)
#'   rpk(n = 10, model = mod)
#'
#' @export
#'
rpk <- function(n, model){
  if (!"pkm" %in% class(model)) stop("model not of class pkm")
  if (anyNA(model$varbeta) || sum(diag(model$varbeta) < 0) > 0){
    stop("Variance in pkm not well-defined. Cannot simulate.")
  }
  if (model$pOnly){
   stop("k not included in 'model'. Cannot simulate pk.")
  } else {
    which_beta_k <- (model$nbeta_p + 1):(model$nbeta_p + model$nbeta_k)
  }

  sim_beta <- mvtnorm::rmvnorm(n,
    mean = c(model$betahat_p, model$betahat_k),
    sigma = model$varbeta,
    method =  "svd")
  sim_p <- as.matrix(alogit(sim_beta[ , 1:model$nbeta_p] %*% t(model$cellMM_p)))
  colnames(sim_p) <- model$cells$CellNames

  if (is.null(model$kFixed) || is.na(model$kFixed)){
    sim_k <- as.matrix(alogit(sim_beta[ , which_beta_k] %*% t(model$cellMM_k)))
  } else {
    sim_k <- matrix(model$kFixed, ncol = model$ncell, nrow = n)
  }
  colnames(sim_k) <- model$cells$CellNames
  output <- lapply(model$cells$CellNames, function(x) cbind(p = sim_p[, x], k = sim_k[, x]))
  names(output) <- model$cells$CellNames

  return(output)
}

#' @title Quantiles of marginal distributions of \eqn{\hat{p}} and \eqn{\hat{k}}
#'
#' @description Calculate quantiles of marginal distributions of \eqn{\hat{p}}
#'  and \eqn{\hat{k}} for a \code{\link{pkm}} model object
#'
#' @param p vector of probabilities
#'
#' @param model A \code{\link{pkm}} object (which is returned from
#'   \code{pkm()})
#'
#' @return either a list of \code{ncell} \eqn{\times} \code{length(p)} matrices
#'  of quantiles for \code{$p} and \code{$k} for cells defined by the
#'  \code{model} object (if \code{model$pOnly == FALSE}) or a \code{ncell}
#'  \eqn{\times} \code{length(p)} matrix of quantiles for \code{p}
#'
#' @seealso \code{\link{rpk}}, \code{\link{pkm}}
#'
#' @examples
#'  # 90% confidence intervals for \code{p} and \code{k}
#'   mod <- pkm(formula_p = p ~ Visibility * Season, formula_k = k ~ Season,
#'    data = wind_cleared$SE)
#'   qpk(p = c(0.05, 0.95), model = mod)
#'
#' @export
#'
qpk <- function(p, model){
  if (!"pkm" %in% class(model)) stop("model not of class pkm")
  if (!is.numeric(p) || !is.vector(p)) stop("p must be a numeric vector")
  if (any(is.na(p))) stop("p must be numeric with no NA")
  if (max(p) >= 1 | min(p) <= 0) stop("p must be in (0, 1)")
  qp <- with(model, {
    varbeta_p <- varbeta[1:nbeta_p, 1:nbeta_p]
    cellMean_p <- cellMM_p %*% betahat_p
    cellVar_p <- cellMM_p %*% varbeta_p %*% t(cellMM_p)
    cellSD_p <- suppressWarnings(sqrt(diag(cellVar_p)))
    probs <- list(0.5, (1 - CL) / 2, 1 - (1 - CL) / 2)
    lapply(lapply(p, qnorm, mean = cellMean_p, sd = cellSD_p), alogit)
  })
  if (!model$pOnly){
    qk <- with(model, {
      if (!exists("kFixed") || is.null(kFixed) || is.na(kFixed)){
        which_k <- (nbeta_p + 1):(nbeta_p + nbeta_k)
        varbeta_k <- varbeta[which_k, which_k]
        cellMean_k <- cellMM_k %*% betahat_k
        cellVar_k <- cellMM_k %*% varbeta_k %*% t(cellMM_k)
        cellSD_k <- suppressWarnings(sqrt(diag(cellVar_k)))
      } else {
        cellMean_k <- rep(logit(kFixed), ncell)
        cellSD_k <- rep(0, ncell)
      }
      lapply(lapply(p, qnorm, mean = cellMean_k, sd = cellSD_k), alogit)
    })
    out <- list(
      p = matrix(unlist(qp), nrow = model$ncell),
      k = matrix(unlist(qk), nrow = model$ncell))
    rownames(out[["p"]]) <- rownames(out[["k"]]) <- model$cells$CellNames
    colnames(out[["p"]]) <- colnames(out[["k"]]) <- paste0("q", p)
    return(out)
  } else {
    qp <- matrix(unlist(qp), nrow = model$ncell)
    rownames(qp) <- model$cells$CellNames
    colnames(qp) <- paste0("q", p)
    return(qp)
  }
}

#' @title Check if a pk model is well-fit
#'
#' @description Run a check the arg is a well-fit pkm object
#'
#' @param pkmod A \code{\link{pkm}} object to test
#'
#' @return logical value indicating a failed fit (TRUE) or successful (FALSE)
#'
#' @export
#'
pkmFail <- function(pkmod){
!("pkm" %in% class(pkmod)) || anyNA(pkmod) || sum(diag(pkmod$varbeta) < 0) > 0
}


#' @title Check if pkm models fail
#' 
#' @description Run a check on each model within a \code{\link{pkmSet}} 
#'   object to determine if it failed or not
#'
#' @param pkmSetToCheck A \code{\link{pkmSet}} object to test
#'
#' @return A vector of logical values indicating if each of the models failed
#'
#' @export
#'
pkmSetFail <- function(pkmSetToCheck){
  unlist(lapply(pkmSetToCheck, pkmFail)) 
}

#' @title Check if all of the pkm models fail
#'
#' @description Run a check on each model within a \code{\link[=pkm]{pkmSetSize}}
#'   object to determine if they all failed or not
#'
#' @param pkmSetSizeToCheck A \code{pkmSetSize} object to test
#'
#' @return A list of logical vectors indicating which models failed
#'
#' @export
#'
pkmSetSizeFail <- function(pkmSetSizeToCheck){
  lapply(pkmSetSizeToCheck, pkmSetFail)
}

#' @title Check if all of the pkm models fail within a given set
#'
#' @description Run a check on each model within a \code{\link{pkmSet}}
#'   object to determine if they all failed or not
#'
#' @param pkmSetToCheck A \code{\link{pkmSet}} object to test
#'
#' @return A logical value indicating if all models failed in the set
#'
#' @export
#'
pkmSetAllFail <- function(pkmSetToCheck){
  modchecks <- unlist(lapply(pkmSetToCheck, pkmFail))
  if (length(modchecks) == sum(modchecks)){
    return(TRUE)
  }
  FALSE
}

#' @title Remove failed pkm models from a \code{pkmSet} object
#'
#' @description Remove all failed models within a \code{\link[=pkm]{pkmSet}} object
#'
#' @param pkmSetToTidy A \code{\link{pkmSet}} object to tidy
#'
#' @return A \code{\link{pkmSet}} object with failed models removed
#'
#' @export
#'
pkmSetFailRemove <- function(pkmSetToTidy){
  out <- pkmSetToTidy[!pkmSetFail(pkmSetToTidy)]
  class(out) <- c("pkmSet", "list")
  return(out)
}

#' @title Remove failed pkm models from a \code{pkmSetSize} object
#'
#' @description Remove failed models from a \code{\link[=pkm]{pkmSetSize}} object
#'
#' @param pkmSetSizeToTidy A list of \code{pkmSetSize} objects to tidy
#'
#' @return A list of \code{pkmSet} objects with failed models removed
#'
#' @export
#'
pkmSetSizeFailRemove <- function(pkmSetSizeToTidy){
  out <- list()
  for (sci in names(pkmSetSizeToTidy)){
    out[[sci]] <- pkmSetFailRemove(pkmSetSizeToTidy[[sci]])
  }
  class(out) <- c("pkmSetSize", "list")
  return(out)
}

#' @title Calculate decayed searcher efficiency
#'
#' @description Calculate searcher efficiency after some searches under 
#'   pk values
#'
#' @param days search days
#'
#' @param pk \code{p} and \code{k} values
#'
#' @return searcher efficiency that matches the output of ppersist
#'
#' @export 
#'
SEsi <- function(days, pk){ 
  if (is.null(dim(pk)) || nrow(pk) == 1) return (SEsi0(days, pk))
  npk <- nrow(pk)
  nsearch <- length(days) - 1
  ind1 <- rep(1:nsearch, times = nsearch:1)
  ind2 <- ind1 + 1
  ind3 <- unlist(lapply(1:nsearch, function(x) x:nsearch)) + 1
  schedule <- cbind(days[ind1], days[ind2], days[ind3])
  schedule.index <- cbind(ind1, ind2, ind3)
  nmiss <- schedule.index[, 3] - schedule.index[, 2]
  maxmiss <- max(nmiss)
  if (maxmiss == 0) {
      pfind.si <- pk[, 1]
  } else if (maxmiss == 1) {
      pfind.si <- cbind(pk[, 1], (1 - pk[, 1]) * pk[, 2] * pk[, 1])
  } else {
      powk <- array(rep(pk[, 2], maxmiss + 1), dim = c(npk, maxmiss + 1))
      powk[ , 1] <- 1
      powk <- matrixStats::rowCumprods(powk)
      pfind.si <- pk[, 1] * powk * cbind(
        rep(1, npk), matrixStats::rowCumprods(1 - (pk[, 1] * powk[, 1:maxmiss]))
      )
  }
  return(t(pfind.si)) 
}

#' @title Calculate decayed searcher efficiency for a single pk
#'
#' @description Calculate searcher efficiency after some searches for a 
#'   single pk combination
#'
#' @param days search days
#'
#' @param pk pk combination
#'
#' @return searcher efficiency that matches the output of ppersist
#'
#' @export 
#'
SEsi0 <- function(days, pk){ 
  nsearch <- length(days) - 1
  ind1 <- rep(1:nsearch, times = nsearch:1)
  ind2 <- ind1 + 1
  ind3 <- unlist(lapply(1:nsearch, function(x) x:nsearch)) + 1
  schedule <- cbind(days[ind1], days[ind2], days[ind3])
  schedule.index <- cbind(ind1, ind2, ind3)
  nmiss <- schedule.index[, 3] - schedule.index[, 2]
  maxmiss <- max(nmiss)
  if (maxmiss == 0) {
    pfind.si <- pk[1]
  } else if (maxmiss == 1) {
    pfind.si <- c(pk[1], (1 - pk[1]) * pk[2] * pk[1])
  } else {
    powk <- rep(pk[2], maxmiss + 1)
    powk[1] <- 1
    powk <- cumprod(powk)
    pfind.si <- pk[1] * powk * c(1, cumprod(1 - (pk[1] * powk[1:maxmiss])))
  }
  return(pfind.si)
}
