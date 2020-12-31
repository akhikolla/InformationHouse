#' @title Fit cp carcass persistence models
#' 
#' @description Carcass persistence is modeled as survival function where the 
#'   one or both parameter(s) can depend on any number of covariates. Format 
#'   and usage parallel that of common \code{R} functions such as \code{lm}, 
#'   \code{glm}, and \code{gam}. However, the input data (\code{data}) are 
#'   structured differently to accommodate the survival model approach (see 
#'   "Details"), and model formulas may be entered for both \code{l} 
#'   ("location") and \code{s} ("scale").
#'
#' @details The probability of a carcass persisting to a particular time is 
#'   dictated by the specific distribution chosen and its underlying location
#'   (l) and  scale (s) parameters (for all models except the exponential,  
#'   which only has a location parameter). Both \code{l} and \code{s} may 
#'   depend on covariates such as ground cover, season, species, etc., and a 
#'   separate model format (\code{formula_l} and \code{formula_s}) may be 
#'   entered for each. The models are entered as they would be in the familiar 
#'   \code{lm} or \code{glm} functions in R. For example, \code{l} might vary
#'   with \code{A}, \code{B}, and \code{C}, while \code{k} varies only with
#'  \code{A}. A user might then enter \code{p ~ A + B + C} for \code{formula_l}
#'   and \code{k ~ A} for \code{formula_s}. Other R conventions for defining
#'   formulas may also be used, with \code{A:B} for the interaction between
#'   covariates A and B and \code{A * B} as short-hand for \code{A + B + A:B}.
#'
#' Carcass persistence \code{data} must be entered in a data frame with data 
#'   in each row giving the fate of a single carcass in the trials. There
#'   must be a column for each of the last time the carcass was observed 
#'   present and the first time the carcass was observed absent (or NA if the
#'   carcass was always present). Additional columns with values for
#'   categorical covariates (e.g., visibility = E, M, or D) may also be 
#'   included.
#'
#' @param formula_l Formula for location; an object of class 
#'  "\code{\link{formula}}" (or one that can be coerced to that class):
#'  a symbolic description of the model to be fitted. Details of model 
#'  specification are given under "Details".
#'
#' @param formula_s Formula for scale; an object of class 
#'   "\code{\link{formula}}" (or one that can be coerced to that class):
#'   a symbolic description of the model to be fitted. Details of model 
#'   specification are given under "Details".
#'
#' @param data Data frame with results from carcass persistence trials and any
#'   covariates included in \code{formula_l} or {formula_s} (required).
#'
#' @param left Name of columns in \code{data} where the time of last present
#'   observation is stored.
#'
#' @param right Name of columns in \code{data} where the time of first absent
#'   observation is stored.
#'
#' @param dist Distribution name ("exponential", "weibull", "loglogistic", or 
#'   "lognormal")
#'
#' @param allCombos logical. If \code{allCombos = FALSE}, then the single model
#'  expressed by \code{formula_l} and \code{formula_s} is fit using a call to
#'  \code{cpm0}. If \code{allCombos = TRUE}, a full set of \code{\link{cpm}}
#'  submodels derived from combinations of the given covariates for \code{p}
#'  and \code{k} is fit. For example, submodels of \code{formula_l = p ~ A * B}
#'  would be \code{p ~ A * B}, \code{p ~ A + B}, \code{p ~ A}, \code{p ~ B},
#'  and \code{p ~ 1}. Models for each pairing of a \code{p} submodel with a
#' \code{k} submodel are fit via \code{cpmSet}, which fits each model
#'  combination using successive calls to \code{cpm0}, which fits a
#'  single model.
#'
#' @param sizeCol character string. The name of the column in \code{data} that
#'  gives the size class of the carcasses in the field trials. If
#'  \code{sizeCol = NULL}, then models are not segregated by size. If a
#'  \code{sizeCol} is provided, then separate models are fit for the \code{data}
#'  subsetted by \code{sizeCol}.
#'
#' @param CL confidence level
#'
#' @param quiet Logical indicator of whether or not to print messsages
#'
#' @return an object of an object of class \code{cpm}, \code{cpmSet},
#'  \code{cpmSize}, or \code{cpmSetSize}.
#' \describe{
#'  \item{\code{cpm0()}}{returns a \code{cpm} object, which is a description
#'    of a single, fitted pk model. Due to the large number and complexity of
#'    components of a\code{cpm} model, only a subset of them is printed
#'    automatically; the rest can be viewed/accessed via the \code{$} operator
#'    if desired. These are described in detail in the '\code{cpm} Components'
#'    section.}
#'  \item{\code{cpmSet()}}{returns a list of \code{cpm} objects, one for each
#'    of the submodels, as described with parameter \code{allCombos = TRUE}.}
#'  \item{\code{cpmSize()}}{returns a list of \code{cpmSet} objects (one for
#'    each 'size') if \code{allCombos = T}, or a list of \code{cpm} objects (one
#'    for each 'size') if \code{allCombos = T}}
#'  \item{\code{cpm}}{returns a \code{cpm}, \code{cpmSet}, \code{cpmSize}, or
#'    \code{cpmSetSize} object:
#'     \itemize{
#'        \item \code{cpm} object if \code{allCombos = FALSE, sizeCol = NULL}
#'        \item \code{cpmSet} object if \code{allCombos = TRUE, sizeCol = NULL}
#'        \item \code{cpmSize} object if \code{allCombos = FALSE, sizeCol != NULL}
#'        \item \code{cpmSetSize} object if \code{allCombos = TRUE, sizeCol != NULL}
#'     }
#'  }
#' }
#'
#' @section \code{cpm} Components:
#'
#' The following components of a \code{cpm} object are displayed automatically:
#'
#' \describe{
#'  \item{\code{call}}{the function call to fit the model}
#'  \item{\code{formula_l}}{the model formula for the \code{p} parameter}
#'  \item{\code{formula_s}}{the model formula for the \code{k} parameter}
#'  \item{\code{distribution}}{distribution used}
#'  \item{\code{predictors}}{list of covariates of \code{l} and/or \code{s}}
#'  \item{\code{AICc}}{the AIC value as corrected for small sample size}
#'  \item{\code{convergence}}{convergence status of the numerical optimization
#'    to find the maximum likelihood estimates of \code{p} and \code{k}. A
#'    value of \code{0} indicates that the model was fit successfully. For
#'    help in deciphering other values, see \code{\link{optim}}.}
#'  \item{\code{cell_ls}}{summary statistics for estimated cellwise
#'    \code{l} and \code{s}, including the medians and upper & lower bounds
#'    on CIs for each parameter, indexed by cell (or combination of
#'    covariate levels).}
#'  \item{\code{cell_ab}}{summary statistics for estimated cellwise
#'    \code{pda} and \code{pdb}, including the medians and upper & lower 
#'    bounds on CIs for each parameter, indexed by cell (or combination of
#'    covariate levels).}
#'  \item{\code{cell_desc}}{Descriptive statistics for estimated
#'    cellwise median persistence time and rI for search intervals of 1, 3, 7
#'    14, and 28 days, where rI is the probability of that carcass that arrives
#'    at a uniform random time in within a search interval of I days persists
#'    until the first search after arrival. }
#' }
#'
#' The following components are not printed automatically but can be accessed
#' via the \code{$} operator:
#' \describe{
#'  \item{\code{data}}{the data used to fit the model}
#'   \item{\code{betahat_l}}{parameter estimates for the terms in the 
#'     regression model for for \code{l}}
#'   \item{\code{betahat_s}}{parameter estimates for the terms in the 
#'     regression model for for \code{s}. If dist = "exponential", \code{s} 
#'     is set at 1 and not calculated.}
#'   \item{\code{varbeta}}{the variance-covariance matrix of the estimators
#'     for \code{c(betahat_l, betahat_s}.}
#'   \item{\code{cellMM_l}}{a cellwise model (design) matrix for covariate 
#'     structure of \code{l_formula}}
#'   \item{\code{cellMM_s}}{a cellwise model(design) matrix for covariate 
#'     structure of \code{s_formula}}
#'   \item{\code{levels_l}}{all levels of each covariate of \code{l}}
#'   \item{\code{levels_s}}{all levels of each covariate of \code{s}}
#'   \item{\code{nbeta_l}}{number of parameters fit for \code{l}}
#'   \item{\code{nbeta_s}}{number of parameters fit for \code{s}}
#'   \item{\code{cells}}{cell structure of the cp-model, i.e., combinations of
#'     all levels for each covariate of \code{p} and \code{k}. For example, if
#'     \code{covar1} has levels \code{"a"}, \code{"b"}, and \code{"c"}, and
#'     \code{covar2} has levels \code{"X"} and \code{"Y"}, then the cells 
#'     would consist of \code{a.X}, \code{a.Y}, \code{b.X}, \code{b.Y}, 
#'     \code{c.X}, and \code{c.Y}.}
#'  \item{\code{ncell}}{total number of cells}
#'  \item{\code{predictors_l}}{list of covariates of \code{l}}
#'  \item{\code{predictors_s}}{list of covariates of \code{s}}
#'  \item{\code{observations}}{observations used to fit the model}
#'  \item{\code{carcCells}}{the cell to which each carcass belongs}
#'  \item{\code{AIC}}{the 
#'    \href{https://en.wikipedia.org/wiki/Akaike_information_criterion}{AIC}
#'    value for the fitted model}
#'  \item{\code{CL}}{the input \code{CL}}
#'}
#' @section Advanced:
#'  \code{cpmSize} may also be used to fit a single model for each size class if
#'  \code{allCombos = FALSE}. To do so, \code{formula_l}, \code{formula_s}, and
#'  \code{dist} be named lists with names matching the sizes listed in
#'  \code{unique(data[, sizeCol])}. The return value is then a list of
#'  \code{cpm} objects, one for each size.
#'
#' @examples
#'  head(data(wind_RP))
#'  mod1 <- cpm(formula_l = l ~ Season, formula_s = s ~ 1, data = wind_RP$CP,
#'    left = "LastPresent", right = "FirstAbsent")
#'  class(mod1)
#'  mod2 <- cpm(formula_l = l ~ Season, formula_s = s ~ 1, data = wind_RP$CP,
#'    left = "LastPresent", right = "FirstAbsent", allCombos = TRUE)
#'  class(mod2)
#'  names(mod2)
#'  class(mod2[[1]])
#'  mod3 <- cpm(formula_l = l ~ Season, formula_s = s ~ 1, data = wind_RP$CP,
#'    left = "LastPresent", right = "FirstAbsent",
#'    allCombos = TRUE, sizeCol = "Size")
#'  class(mod3)
#'  names(mod3)
#'  class(mod3[[1]])
#'  class(mod3[[1]][[1]])
#'
#' @export
#'
cpm <- function(formula_l, formula_s = NULL, data, left, right,
    dist = "weibull", allCombos = FALSE, sizeCol = NULL,
    CL = 0.90, quiet = FALSE){

  if (is.null(allCombos) || is.na(allCombos) || !is.logical(allCombos)){
    stop("allCombos must be TRUE or FALSE")
  }
  if (is.null(sizeCol) || is.na(sizeCol) || sizeCol == F){
    if (!allCombos){ # single model
      out <- cpm0(formula_l = formula_l, formula_s = formula_s, data = data,
        left = left, right = right, dist = dist, CL = CL, quiet = quiet)
    } else { # allCombos of l and s subformulas
      out <- cpmSet(formula_l = formula_l, formula_s = formula_s, data = data,
        left = left, right = right, dist = dist, CL = CL, quiet = quiet)
    }
  } else { # specified formula for l and s, split by size class
    out <- cpmSize(formula_l = formula_l, formula_s = formula_s, data = data,
      left = left, right = right, dist = dist, sizeCol = sizeCol,
      allCombos = allCombos, CL = CL, quiet = quiet)
   }
  return(out)
}
#' @rdname cpm
#' @export
#'
cpm0 <- function(formula_l, formula_s = NULL, data = NULL, left = NULL,
    right = NULL, dist = "weibull", CL = 0.90, quiet = FALSE){
  dist <- tolower(dist)
  # initial error-checking
  ind <- sapply(data, is.factor)
  data[ind] <- lapply(data[ind], as.character)
  if (length(left) == 0){
    left <- "left"
    if (!"left" %in% colnames(data)){
      stop("No column name provided for first time observed ",
           "(left) and no column in data is named \"left\".")
    }
  } else if (length(left) > 1){
    stop("Input for first time observed column can only be length 0 or 1.")
  }
  if (!left %in% colnames(data)){
    stop("Column for last time observed (left) missing from data.")
  }
  if (length(right) == 0){
    right <- "right"
    if (!"right" %in% colnames(data)){
      stop("No column name provided for last time observed ",
           "(right) and no column in data is named \"right\".")
    }
  } else if (length(right) > 1){
    stop("Input for last time absent column can only be length 0 or 1.")
  }
  if (!right %in% colnames(data)){
    stop("Column name for last time absent (right) is not in the data.")
  }
  if (dist == "exponential"){
    if (!is.null(formula_s) && !quiet){
      message("Exponential distribution does not have a scale parameter. ",
              "formula_s ignored.")
    }
    formula_s <- formula(s ~ 1) # not used but s must have no predictors
  }
  formulaRHS_l <- formula(delete.response(terms(formula_l)))
  preds_l <- all.vars(formulaRHS_l)

  if (is.null(formula_s)) formula_s <- formula(s ~ 1)
  formulaRHS_s <- formula(delete.response(terms(formula_s)))
  preds_s <- all.vars(formulaRHS_s)

  if (!all(c(preds_l, preds_s) %in% colnames(data))){
    stop("Predictor in CP formula not found in CP data.")
  }
  preds <- unique(c(preds_l, preds_s))
  if (grepl("[-.]", paste0(preds, collapse = ''))){
    stop("Hyphen ( - ) and dot ( . ) not allowed in predictor names")
  }
  for (pri in preds){
    if (grepl("[-.]", paste0(data[, pri], collapse = '')))
      stop("Hyphen ( - ) and dot ( . ) not allowed in predictor levels")
  }
  if (anyNA(data[, left])){
    stop("NA not allowed for 'last time present' in CP data.")
  }
  if (any(data[ , left] > data[ , right], na.rm = TRUE)){
    stop("Some carcasses are observed after they have been removed?")
  }

  # in all cases, formula_l is used (with formula_s appended in some cases)
#  if (length(all.vars(formula_l)) == 1) mod_l <- l ~ 1
  if (length(all.vars(formula_l)) == 1) mod_l <- "1"
  mod_l <- paste(attr(terms(formula_l), "term.labels"), collapse = " + ")

  # parsing the formulas and covariates
  if (length(preds_l) > 0){
    for (predi in 1:length(preds_l)){
      data[ , preds_l[predi]] <- as.character(data[ , preds_l[predi]])
    }
  }
  levels_l <- .getXlevels(terms(formulaRHS_l), data)
  if (length(preds_s) > 0){
    for (predi in 1:length(preds_s)){
      data[ , preds_s[predi]] <- as.character(data[ , preds_s[predi]])
    }
  }
  levels_s <- .getXlevels(terms(formulaRHS_s), data)

  cells <- combinePreds(preds, data)
  ncell <- nrow(cells)
  cellNames <- cells$CellNames
  cellMM_l <- model.matrix(formulaRHS_l, cells)
  cellMM_s <- model.matrix(formulaRHS_s, cells)
  cellMM <- cbind(cellMM_l, cellMM_s)

  nbeta_l <- ncol(cellMM_l)
  nbeta_s <- ncol(cellMM_s)
  nbeta <- nbeta_l + nbeta_s

### 2n fix -->
  data0 <- data
  data0[is.na(data0[ , right]), right] <- Inf
  data00 <- data0
  fixBadCells <- NULL
  # if it is a full cell model (i.e., every combination of levels is [essentially]
  # fit indepedendently or model is 1, A, or A * B and not A + B), then apply the
  # 2n fix for any cell that has right == Inf.
  # NOTE: the "l" and the "s" are automatically crossed with each other, so we
  # just need to be sure that each is full cell within itself
  if ( # no predictors or all factor combinations are included in model
    (length(preds_l) == 0 || NCOL(cellMM_l) == prod(unlist(lapply(levels_l, length)))) &
    (length(preds_s) == 0 || NCOL(cellMM_s) == prod(unlist(lapply(levels_s, length))))){
    # then full cell model for both "l" and "s", so can check for right == Inf
    if (length(preds_l) == 0 & length(preds_s) == 0){# no predictors, no subsetting
      if (all(data[, right] == Inf)){ # 2n fix necessary
        if (dist != "exponential")
          stop("Must use exponential model when all data are right-censored")
        data0 <- rbind(data0, data0)
        dmed <- abs(data0[, left] - median(data0[, left]))
        i <- max(which(dmed == min(dmed)))
        data0[i, right] <- data0[i, left]
        fixBadCells <- "all"
      }
    } else { # have to check all cells
      for (ci in 1:nrow(cells)){
        ind <- which(matrixStats::colProds( # factor levels match cell
          t(data0[ , colnames(cells)[-ncol(cells)]]) ==
          as.character(cells[ci, -ncol(cells)])) == 1)
        if (all(data0[ind, right] == Inf)){ # then need to apply 2n fix for cell ci
          if (dist != "exponential")
            stop("Must use exponential model when ",
                 "all data in a cell are right-censored")
          tmp <- data0[ind, ]
          dmed <- abs(tmp[, left] - median(tmp[, left]))
          i <- max(which(dmed == min(dmed)))
          tmp[i, right] <- tmp[i, left]
          data0 <- rbind(data0, tmp)
          fixBadCells <- c(fixBadCells, cells$CellNames[ci])
        }
      }
    }
  } else {# additive model:
    # check factor levels and abort if right === Inf for any level
    for (pri in preds){
      for (li in unique(as.character(data[ ,pri]))){
        ind <- which(data0[, pri] == li)
        if (all(data0[ind, right] == Inf)){
          stop("Cannot fit additive model when ", right, " = Inf for all ",
               "carcasses in one level (", li, ") of a predictor (", pri, ")")
        }
      }
    }
  }
### <-- 2n fix
  # build the response variable (Surv object)
  t1 <- data0[ , left]
  t2 <- data0[ , right]
  event <- rep(3, length(t1)) # interval censored as default, but...
  event[is.na(t2) | is.infinite(t2)] <- 0 # study ends before removal
  event[round(t1, 3) == round(t2, 3)] <- 1 # carcass removal observed
  t1 <- pmax(t1, 0.0001)
  tevent <- survival::Surv(time = t1, time2 = t2, event = event, type = "interval")

  # traffic directing:
  #  use survreg if:
  #   1. dist == "exponential", or
  #   2. there is no "+" in formula_s and length(preds_s) <= 2
  #  otherwise, use custom MLE fitting function with optim
  if (dist == "exponential"){
    # use survreg with formula_l
    use_survreg <- T
    formula_cp <- reformulate(as.character(formulaRHS_l[-1]), response = "tevent")
  } else if ("+" %in% all.names(formula_s) || length(preds_s) > 2){
    # use custom fitting
    use_survreg <- F
  } else {
    # use survreg and convert formula_s to 'strata' format
    use_survreg <- T
    if (!is.language(formula_s) || length(all.vars(formula_s)) == 1){
      # scale formula not provided (assumed constant) or is set to constant
      formula_cp <- reformulate(as.character(formulaRHS_l[-1]),
        response = "tevent")
    } else {
      mod_s <- paste0(
        "strata(", paste(all.vars(formula_s)[-1], collapse = ", "), ")")
      formula_cp <- reformulate(paste(mod_l, mod_s, sep = " + "),
        response = "tevent")
    }
  }

  if (use_survreg){
    cpmod <- tryCatch(
      survival::survreg(formula = formula_cp, data = data0, dist = dist),
      error = function(x) NA, warning = function (x) NA
    )
    if (length(cpmod) == 1){
      stop("Failed CP model optimization.")
    }
    betahat_l <- cpmod$coefficients
    npreds_s <- length(all.vars(formula_s)) - 1
    if (npreds_s == 0){
      betahat_s <- log(cpmod$scale)
      varbeta <- cpmod$var
      carcCells <- rep('all', nrow(data0))
    } else if (npreds_s == 1){
      tmat <- diag(nbeta)
      tmat[-(1:(nbeta_l+1)), nbeta_l + 1] <- -1
      varbeta <- tmat %*% cpmod$var %*% t(tmat)
      betahat_s <- as.vector(tmat[-(1:nbeta_l), -(1:nbeta_l)] %*% log(cpmod$scale))
    } else if (npreds_s == 2) {
      nlev1 <- length(levels_s[[1]])
      nlev2 <- length(levels_s[[2]])
      nbeta_s <- nlev1 * nlev2
      tmat <- matrix(0, nrow = nbeta_s, ncol = nbeta_s)
      tmat[, 1] <- 1
      matInd1 <- cbind(rep(1:nlev1, each = nlev2), rep(1:nlev2, nlev1))
      interactInd <- expand.grid(1:(nlev1-1), 1:(nlev2-1)) + 1
      for (vvi in 1:nbeta_s){
        # additive terms
        i1 <- matInd1[vvi, 1]
        i2 <- matInd1[vvi, 2]
        tmat[vvi, i1] <- 1
        if (i2 > 1){
          tmat[vvi, nlev1 - 1 + i2] <- 1
          # interactions
          if (i1 > 1){
            tmat[vvi, nlev1 + nlev2 -1 + which(interactInd[,1] ==i1 & interactInd[,2] == i2)] <- 1
          }
        }
      }
      tmat <- solve(tmat)
      tM <- diag(nbeta)
      tM[-(1:nbeta_l), -(1:nbeta_l)] <- tmat
      varbeta <- tM %*% cpmod$var %*% t(tM)
      betahat_s <- as.vector(tmat %*% log(cpmod$scale))
    }
    varbeta_l <- varbeta[1:nbeta_l, 1:nbeta_l]
    cellMean_l <- cellMM_l %*% betahat_l
    cellVar_l <- cellMM_l %*% varbeta_l %*% t(cellMM_l)
    cellSD_l <- sqrt(diag(cellVar_l))

    if (dist == "exponential"){
      cellMean_s <- 0
      cellSD_s <- 0
    } else {
      which_s <- (nbeta_l + 1):nbeta
      varbeta_s <- varbeta[which_s, which_s]
      cellMean_s <- cellMM_s %*% betahat_s
      cellVar_s <- cellMM_s %*% varbeta_s %*% t(cellMM_s)
      cellSD_s <- sqrt(diag(cellVar_s))
    }

    probs <- data.frame(c(0.5, (1 - CL) / 2, 1 - (1 - CL) / 2))
    if (length(preds) > 0){
      carcCells <- data0[, preds[1]]
      if (length(preds) > 1){
        for (pri in 2:length(preds)){
          carcCells <- paste(carcCells, data0[, preds[pri]], sep = '.')
        }
      }
    }
    cell_n <- as.numeric(table(carcCells)[cellNames])
    cellTable_l <- apply(probs, 1, qnorm, mean = cellMean_l, sd = cellSD_l)
    cellTable_l <- round(matrix(cellTable_l, nrow = ncell, ncol = 3), 3)
    colnames(cellTable_l) <- c("l_median", "l_lwr", "l_upr")
    cellTable_s <- exp(apply(probs, 1, qnorm, mean = cellMean_s, sd = cellSD_s))
    cellTable_s <- round(matrix(cellTable_s, nrow = ncell, ncol = 3), 3)
    colnames(cellTable_s) <- c("s_median", "s_lwr", "s_upr")
    cellTable_ls <- data.frame(cell = cellNames, n = cell_n,
      cellTable_l, cellTable_s)
    if (dist == "exponential"){
      cellTable_a <- matrix("-", nrow = ncell, ncol = 3)
      colnames(cellTable_a) <- c("pda_median", "pda_lwr", "pda_upr")
      cellTable_b <- round(exp(cellTable_l), 3)
      colnames(cellTable_b) <- c("pdb_median", "pdb_lwr", "pdb_upr")
    }
    if (dist == "weibull"){
      cellTable_a <- round(1/cellTable_s, 3)[ , c(1, 3, 2)]
      dim(cellTable_a) <- dim(cellTable_s)
      #note: taking the reciprocal swaps the positions of the lwr and upr values
      colnames(cellTable_a) <- c("pda_median", "pda_lwr", "pda_upr")
      cellTable_b <- round(exp(cellTable_l), 3)
      dim(cellTable_b) <- dim(cellTable_l)
      colnames(cellTable_b) <- c("pdb_median", "pdb_lwr", "pdb_upr")
    }
    if (dist == "lognormal"){
      cellTable_a <- round(cellTable_s^2, 3)
      colnames(cellTable_a) <- c("pda_median", "pda_lwr", "pda_upr")
      cellTable_b <- round(cellTable_l, 3)
      colnames(cellTable_b) <- c("pdb_median", "pdb_lwr", "pdb_upr")
    }
    if (dist == "loglogistic"){
      cellTable_a <- round(1/cellTable_s, 3)[ , c(1, 3, 2)]
      dim(cellTable_a) <- dim(cellTable_s)
      #note: taking the reciprocal swaps the positions of the lwr and upr values
      colnames(cellTable_a) <- c("pda_median", "pda_lwr", "pda_upr")
      cellTable_b <- round(exp(cellTable_l), 3)
      dim(cellTable_b) <- dim(cellTable_l)
      colnames(cellTable_b) <- c("pdb_median", "pdb_lwr", "pdb_upr")
    }
    cellTable_ab <- data.frame(cell = cellNames, n = cell_n,
      cellTable_a, cellTable_b)
    if (dist == "exponential"){
      nbeta_s <- 0
    }
    output <- list()
    output$call <- match.call()
    output$data <- data
    output$formula_l <- formula_l
    if(dist != "exponential"){
      output$formula_s <- formula_s
    } else {
      output$formula_s <- NULL
    }
    output$distribution <- dist
    output$predictors <- preds
    output$predictors_l <- preds_l
    output$predictors_s <- preds_s
    output$AIC <- stats::AIC(cpmod)
    k <- ifelse(is.null(dim(cpmod$var)), 1, dim(cpmod$var)[1])
    n <- dim(cpmod$y)[1]
    output$AICc <- round(output$AIC + (2*k*(k + 1))/(n - k - 1), 2)
    output$convergence <- 0 # cpmod previously error-checked
    output$varbeta <- varbeta
    output$cellMM_l <- cellMM_l
    output$cellMM_s <- cellMM_s
    output$nbeta_l <- nbeta_l
    output$nbeta_s <- nbeta_s
    output$betahat_l <- betahat_l
    output$betahat_s <- betahat_s
    output$levels_l <- levels_l
    output$levels_s <- levels_s
    output$cells <- cells
    output$ncell <- ncell
    output$cell_ls <- cellTable_ls
    output$cell_ab <- cellTable_ab
    output$CL <- CL
    output$observations <- data[ , c(left, right)]
    output$carcCells <- carcCells
    output$loglik <- cpmod$loglik[2]
  } else if (!use_survreg){
    dataMM_l <- model.matrix(formulaRHS_l, data0)
    dataMM_s <- model.matrix(formulaRHS_s, data0)
    dataMM <- t(cbind(dataMM_l, dataMM_s))
    ncarc <- nrow(data)
    cellByCarc <- numeric(ncarc)
    for (celli in 1:ncell){
      groupPattern <- cellMM[celli, ]
      matchingMatrix <- dataMM == groupPattern
      matchingParts <- apply(matchingMatrix, 2, sum)
      matchingTotal <- matchingParts == ncol(cellMM)
      cellByCarc[matchingTotal] <- celli
    }
    carcCells <- cellNames[cellByCarc]
    cell_n <- as.numeric(table(carcCells)[cellNames])
    init_formRHS <- as.character(formulaRHS_l)[-1]
    init_form <- reformulate(init_formRHS, response = "tevent")
    init_mod <- survival::survreg(formula = init_form, data = data0, dist = dist)
    init_l <- init_mod$coef
    names(init_l) <- paste("l_", names(init_l), sep = "")
    init_s <- rep(init_mod$scale, nbeta_s)
    names(init_s) <- paste("s_", colnames(cellMM_s), sep = "")
    betaInit <- c(init_l, log(init_s))

    MLE <- tryCatch(
             optim(par = betaInit, fn = cpLogLik, method = "BFGS",
               t1 = t1, t2 = t2, cellMM = cellMM, dist = dist, hessian = TRUE,
               nbeta_l = nbeta_l, cellByCarc = cellByCarc, dataMM = dataMM,
               control = list(maxit = 1000)
             ), error = function(x) {NA}
           )

    if (length(MLE) == 1 && is.na(MLE)){
      stop("Failed optimization. Consider simplifying predictors.")
    }

    betahat <- MLE$par
    convergence <- MLE$convergence
    betaHessian <- MLE$hessian
    if (dist == "exponential"){
      which_s <- (nbeta_l + 1):nbeta
      betaHessian <- betaHessian[-which_s, -which_s]
    }
    llik <- -MLE$value

    nparam <- length(betahat)
    if (dist == "exponential"){
      nparam <- length(betahat) - 1
    }
    AIC <- 2 * nparam - 2 * llik
    AICc <- AIC + (2 * nparam * (nparam + 1)) / (ncarc - nparam - 1)
    AICc <- round(AICc, 2)
    AIC <- round(AIC, 2)

    betahat_l <- betahat[1:nbeta_l]
    names(betahat_l) <- colnames(dataMM_l)
    betahat_s <- betahat[(nbeta_l + 1):(nbeta)]
    names(betahat_s) <- colnames(dataMM_s)

    varbeta <- tryCatch(solve(betaHessian), error = function(x) {NA})
    if (is.na(varbeta)[1]){
      stop("Model generates unstable variance estimate.")
    }
    varbeta_l <- varbeta[1:nbeta_l, 1:nbeta_l]
    cellMean_l <- cellMM_l %*% betahat_l
    cellVar_l <- cellMM_l %*% varbeta_l %*% t(cellMM_l)
    cellSD_l <- sqrt(diag(cellVar_l))

    if (dist == "exponential"){
      cellMean_s <- 1
      cellSD_s <- 0
    } else {
      which_s <- (nbeta_l + 1):(nbeta)
      varbeta_s <- varbeta[which_s, which_s]
      cellMean_s <- cellMM_s %*% betahat_s
      cellVar_s <- cellMM_s %*% varbeta_s %*% t(cellMM_s)
      cellSD_s <- sqrt(diag(cellVar_s))
    }

    probs <- data.frame(c(0.5, (1 - CL)/2, 1 - (1 - CL)/2))
    cellTable_l <- apply(probs, 1, qnorm, mean = cellMean_l, sd = cellSD_l)
    cellTable_l <- round(matrix(cellTable_l, nrow = ncell, ncol = 3), 3)
    colnames(cellTable_l) <- c("l_median", "l_lwr", "l_upr")
    cellTable_s <- exp(apply(probs, 1, qnorm, mean = cellMean_s, sd = cellSD_s))
    cellTable_s <- round(matrix(cellTable_s, nrow = ncell, ncol = 3), 3)
    colnames(cellTable_s) <- c("s_median", "s_lwr", "s_upr")
    cellTable_ls <- data.frame(cell = cellNames, n = cell_n, cellTable_l, cellTable_s)

    if (dist == "exponential"){
      cellTable_a <- matrix("-", nrow = ncell, ncol = 3)
      colnames(cellTable_a) <- c("pda_median", "pda_lwr", "pda_upr")
      cellTable_b <- round(exp(cellTable_l), 3)
      colnames(cellTable_b) <- c("pdb_median", "pdb_lwr", "pdb_upr")
    }
    if (dist == "weibull"){
      cellTable_a <- round(1/cellTable_s, 3)
      colnames(cellTable_a) <- c("pda_median", "pda_lwr", "pda_upr")
      cellTable_b <- round(exp(cellTable_l), 3)
      colnames(cellTable_b) <- c("pdb_median", "pdb_lwr", "pdb_upr")
    }
    if (dist == "lognormal"){
      cellTable_a <- round(cellTable_s^2, 3)
      colnames(cellTable_a) <- c("pda_median", "pda_lwr", "pda_upr")
      cellTable_b <- round(cellTable_l, 3)
      colnames(cellTable_b) <- c("pdb_median", "pdb_lwr", "pdb_upr")
    }
    if (dist == "loglogistic"){
      cellTable_a <- round(1/cellTable_s, 3)
      colnames(cellTable_a) <- c("pda_median", "pda_lwr", "pda_upr")
    cellTable_b <- round(exp(cellTable_l), 3)
      colnames(cellTable_b) <- c("pdb_median", "pdb_lwr", "pdb_upr")
    }
    cellTable_ab <- data.frame(cell = cellNames, n = cell_n, cellTable_a, cellTable_b)


    if (dist == "exponential"){
      nbeta_s <- 0
    }
    output <- list()
    output$call <- match.call()
    output$data <- data
    if (!is.null(fixBadCells)){
      output$data_adj <- data0
      output$fixBadCells <- fixBadCells
    }
    output$formula_l <- formula_l
    if(dist != "exponential"){
      output$formula_s <- formula_s
    } else {
      output$formula_s <- NULL
    }
    output$distribution <- dist
    output$predictors <- preds
    output$predictors_l <- preds_l
    output$predictors_s <- preds_s
    output$AIC <- AIC
    output$AICc <- AICc
    output$convergence <- 0
    output$varbeta <- varbeta
    output$cellMM_l <- cellMM_l
    output$cellMM_s <- cellMM_s
    output$nbeta_l <- nbeta_l
    output$nbeta_s <- nbeta_s
    output$betahat_l <- betahat_l
    output$betahat_s <- betahat_s
    output$levels_l <- levels_l
    output$levels_s <- levels_s
    output$cells <- cells
    output$ncell <- ncell
    output$cell_ls <- cellTable_ls
    output$cell_ab <- cellTable_ab
    output$CL <- CL
    output$observations <- data0[ , c(left, right)]
    output$carcCells <- carcCells
    output$loglik <- llik

  }

  Ir <- c(1, 3, 7, 14, 28) # search intervals for calculating r
  t0 <- numeric(length(Ir))
  t1 <- c(1, 3, 7, 14, 28)
  tf <- t1
  cell_desc <- matrix(nrow = ncell, ncol = 2 + length(Ir))
  colnames(cell_desc) <- c("cell", "medianCP", paste0("r", Ir))
  pda <- output$cell_ab[ , "pda_median"]
  pdb <- output$cell_ab[ , "pdb_median"]
  cell_desc[ , 2 + 1:length(Ir)] <- t(ppersist(pda = pda, pdb = pdb,
    dist = dist, t_arrive0 = t0, t_arrive1 = t1, t_search = tf))
  if (dist == "weibull"){
    cell_desc[ , "medianCP"] <- pdb * log(2)^(1/pda)
  } else if (dist == "lognormal"){
    cell_desc[ , "medianCP"] <- exp(pdb)
  } else if (dist == "loglogistic"){
    cell_desc[ , "medianCP"] <- pdb
  } else if (dist == "exponential"){
    cell_desc[ , "medianCP"] <- log(2) * pdb
  }
  output$cell_desc <- data.frame(cell_desc)
  output$cell_desc[ , "cell"] <- cellNames
   class(output) <- c("cpm", "list")
  attr(output, "hidden") <- c("data", "predictors_l", "predictors_s",
    "betahat_l", "betahat_s", "cellMM_l", "cellMM_s", "nbeta_l", "nbeta_s",
    "varbeta", "levels_l", "levels_s", "carcCells", "AIC", "cells", "ncell",
    "observations", "loglik")
  return(output)
}

#' @rdname cpm
#' @export
cpmSet <- function(formula_l, formula_s = NULL, data, left, right,
    dist = c("exponential", "weibull", "lognormal", "loglogistic"),
    CL = 0.90, quiet = FALSE){

  if (length(formula_s) == 0){
    formula_s <- formula(s ~ 1)
  }

  terms_l <- attr(terms(formula_l), "term.labels")
  terms_s <- attr(terms(formula_s), "term.labels")
  nterms_l <- length(terms_l)
  nterms_s <- length(terms_s)
  nformula_l <- 2^(nterms_l)
  nformula_s <- 2^(nterms_s)

  dropComplex_l <- rep(1:nterms_l, choose(nterms_l, 1:nterms_l))
  dropWhich_l <- numeric(0)
  if (nterms_l > 0){
    for (termi in 1:nterms_l){
      specificDrop <- seq(1, choose(nterms_l, (1:nterms_l)[termi]))
      dropWhich_l <- c(dropWhich_l, specificDrop)
    }
  }
  optionFormula_l <- vector("list", nformula_l)
  optionFormula_l[[1]] <- formula_l
  keepFormula_l <- rep(TRUE, nformula_l)
  if (nformula_l > 1){
    for (formi in 2:nformula_l){
      termDropComplex <- combn(terms_l, dropComplex_l[formi - 1])
      termDropSpec <- termDropComplex[ , dropWhich_l[formi - 1]]
      termDrop <- paste(termDropSpec, collapse = " - ")
      formulaUpdate <- paste(format(~.), "-", termDrop)
      updatedFormula <- update.formula(formula_l, formulaUpdate)
      optionFormula_l[[formi]] <- updatedFormula
      keepFormula_l[formi] <- checkComponents(updatedFormula)
    }
    nkeepFormula_l <- sum(keepFormula_l)
    whichKeepFormula_l <- which(keepFormula_l == TRUE)
    keptFormula_l <- vector("list", nkeepFormula_l)
    for (kepti in 1:nkeepFormula_l){
      keptFormula_l[[kepti]] <- optionFormula_l[[whichKeepFormula_l[kepti]]]
    }
  }else{
    keptFormula_l <- optionFormula_l
  }

  dropComplex_s <- rep(1:nterms_s, choose(nterms_s, 1:nterms_s))
  dropWhich_s <- numeric(0)
  if (nterms_s > 0){
    for (termi in 1:nterms_s){
      specificDrop <- seq(1, choose(nterms_s, (1:nterms_s)[termi]))
      dropWhich_s <- c(dropWhich_s, specificDrop)
    }
  }
  optionFormula_s <- vector("list", nformula_s)
  optionFormula_s[[1]] <- formula_s
  keepFormula_s <- rep(TRUE, nformula_s)
  if (nformula_s > 1){
    for (formi in 2:nformula_s){
      termDropComplex <- combn(terms_s, dropComplex_s[formi - 1])
      termDropSpec <- termDropComplex[ , dropWhich_s[formi - 1]]
      termDrop <- paste(termDropSpec, collapse = " - ")
      formulaUpdate <- paste(format(~.), "-", termDrop)
      updatedFormula <- update.formula(formula_s, formulaUpdate)
      optionFormula_s[[formi]] <- updatedFormula
      keepFormula_s[formi] <- checkComponents(updatedFormula)
    }
    nkeepFormula_s <- sum(keepFormula_s)
    whichKeepFormula_s <- which(keepFormula_s == TRUE)
    keptFormula_s <- vector("list", nkeepFormula_s)
    for (kepti in 1:nkeepFormula_s){
      keptFormula_s[[kepti]] <- optionFormula_s[[whichKeepFormula_s[kepti]]]
    }
  }else{
    keptFormula_s <- optionFormula_s
  }

  expandedKeptFormulae <- expand.grid(keptFormula_l, keptFormula_s, dist)
  keptFormula_l <- expandedKeptFormulae[ , 1]
  keptFormula_s <- expandedKeptFormulae[ , 2]
  dist <- as.character(expandedKeptFormulae[ , 3])
  nmod <- nrow(expandedKeptFormulae)
  preoutput <- vector("list", nmod)
  for (modi in 1:nmod){
    formi_l <- keptFormula_l[modi][[1]]
    formi_s <- keptFormula_s[modi][[1]]
    disti <- dist[modi]
    if (disti == "exponential"){
      formi_s <- NULL
    }
    cpm_i <- tryCatch(
      cpm0(formi_l, formi_s, data = data, left = left, right = right,
        dist = disti, CL = CL, quiet = quiet),
      error = function(x) {
        paste("Failed model fit: ", geterrmessage(), sep = "")
      }
    )
    name_d <- disti
    name_l <- paste(format(formi_l), collapse = "")
    name_l <- gsub("    ", "", name_l)
    name_s <- paste(format(formi_s), collapse = "")
    name_s <- gsub("    ", "", name_s)
    modName <- paste("dist: ", name_d, "; ", name_l, "; ", name_s, sep = "")

    preoutput[[modi]] <- cpm_i
    names(preoutput)[modi] <- modName
  }
  uniqueMods <- unique(names(preoutput))
  nuniqueMods <- length(uniqueMods)
  output <- vector("list", nuniqueMods)
  names(output) <- uniqueMods
  for (modi in 1:nuniqueMods){
    output[[modi]] <- preoutput[[uniqueMods[modi][1]]]
  }
  class(output) <- c("cpmSet", "list")
  return(output)
}
#' @rdname cpm
#' @export
cpmSize <- function(formula_l, formula_s = NULL, data, left, right,
  dist = c("exponential", "weibull", "lognormal", "loglogistic"),
  sizeCol = NULL, allCombos = FALSE, CL = 0.90, quiet = FALSE){

  if (length(sizeCol) == 0 || anyNA(sizeCol)){
    out <- list()
    if (!allCombos){
      if (length(dist) > 1){
        stop("Ambiguous cpm call: length(dist) > 1)")
      }
      out[["all"]] <- cpm0(formula_l = formula_l, formula_s = formula_s,
        data = data, left = left, right = right, dist = dist,
        CL = CL, quiet = quiet)
    } else {
      out[["all"]] <- cpmSet(formula_l = formula_l, formula_s = formula_s,
        data = data, left = left, right = right, dist = dist,
        CL = CL, quiet = quiet)
    }
    class(out) <- c(ifelse(allCombos, "cpmSetSize", "cpmSize"), "list")
    return(out)
  }

  if (!(sizeCol %in% colnames(data))){
    stop("sizeCol not in data set.")
  }


  sizeclasses <- sort(unique(as.character(data[ , sizeCol])))
  if (!allCombos){
    if ("list" %in% c(class(formula_l), class(formula_s), class(dist))){
      if (!("list" %in% intersect(intersect(class(formula_l), class(formula_s)), class(dist)))){
        stop("formula_l, formula_s, and dist must be parallel lists, or ",
             "formula_l and formula_s must be 'scalars'")
      } else {
      # parallel lists of formulas and distributions
      # fit one model for each carcass class
      # then fit the specific models for each formula and corresponding size
        if (!setequal(names(formula_l), names(formula_s)) ||
            !setequal(names(formula_l), sizeclasses) ||
            !setequal(names(formula_l), names(dist))){
          stop("l and s formula names and list of distributions names ",
              "must match carcass classes")
        }
        out <- list()
        for (szi in sizeclasses){
          out[[szi]] <- cpm0(
            formula_l = formula_l[[szi]], formula_s = formula_s[[szi]],
            dist = dist[[szi]],
            data = data[data[ , sizeCol] == szi, ], left = left, right = right,
            CL = CL, quiet = quiet
          )
        }
        class(out) <- c("cpmSize", "list")
        return(out)
      }
    } else { # just one formula and one distribution shared by all sizes
      if (is.null(dist)){
        dist <- "weibull"
      } else if (length(dist) > 1){
        stop("If allCombos == FALSE, then dist must be NULL or a single value ",
             "if there is only one formula_l, or dist must be a list if ",
             "formula_l is."
        )
      } else {
        out <- list()
        for (szi in sizeclasses){
          out[[szi]] <- cpm0(formula_l = formula_l, formula_s = formula_s,
            dist = dist, data = data[data[ , sizeCol] == szi, ],
            left = left, right = right, CL = CL, quiet = quiet)
        }
        class(out) <- c("cpmSize", "list")
        return(out)
      }
    }
  } else {
    if ("list" %in% c(class(formula_l), class(formula_s), class(dist))){
      stop(
        "If allCombos == TRUE, then formula_l and formula_s must be 'scalars'"
      )
    }
    out <- list()
    for (sci in sizeclasses){
      out[[sci]] <- cpmSet(formula_l = formula_l, formula_s = formula_s,
        data = data[data[, sizeCol] == sci, ], left = left, right = right,
        dist = dist, CL = CL, quiet = quiet
      )
    }
    class(out) <- c("cpmSetSize", "list")
    return(out)
  }
}

#' @title Print a \code{\link{cpm}} model object
#'
#' @description Print a \code{\link{cpm}} model object
#'
#' @param x a \code{\link{cpm}} model object
#'
#' @param ... to be passed down
#'
#' @export
#'
print.cpm <- function(x, ...){
  hid <- attr(x, "hidden")
  notHid <- !names(x) %in% hid
  print(x[notHid])
}
 
#' @title Calculate the negative log-likelihood of a carcass persistence model
#' 
#' @description The function used to calculate the negative-loglikelihood of
#'   a given carcass persistence model (\code{\link{cpm}}) with a given data
#'   set
#'
#' @param t1 last times observed present
#'
#' @param t2 first times observed absent
#'
#' @param beta Parameters to be optimized.
#'
#' @param nbeta_l Number of parameters associated with l.
#'
#' @param cellByCarc Which cell each observation belongs to.
#'
#' @param cellMM Combined model matrix.
#'
#' @param dataMM Combined model matrix expanded to the data.
#'
#' @param dist Name of distribution. 
#'
#' @return Negative log likelihood of the observations, given the parameters.
#'
#' @export 
#'
cpLogLik <- function(t1, t2, beta, nbeta_l, cellByCarc, cellMM, dataMM, dist){
  # scaling for beta 
  t2[which(is.na(t2))] <- Inf
  ncell <- nrow(cellMM)
  nbeta <- length(beta)
  nbeta_s <- nbeta - nbeta_l
  which_l <- 1:nbeta_l
  which_s <- (nbeta_l + 1):nbeta
  beta_l <- beta[which_l]
  beta_s <- beta[which_s]
  if (dist == "exponential"){
    beta_s <- c(log(1), rep(0, nbeta_s - 1))
  } 
  dataMM_l <- matrix(dataMM[which_l, ], ncol = nbeta_l, byrow = TRUE)
  dataMM_s <- matrix(dataMM[which_s, ], ncol = nbeta_s, byrow = TRUE)
  Beta_l <- dataMM_l %*% beta_l 
  Beta_s <- dataMM_s %*% beta_s  
  psurv_t1 <- survival::psurvreg(t1, Beta_l, exp(Beta_s), dist)
  psurv_t2 <- survival::psurvreg(t2, Beta_l, exp(Beta_s), dist)
  psurv_t2[which(is.na(psurv_t2))] <- 1
  lik <- psurv_t2 - psurv_t1
  too_small <- (t1 + 0.0001) >= t2
  if (any(too_small)){
    lik[too_small] <- survival::dsurvreg(t2[too_small], Beta_l[too_small],
      exp(Beta_s)[too_small], dist)
  }
  lik <- pmax(lik, .Machine$double.eps) 
  nll <- -sum(log(lik))
  return(nll)
}


#' @title Simulate parameters from a fitted cp model
#'
#' @description Simulate parameters from a \code{\link{cpm}} model object, and 
#'   format them as either type \code{"survreg"} or \code{"ppersist"}
#'
#' @param n the number of simulation draws
#'
#' @param model A \code{cpm} object (which is returned from 
#'   \code{\link{cpm}})
#'
#' @param type The type of parameters requested. \code{"survreg"} or 
#'   \code{"ppersist"}
#'
#' @return list of two matrices of \code{n} simulated \code{l} and \code{s}
#'   (if \code{type = "survreg"}) or \code{a} and \code{b} (if \code{type = 
#'   "ppersist"})for cells defined by the \code{model} object. 
#'
#' @examples
#'   data(wind_RP)
#'   mod <- cpm(formula_l = l ~ 1, data = wind_RP$CP, left = "LastPresent",
#'            right = "FirstAbsent"
#'          )
#'   rcp(n = 10, model = mod, type = "survreg")
#'   rcp(n = 10, model = mod, type = "ppersist")
#'
#' @export
#'
rcp <- function(n, model, type = "survreg"){

  if (!"cpm" %in% class(model)){
    stop("model not of class cpm.")
  }
  if (!type %in% c("survreg", "ppersist")){
    stop(paste("type ", type, " is not supported.", sep = ""))
  }
  dist <- model$dist
  nbeta_l <- model$nbeta_l 
  nbeta_s <- model$nbeta_s
  if (dist == "exponential"){
    nbeta_s <- 1
  } 
  which_beta_s <- (nbeta_l + 1):(nbeta_l + nbeta_s)
  cellMM_l <- model$cellMM_l
  cellMM_s <- model$cellMM_s
  ncell <- model$ncell
  cellNames <- model$cells[ , "CellNames"]
  meanbeta <- c(model$betahat_l, model$betahat_s)
  if (dist == "exponential"){
    varbeta <- rbind(cbind(model$varbeta, 0), 0)
  } else{
    varbeta <- model$varbeta
  }
  method <-  "svd"

  sim_beta <- mvtnorm::rmvnorm(n, mean = meanbeta, sigma = varbeta, method)

  sim_l <- as.matrix(sim_beta[ , 1:nbeta_l] %*% t(cellMM_l))  # coef
  sim_s <- exp(as.matrix(sim_beta[ , which_beta_s] %*% t(cellMM_s))) # scale

  if (type == "ppersist"){
    if (dist == "exponential"){
      sim_a <- matrix(NA, nrow = n, ncol = ncell)
      sim_b <- exp(sim_l)
    }
    if (dist == "weibull"){
      sim_a <- 1 / sim_s
      sim_b <- exp(sim_l)
    }
    if (dist == "lognormal"){
      sim_a <- sim_s^2
      sim_b <- sim_l
    }
    if (dist == "loglogistic"){
      sim_a <- 1 / sim_s
      sim_b <- exp(sim_l)
    } 
    sim_p1 <- sim_a
    sim_p2 <- sim_b
  } else{
    sim_p1 <- sim_l
    sim_p2 <- sim_s    
  }

  colnames(sim_p1) <- cellNames
  colnames(sim_p2) <- cellNames


  paramNames <- switch(type, 
                  "survreg" = c("l", "s"), "ppersist" = c("pda", "pdb"))

  output <- vector("list", ncell)
  names(output) <- cellNames
  for (celli in 1:ncell){
    cellp12 <- cbind(sim_p1[ , celli], sim_p2[ , celli])
    colnames(cellp12) <- paramNames
    output[[celli]] <-  cellp12
  }

  return(output)
}

#' @title Create the AICc tables for a set of carcass persistence models
#' 
#' @description S3 function to generate model comparison tables based on AICc
#'  values for a set of CP models generated by \code{\link{cpmSet}}
#'
#' @param x Set of carcass persistence models fit to the same
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
#'   mod <- cpmSet(formula_l = l ~ Season * Visibility, formula_s = s ~ Season,
#'            data = wind_RP$CP, left = "LastPresent", right = "FirstAbsent")
#'  aicc(mod)
#'
#' @export 
#'
aicc.cpmSet <- function(x, ... , quiet = FALSE, app = FALSE){
  cpmset <- x
  nmod <- length(cpmset)
  formulas <- names(cpmset)
  dist <- rep(NA, nmod)
  formulas_l <- rep(NA, nmod)
  formulas_s <- rep(NA, nmod)
  AICc <- rep(NA, nmod)
  deltaAICc <- rep(NA, nmod)

  if (nmod == 1){
    splitFormulas <- strsplit(formulas, "; ")[[1]]
    dist <- strsplit(splitFormulas[1], "dist: ")[[1]][2]
    formulas_l <- splitFormulas[2] 
    formulas_s <- splitFormulas[3]
    AICc <- tryCatch(cpmset[[1]]$AICc, error = function(x) {1e7})
    deltaAICc <- 0    
    AICcOrder <- 1
  } else {
    for (modi in 1:nmod){
      splitFormulas_i <- strsplit(formulas[modi], "; ")[[1]]
      dist[modi] <- strsplit(splitFormulas_i, "dist: ")[[1]][2]
      formulas_l[modi] <- splitFormulas_i[2] 
      formulas_s[modi] <- splitFormulas_i[3]
      AICc[modi] <- tryCatch(cpmset[[modi]]$AICc, error = function(x) {1e7})
    }
    AICc <- round(AICc, 2)
    AICcOrder <- order(AICc)
    deltaAICc <- round(AICc - min(AICc), 2)
    which_fails <- which(AICc == 1e7)
    AICc[which_fails] <- NA
    deltaAICc[which_fails] <- NA
  }

  if (app){
    formulas_l <- gsub("~ 1", "~ constant", formulas_l)
    formulas_s <- gsub("~ 1", "~ constant", formulas_s)
  }
  output <- data.frame(dist, formulas_l, formulas_s, AICc, deltaAICc)
  output <- output[AICcOrder, ]
  colnames(output) <- c("Distribution", "Location Formula", "Scale Formula", 
                        "AICc", "\u0394AICc"
                       )
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
  return(output)
}

#' @title Create the AICc tables for a list of sets of searcher efficiency models
#'
#' @description S3 function to generate model comparison tables for lists of
#'  of sets of CP models of class \code{\link[=cpm]{cpmSetSize}}
#'
#' @param x List of sets of CP models fit to the same observations
#'
#' @param ... further arguments passed to or from other methods
#'
#' @return AICc table
#'
#' @examples
#'  cpmods <- cpm(formula_l = l ~ Visibility, data = wind_RP$CP,
#'    left = "LastPresent", right = "FirstAbsent", sizeCol = "Size",
#'    allCombos = TRUE)
#'  aicc(cpmods)
#'
#' @export
#'
aicc.cpmSetSize <- function(x, ... ){
  return(lapply(x, FUN = function(y){
    class(y) <- c("cpmSet", "list")
    aicc(y)
  }))
}

#' @title Extract AIC and AICc for a carcass persistence model
#'
#' @description S3 function for generating AIC for \code{\link{cpm}} objects
#'
#' @param x Carcass persistence model (\code{cpm} objects)
#'
#' @param ... further arguments passed to or from other methods
#'
#' @return AIC, AICc vector
#'
#' @examples
#'   data(wind_RP)
#'   mod <- cpm(formula_l = l ~ Season, formula_s = s ~ Season,
#'            data = wind_RP$CP, left = "LastPresent", right = "FirstAbsent")
#'  aicc(mod)
#'
#' @export
#'
aicc.cpm <- function(x,...){
  if (x$distribution == "exponential"){
    formula_s = NULL
  } else {
    formula_s = deparse(x$formula_s)
  }
  return(
    data.frame(cbind(
      "formula_l" = deparse(x$formula_l),
      "formula_s" = formula_s,
      "dist" = x$distribution,
      "AICc" = x$AICc
    ))
  )
}

#' @title Calculate the probability of persistence to detection
#' 
#' @description Given a set of CP parameters (of \code{"ppersist"} type), 
#'   calculate the probability of persistence to detection for a carcass.
#' 
#' @param pda parameter a.
#' 
#' @param pdb parameter b.
#' 
#' @param dist Distribution used.
#' 
#' @param t_arrive0 Beginning of arrival window.
#' 
#' @param t_arrive1 End of arrival window.
#' 
#' @param t_search Search time.
#' 
#' @return Probability of persistence of detection to at t_search, 
#'          given arrival between t_arrive0 and t_arrive1
#' 
#' @export 
#'
ppersist <- function(pda, pdb, dist, t_arrive0, t_arrive1, t_search){
  dist <- tolower(dist)
  if (dist == "weibull"){
    pda[log(pda) < -5] <- exp(-5) # adjustment to avoid overflow errors
    pda[log(pda) > 5] <- exp(5)   # adjustment to avoid overflow errors
    sa0 <- pgamma(outer(1/pdb, t_search - t_arrive0)^pda, 1/pda, log.p = TRUE)
    sa1 <- pgamma(outer(1/pdb, t_search - t_arrive1)^pda, 1/pda, log.p = TRUE)
    a1a0 <- outer(pdb, 1/(t_arrive1 - t_arrive0))
    probs <- (exp(sa0) - exp(sa1)) * gamma(1 + 1/pda) * a1a0
    probs <- t(probs)
  } else if (dist == "exponential"){
    a1a0 <- outer(t_arrive1 - t_arrive0, 1/pdb)
    a0s <- outer(t_arrive0 - t_search, 1/pdb)
    a1s <- outer(t_arrive1 - t_search, 1/pdb)
    probs <- (exp(a1s) - exp(a0s))/(a1a0)
  } else if (dist == "lognormal"){
    root_pda <- sqrt(pda)
    exp_value <- exp((pda / 2) + pdb)
    tt <- t_search - t_arrive0
    p1 <- exp(pnorm(outer(pdb, -log(tt), "+") / root_pda, log.p = TRUE))
    p2 <- exp(pnorm(outer(-pdb, log(tt), "+") / root_pda - root_pda, log.p = TRUE)) * exp_value
    part0 <- t(p1) * tt + t(p2)
    tt <- t_search - t_arrive1
    p1 <- exp(pnorm(outer(pdb, -log(tt), "+") / root_pda, log.p = TRUE))
    p2 <- exp(pnorm(outer(-pdb, log(tt), "+") / root_pda - root_pda, log.p = TRUE)) * exp_value
    part1 <- t(p1) * tt + t(p2)
    probs <- -(part1 - part0) / (t_arrive1 - t_arrive0)
    probs[is.na(probs) & pdb > 0] <- 1
    probs[is.na(probs) & pdb <= 0] <- 0
  } else if (dist == "loglogistic" | dist == "log-logistic"){
    yox <- function(x, y) y/x
    t1 <- t_search-t_arrive1
    t0 <- t_search-t_arrive0
    tob <- outer(pdb, t1, "yox")
    part1 <- t1/t(1 + tob^pda) * 
      t(gsl::hyperg_2F1(1, 1, 1 + 1/pda, 1/(1 + tob^(-pda))))
    tob <- outer(pdb, t0, "yox")
    part0 <- t0 / t(1 + tob^pda) *
      t(gsl::hyperg_2F1(1, 1, 1 + 1/pda, 1/(1 + tob^(-pda))))
    probs <- (part0 - part1)/(t_arrive1 - t_arrive0)
    # correction for overflow errors
    probs[pllogis(q = rep(t1, length(pda)),
      pda = rep(pda, each = length(t1)),
      pdb = rep(pdb, each = length(t1))) > 1 - 1e-7] <- 0
    probs[is.na(probs)] <- 0
  }
  return(probs)
}

#' @title Check if a CP model is well-fit
#'
#' @description Run a check the arg is a well-fit cpm object
#'
#' @param cpmod A \code{\link{cpm}} object to test
#'
#' @return logical value indicating a failed fit (TRUE) or successful (FALSE)
#'
#' @export
#'
cpmFail <- function(cpmod){
!("cpm" %in% class(cpmod)) || anyNA(cpmod) || sum(diag(cpmod$varbeta) < 0) > 0
}

#' @title Check if cpm models fail
#' 
#' @description Run a check on each model within a \code{\link{cpmSet}} object
#'   to determine if it failed or not
#'
#' @param cpmSetToCheck A \code{\link{cpmSet}} object to test
#'
#' @return A vector of logical values indicating if each of the models failed
#'
#' @export
#'
cpmSetFail <- function(cpmSetToCheck){
  unlist(lapply(cpmSetToCheck, cpmFail)) 
}

#' @title Check if all of the cpm models fail
#'
#' @description Run a check on each model within a \code{\link[=cpm]{cpmSetSize}}
#'   object to determine if they all failed or not
#'
#' @param cpmSetSizeToCheck A \code{cpmSetSize} object to test
#'
#' @return A list of vectors of logical values indicating if each of the 
#'   models failed
#'
#' @export
#'
cpmSetSizeFail <- function(cpmSetSizeToCheck){
  lapply(cpmSetSizeToCheck, cpmSetFail)
}

#' @title Remove failed cpm models from a \code{\link{cpmSet}} object
#'
#' @description Remove all failed models within a \code{\link{cpmSet}} object
#'
#' @param cpmSetToTidy A \code{\link{cpmSet}} object to tidy
#'
#' @return A \code{\link{cpmSet}} object with failed models removed
#'
#' @export
#'
cpmSetFailRemove <- function(cpmSetToTidy){
  out <- cpmSetToTidy[!cpmSetFail(cpmSetToTidy)]
  class(out) <- c("cpmSet", "list")
  return(out)
}

#' @title Remove failed cpm models from a \code{cpmSetSize} object
#'
#' @description Remove failed models from a \code{\link[=cpm]{cpmSetSize}} object
#'
#' @param cpmSetSizeToTidy A list of \code{cpmSetSize} objects to tidy
#'
#' @return A list of \code{\link{cpmSet}} objects with failed models removed
#'
#' @export
#'
cpmSetSizeFailRemove <- function(cpmSetSizeToTidy){
  out <- list()
  for (sci in names(cpmSetSizeToTidy)){
    out[[sci]] <- cpmSetFailRemove(cpmSetSizeToTidy[[sci]])
  }
  class(out) <- c("cpmSetSize", "list")
  return(out)
}

#' @export
aicc.cpmSize <- function(x, ... ){
  return(lapply(x, FUN = function(y){
    class(y) <- c("cpm", "list")
    aicc(y)
  }))
}

#' @title Descriptive statistics for a fitted CP model
#'
#' @description Given a \code{cpm} object, calculate convenient descriptive statistics,
#'  including the median CP, specified \code{rI} statistics, and \code{pda}
#'  and \code{pdb} statistics for the fitted model (EoA parameterization), and
#'  location and scale parameters for the fitted model (\code{survival} package
#'  parameterization) along with estimated CIs.
#'
#' @details The CIs for the r statistics (and the medianCP for the Weibull) ara
#'  based on simulation of the \code{pda} and \code{pdb} parameters, calculation
#'  of the statistics, and taking the empirical distribution of the simulated
#'  values. Other CIs are based on the assumed bivariate normal distributions of
#'  the appropriately transformed \code{l} and \code{s} parameters in the fitted
#'  model using \code{beta_hat} and \code{varbeta}.
#'
#'  NOTE: \code{rI} is the probability that a carcass that arrives at a uniform random
#'  time in an interval of \code{I} days will persist until the first search after
#'  arrival.
#'
#' @param model_CP A fitted CP model (\code{cpm} object)
#'
#' @param Ir The intervals for which to calculate the r statistics
#'
#' @param CL The confidence level for the CIs.
#'
#' @param nsim Number of simulation draws for estimating CIs
#'
#' @return Matrix of point and interval estimates for the median CP and the r
#'  statistics for the specified intervals. The matrix is assigned to class
#'  \code{descCP} that is simply a matrix with dimensions
#'  \code{ncell x (1 + 3*(5 + length(Ir)))}, column names that give the number of
#'  observations in each cell, statistic name and upper and lower bounds
#' (in triplets), and row names giving the names of the cells. \code{CL}, \code{nsim},
#'  and the name of the fitted model (\code{model_CP}) are included as object
#'  attributes.
#'
#' @seealso \code{\link{cpm}}, \code{\link{rcp}}, \code{\link{ppersist}}
#'
#' @export
#'
desc <- function(model_CP, Ir = c(1, 3, 7, 14, 28), CL = 0.9, nsim = 10000){
  if (!"cpm" %in% class(model_CP)) stop("model_CP must be a cpm object.")
  # set up summary table
  Ir <- sort(Ir)
  t0 <- numeric(length(Ir))
  t1 <- Ir
  tf <- t1
  Irv <- paste0("r", Ir)
  cell_desc <- matrix(nrow = model_CP$ncell, ncol = 1 + 3 * (5 + length(Ir)))
  Irvec <- c(Irv, paste0(Irv, "_lwr"), paste0(Irv, "_upr"))
  colnames(cell_desc) <- c("n", "medianCP", "CP_lwr", "CP_upr", Irvec,
    paste0("pda_", c("median", "lwr", "upr")), paste0("pdb_", c("median", "lwr", "upr")),
    paste0("l_", c("median", "lwr", "upr")), paste0("s_", c("median", "lwr", "upr")))

  rownames(cell_desc) <- model_CP$cells$CellNames
  cell_desc[ , "n"] <- model_CP$cell_ab$n
  # fill in the MLE's as the point estimates
  if (model_CP$distribution != "exponential")
    pda <- model_CP$cell_ab[ , "pda_median"]
    pdb <- model_CP$cell_ab[ , "pdb_median"]
  if (model_CP$distribution == "weibull"){
    cell_desc[ , "medianCP"] <- pdb * log(2)^(1/pda)
  } else if (model_CP$distribution == "lognormal"){
    cell_desc[ , "medianCP"] <- exp(pdb)
  } else if (model_CP$distribution == "loglogistic"){
    cell_desc[ , "medianCP"] <- pdb
  } else if (model_CP$distribution == "exponential"){
    cell_desc[ , "medianCP"] <- log(2) * pdb
  }
  cell_desc[ , Irv] <- t(ppersist(pda = pda, pdb = pdb,
    dist = model_CP$distribution, t_arrive0 = t0, t_arrive1 = t1, t_search = tf))
  # simulate pda, pdb for estimating CIs
  absim <- rcp(n = nsim, model = model_CP, type = "ppersist")
  # calculate r statistics
  rstat <- lapply(absim, function(xx){
    ppersist(xx[,1], xx[,2], dist = model_CP$distribution,
      t_arrive0 = t0, t_arrive1 = t1, t_search = tf)
  })
  ci_lu <- c((1 - CL)/2, 1 - (1 - CL)/2)
  rsum <- lapply(rstat, function(xx)
    matrixStats::rowQuantiles(xx, probs = ci_lu))
  for (ci in rownames(cell_desc)){
    cell_desc[ci, c(paste0(Irv, "_lwr"), paste0(Irv, "_upr"))]  <-
      t(rsum[[ci]])
  }
  # calculate cp statistics
  if (model_CP$distribution == "weibull"){
    mstat <- t(array(unlist(lapply(absim, function(xx)
        quantile(xx[, "pdb"] * log(2)^(1/xx[ , "pda"]), probs = ci_lu))),
        dim = c(2, model_CP$ncell)
    ))
  } else if (model_CP$distribution == "lognormal"){
    mstat <- unlist(exp(model_CP$cell_ab[ , c("pdb_lwr", "pdb_upr")]))
  } else if (model_CP$distribution == "loglogistic"){
    mstat <- unlist(model_CP$cell_ab[ , c("pdb_lwr", "pdb_upr")])
  } else if (model_CP$distribution == "exponential"){
    mstat <- unlist(log(2) * model_CP$cell_ab[ , c("pdb_lwr", "pdb_upr")])
  } else {
    stop("invalid CP distribution")
  }
  # write results into a table
  cell_desc[ , c("CP_lwr", "CP_upr")] <- mstat
  if (model_CP$distribution != "exponential")
    cell_desc[ , paste0("pda_", c("median", "lwr", "upr"))] <-
      as.matrix(model_CP$cell_ab[ , paste0("pda_", c("median", "lwr", "upr"))])

  cell_desc[ , paste0("pdb_", c("median", "lwr", "upr"))] <-
    as.matrix(model_CP$cell_ab[ , paste0("pdb_", c("median", "lwr", "upr"))])

  cell_desc[ , paste0("l_", c("median", "lwr", "upr"))] <-
    as.matrix(model_CP$cell_ls[ , paste0("l_", c("median", "lwr", "upr"))])

  cell_desc[ , paste0("s_", c("median", "lwr", "upr"))] <-
    as.matrix(model_CP$cell_ls[ , paste0("s_", c("median", "lwr", "upr"))])

  class(cell_desc) <- c("descCP", class(cell_desc))
  attr(cell_desc, "CL") <- CL
  attr(cell_desc, "model") <-
    paste0("dist: ", model_CP$distribution, "; ",
    deparse(model_CP$formula_l), "; ", deparse(model_CP$formula_s))
  return(cell_desc)
}