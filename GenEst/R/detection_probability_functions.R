#' @title Estimate all carcass-level detection rates and arrival intervals
#'
#' @description Estimate g values and arrival intervals for a set of carcasses
#'   from fitted pk and cp models and search data
#'
#' @param data_CO Carcass Observation data
#'
#' @param data_SS Search Schedule data
#'
#' @param COdate Column name for the date found data
#'
#' @param model_SE Searcher Efficiency model (or list of models if there are
#'   multiple carcass classes)
#'
#' @param model_CP Carcass Persistence model (or list of models if there are
#'   multiple carcass classes)
#'
#' @param model_DWP Density weighted proportion model (or list of models if
#'  there are multiple carcass classes)
#'
#' @param unitCol Column name for the unit indicator
#'
#' @param IDcol Column name for unique carcass IDs (required)
#'
#' @param SSdate Column name for the date searched data. Optional.
#'   If not provided, \code{estg} will try to find the SSdate among
#'   the columns in data_SS. See \code{\link{prepSS}}.
#'
#' @param sizeCol Name of column in \code{data_CO} where the carcass classes
#'   are recorded. Optional. If not provided, no distinctions are made among
#'   sizes. \code{sizeCol} not only identifies what the name of the size
#    column is, it also identifies that the model should include size as a 
#'   segregating class
#'
#' @param nsim the number of simulation draws
#'
#' @param max_intervals maximum number of arrival interval intervals to 
#'   consider for each carcass. Optional. Limiting the number of search 
#'   intervals can greatly increase the speed of calculations with only a 
#'   slight reduction in accuracy in most cases.
#'
#' @return list of [1] g estimates (\code{ghat}) and [2] arrival interval
#'  estimates (\code{Aj}) for each of the carcasses. The row names of the
#'  \code{Aj} matrix are the units at which carcasses were found. Row names of
#'  \code{ghat} are the carcass IDs (in \code{data_CO}).
#'
#' @examples
#'  data(mock)
#'  model_SE <- pkm(formula_p = p ~ HabitatType, formula_k = k ~ 1,
#'               data = mock$SE)
#'  model_CP <- cpm(formula_l = l ~ Visibility, formula_s = s ~ Visibility, 
#'                data = mock$CP, dist = "weibull",
#'                left = "LastPresentDecimalDays", 
#'                right = "FirstAbsentDecimalDays"
#'              )
#'  ghat <- estg(data_CO = mock$CO, COdate = "DateFound",  data_SS = mock$SS,
#'       model_SE = model_SE, model_CP = model_CP, unitCol = "Unit", nsim = 100)
#'
#' @export
#'
estg <- function(data_CO, COdate, data_SS, SSdate = NULL,
                 model_SE, model_CP, model_DWP = NULL,
                 sizeCol = NULL, unitCol = NULL, IDcol = NULL,
                 nsim = 1000, max_intervals = 8){

  i <- sapply(data_CO, is.factor)
  data_CO[i] <- lapply(data_CO[i], as.character)
  SSdat <- prepSS(data_SS) # SSdat name distinguishes this as pre-formatted
# error-checking
  if (is.null(unitCol))
    unitCol <- defineUnitCol(data_CO = data_CO, data_SS = data_SS)
  if (is.null(IDcol)){
    IDcol <- names(which(
      apply(data_CO, FUN = function(x) length(unique(x)), MARGIN = 2) ==
        apply(data_CO, FUN = length, MARGIN = 2)))
    IDcol <- IDcol[!IDcol %in% COdate]
    if (length(IDcol) == 0){
      stop("CO data must include unique identifier for each caracass")
    } else if (length(IDcol) > 1) {
      warning(paste(
        "No carcass ID column was provided. Taking ", IDcol[1], "as carcass ID"))
      IDcol <- IDcol[1]
#      stop("Carcass ID column must be specified in CO data")
    }
  } else {
    if (length(data_CO[ , IDcol]) != unique(length(data_CO[ , IDcol])))
      stop(paste0("Carcass IDs are not unique (in CO, column ", IDcol, ")"))
  }
  if (any(!data_CO[, unitCol] %in% SSdat$unit))
    stop("carcasses found (CO) at units not properly formatted (or missing) in SS")
  if (is.null(SSdate)) SSdate <- SSdat$SSdate
  SSdat$searches_unit[ , 1] <- 1 # set t0 as start of period of inference
  t0date <- SSdat$date0
  dates_CO <- suppressWarnings(checkDate(data_CO[ , COdate]))
  if (is.null(dates_CO))
    stop("data_CO[ , COdate] values  not properly formatted as dates")
  if (t0date > min(dates_CO))
    stop("first carcass discovered before first search date")
  rind <- match(dates_CO, SSdat[[SSdate]])
  cind <- match(data_CO[, unitCol], SSdat[["unit"]])
  if (anyNA(rind)){
    stop("carcasses ", paste0(which(is.na(rind)), collapse = ', '),
         " in CO found on dates not listed in SS")
  }
  if (anyNA(cind)){
    stop("carcasses discovered at unit(s) ",
      paste0(data_CO[rind, unitCol], collapse =", "),
      "in CO...not represented in SS. Cannot estimate g or M.")
  }
#  if (any(as.numeric(data_SS[cbind(rind, cind + 1)]) == 0)){
  if (any(as.numeric(SSdat$searches_unit[cbind(cind, rind)]) == 0)){
    stop("some carcasses (CO) were found at units that were not searched (SS) ",
         "on the date recorded for the carcass discovery")
  }
  COdat <- data_CO # format data_CO
  COdat[ , COdate] <- dateToDay(dates_CO, t0date)
  names(COdat)[names(COdat) == COdate] <- "day" # distinguish integers
  if (is.null(sizeCol) || is.na(sizeCol)){
    sizeCol <- "placeholder"
    COdat[ , sizeCol] <- "value"
    model_SE <- list("value" = model_SE)
    model_CP <- list("value" = model_CP)
  } else {
    if (!(sizeCol %in% colnames(COdat))){
      stop("carcass class column not in carcass data.")
    }
    if (length(setdiff(names(model_SE), names(model_CP))) > 0) {
      stop("model_SE and model_CP must encompass the same carcass classes")
    }
    if (!all(COdat[, sizeCol] %in% names(model_SE))){
      stop("no SE model for some carcass class represented in data_CO")
    }
  }
  sizeclass <- as.list(as.character(COdat[, sizeCol]))
  sizeclasses <- sort(unique(unlist(sizeclass)))
  nsizeclass <- length(sizeclasses)
# data pre-processing
# create lists of arrays for SS (days) and cells (SE and CP)
  for (sc in sizeclasses){
    if (!("pkm" %in% class(model_SE[[sc]]))){
      stop("Invalid pk model ")
    }
    if (sum(diag(model_SE[[sc]]$varbeta) < 0) > 0){
      stop("
        Cannot estimate variance in user-supplied pk model for carcass class '", sc,
        "' Aborting calculation of ghat."
      )
    }
    if (any(unlist(lapply(model_SE, function(x) x$pOnly)))){
      nok <- names(which(unlist(lapply(model_SE, function(x) x$pOnly))))
      stop("valid k not provided for ", paste(nok, collapse = " "))
    }
  }

  preds_SE <- lapply(model_SE, function(x) x$predictors)
  preds_CP <- lapply(model_CP, function(x) x$predictors)
  preds <- unlist(mapply(function(x, y) unique(c(x, y)), preds_SE, preds_CP))
  if (!all(preds %in% c(names(COdat), names(data_SS)))){
    bad_pr <- unique(preds[!preds %in% c(names(data_CO), names(data_SS))])
    bad_pr <- paste(bad_pr, collapse = " and ")
     stop(paste(bad_pr,
      "included in CP or SE model but not in CO or SS data. ",
      "Cannot assign detection probabilities to carcasses."
    ))
  }
  dist <- unlist(lapply(model_CP, function(x) x$dist))
  COpreds <- lapply(preds, function(x) x[x %in% names(COdat)])
  SSpreds <- lapply(preds, function(x) x[!(x %in% names(COdat))])
  if (length(SSpreds) && max(unlist(lapply(SSpreds, length))) > 1){
    stop("At most 1 SS predictor is allowed per carcass class.")
  }
  if (length(unlist(SSpreds)) > 0 && !all(unlist(SSpreds) %in% names(SSdat))){
    stop("Model predictor missing from both CO and SS data.")
  }
  #############################################
  # error check for matching covariate levels
  # SE
  for (sc in sizeclasses){
    if (length(preds_SE[[sc]]) > 0){
      for (pri in preds_SE[[sc]]){
        if (pri %in% names(data_CO)){
          if (!all(unique(data_CO[ , pri]) %in% model_SE[[sc]]$data[, pri])){
            badind <- which(!(unique(data_CO[ , pri]) %in% model_SE[[sc]]$data[, pri]))
            stop(paste0(
              pri, ' = "', paste(unique(data_CO[ , pri])[badind], collapse = '", '),
              '" found in CO data but not present in SE data. Cannot estimate ',
              "detection probability or mortality."))
          }
        }
      }
    }
    if (length(preds_CP[[sc]]) > 0){
      for (pri in preds_CP[[sc]]){
        if (pri %in% names(data_CO)){
          if (!all(unique(data_CO[ , pri]) %in% model_CP[[sc]]$data[, pri])){
            badind <- which(!(unique(data_CO[ , pri]) %in% model_CP[[sc]]$data[, pri]))
            stop(paste0(
              pri, ' = "', paste(unique(data_CO[ , pri])[badind], collapse = '", '),
              '" found in CO data but not present in CP data. Cannot estimate ',
              "detection probability or mortality."))
          }
        }
      }
    }
  }

  #############################################

  pksim <- lapply(model_SE, function(x) rpk(nsim, x))
  names(pksim) <- names(model_SE)

  cpsim <- lapply(model_CP, rcp, n = nsim, type = "ppersist")

  X <- dim(COdat)[1]
  days <- list()
  for (xi in 1:X){
    toi <- COdat$day[xi]
    SSxi <- SSdat$searches_unit[COdat[xi, unitCol], ] * SSdat$days
    SSxi <- c(0, SSxi[SSxi > 0])
    days[[xi]] <- SSxi[SSxi <= toi]
    if (!is.null(max_intervals)){
      dlen <- length(days[[xi]])
      if (dlen > max_intervals + 1)
        days[[xi]] <- days[[xi]][(dlen - max_intervals):dlen]
    }
  }

  cells <- list()
  # take care of the SS preds
  if (sum(unlist(lapply(SSpreds, length))) == 0){
    for (xi in 1:X){
      cells[[xi]] <- list()
      cells[[xi]][[sizeCol]] <- COdat[xi, sizeCol]
      pcol <- preds_SE[[COdat[xi, sizeCol]]]
      if (length(pcol) == 0) {
        cells[[xi]]$SEcell <- "all"
      } else {
        cells[[xi]]$SEcell <- paste(COdat[xi, pcol], collapse = ".")
      }
      cells[[xi]]$SErep <- length(days[[xi]]) - 1
      pcol <- preds_CP[[COdat[xi, sizeCol]]]
      if (length(pcol) == 0) {
        cells[[xi]]$CPcell <- "all"
      } else {
        cells[[xi]]$CPcell <- paste(COdat[xi, pcol], collapse = ".")
      }
      cells[[xi]]$CPrep <- length(days[[xi]]) - 1
    }

  } else {

    for (xi in 1:X){
      cells[[xi]] <- list()
      SEc <- SEr <- CPc <- CPr <- NULL
      sz <- as.character(COdat[xi, sizeCol])
      cells[[xi]][[sizeCol]] <- sz
      # interpret the SE predictors
      nse <- length(preds_SE[[sz]])
      if (nse == 0){
        cells[[xi]]$SEcell <- "all"
        cells[[xi]]$SErep <- length(days[[xi]]) - 1
      } else {
        for (sei in 1:nse){
          predi <- preds_SE[[sz]][sei]
          if (predi %in% SSpreds){
            ssvec <- (SSdat[[predi]][which(SSdat[["days"]] %in% days[[xi]])])
            ssvec <- ssvec[-1]
            SEc <- paste(SEc, unique(ssvec),
              sep = ifelse(is.null(SEc), "", "."))
            SEr <- table(ssvec)[unique(ssvec)]
          } else {
            SEc <- paste(SEc, COdat[xi, predi],
              sep = ifelse(is.null(SEc), "", "."))
          }
        }
        cells[[xi]]$SEcell <- SEc
        if (is.null(SEr)){
          cells[[xi]]$SErep <- length(days[[xi]]) - 1
        } else {
          cells[[xi]]$SErep <- SEr
        }
      }
      # interpret the CP predictors
      nse <- length(preds_CP[[sz]])
      if (nse == 0){
        cells[[xi]]$CPcell <- "all"
        cells[[xi]]$CPrep <- length(days[[xi]]) - 1
      } else {
        for (sei in 1:nse){
          predi <- preds_CP[[sz]][sei]
          if (predi %in% SSpreds){
            ssvec <- (SSdat[[predi]][which(SSdat[["days"]] %in% days[[xi]])])
            ssvec <- ssvec[-1]
            CPc <- paste(CPc, unique(ssvec),
              sep = ifelse(is.null(CPc), "", "."))
            CPr <- table(ssvec)[unique(ssvec)]
          } else {
            CPc <- paste(CPc, COdat[xi, predi],
              sep = ifelse(is.null(CPc), "", "."))
          }
        }
        cells[[xi]]$CPcell <- CPc
        if (is.null(CPr)){
          cells[[xi]]$CPrep <- length(days[[xi]]) - 1
        } else {
          cells[[xi]]$CPrep <- CPr
        }
      }
    }
  }
  # the calculation
  ghat <- matrix(0, nrow = X, ncol = nsim)
  Aj <- matrix(0, nrow = X, ncol = nsim)
  for (xi in 1:X){
    if (COdat$day[xi] == 0) next # cleanout: leaves initial 0s in ghat and Aj
    SSxi <- SSdat$searches_unit[COdat[xi, unitCol], ] * SSdat$days
    SSxi <- c(0, SSxi[SSxi > 0])
    # calculate SE
    sz <- cells[[xi]][[sizeCol]]
    SEr <- cells[[xi]]$SErep
    oi <- length(days[[xi]]) - 1
    rng <- 0
    pOigAj <- NULL # virtually identical calcs done for pAjgOi, but in reverse
    for (sei in 1:length(SEr)){
      rng <- max(rng) + 1:SEr[sei]
      pOigAj <- cbind(pOigAj, SEsi_left(
        oi = oi,
        pk = pksim[[sz]][[cells[[xi]]$SEcell[sei]]],
        rng = rng
      ))
    }
    # multiply by ppersist
    CPr <- cells[[xi]]$CPrep
    rng <- 0
    for (cpi in 1:length(CPr)){
      rng <- max(rng) + 1:CPr[cpi]
      pOigAj[, rng] <- pOigAj[, rng] * t(ppersist(
        pda = cpsim[[sz]][[cells[[xi]]$CPcell[cpi]]][ , "pda"],
        pdb = cpsim[[sz]][[cells[[xi]]$CPcell[cpi]]][ , "pdb"],
        dist = model_CP[[sz]]$distribution,
        t_arrive0 = days[[xi]][rng],
        t_arrive1 = days[[xi]][rng + 1],
        t_search = rep(max(days[[xi]]), length(rng))
      ))
    }

    parrive <- diff(days[[xi]][1:(oi+1)])/days[[xi]][oi+1]
    pAjgOi <- t(pOigAj) * parrive; pAjgOi <- t(t(pAjgOi)/colSums(pAjgOi))
    Aj[xi, ] <- # sim arrival intervals (relative to cind's ss)
       rowSums(matrixStats::rowCumsums(t(pAjgOi)) < runif(nsim)) +
         (sum(SSxi <= min(days[[xi]])))
    xuint <- unique(Aj[xi, ]) # unique xi arrival intervals (in SSxi)
    for (aj in xuint){
      # calculate simulated ghat associated with the given carcass and 
      #   interval (there is much redundant calculation here that could be sped
      #   up substantially with clever bookkeeping)
      simInd <- which(Aj[xi, ] == aj)
      top <- length(SSxi)
      if (!is.null(max_intervals)){
        # the calculations on RHS are more critical and less time consuming
        # calculation-wise, so for now...
        #top <- min(aj + max_intervals, top)
      }
      # use an adjusted search schedule because we "know" when carcass arrived
      # which cell is "active" for the given arrival interval?
      cpi <- findInterval(aj, c(0, min(xuint) + cumsum(cells[[xi]]$CPrep)),
        rightmost.closed = T)
      pda <- cpsim[[sz]][[cells[[xi]]$CPcell[cpi]]][simInd , "pda"]
      pdb <- cpsim[[sz]][[cells[[xi]]$CPcell[cpi]]][simInd , "pdb"]
      ppersu <- ppersist(
        pda = pda,
        pdb = pdb,
        dist = model_CP[[sz]]$distribution,
        t_arrive0 = rep(SSxi[aj], top - aj),
        t_arrive1 = rep(SSxi[aj + 1], top - aj),
        t_search = SSxi[(aj + 1):top]
      )
      pki <- findInterval(aj, c(0, min(xuint) + cumsum(cells[[xi]]$SErep)),
        rightmost.closed = T)
      SE <- t(SEsi_right(
        top - aj,
        pksim[[sz]][[cells[[xi]]$SEcell[pki]]][simInd , ]
      ))
      if (aj < top - 1){
        ghat[xi, simInd] <- colSums(SE * ppersu)
      } else {
        ghat[xi, simInd] <- as.vector(SE) * as.vector(ppersu)
      }
    }
  }
  if (is.null(model_DWP)) {
    DWP <- 1
  } else {
    dwpsim <- rdwp(n = nsim, model = model_DWP)
    if (length(dwpsim) == 1){
      DWP <- 1
    } else {
      DWP <- CO_DWP(dwpsim = dwpsim, data_CO = data_CO, unitCol = unitCol, sizeCol = sizeCol)
    }
  }
  rownames(Aj) <- COdat[ , unitCol]
  ghat[ghat < 1e-6 & ghat > 1e-15] <- 1e-6 # prevents 1/0 in estM
  rownames(ghat) <- COdat[ , IDcol]
  out <- list("ghat" = ghat, "Aj" = Aj, "DWP" = DWP) # ordered by relevance to user
  return(out)
}

#' @title Associate CO carcasses with appropriate DWP values (by unit and carcass class)
#'
#' @description Calculate the conditional probability of observing a carcass
#'   at search oi as a function arrival interval (assuming carcass is not
#'   removed by scavengers before the time of the final search)
#'
#' @param dwpsim \code{rdwp} object
#'
#' @param data_CO data frame with results from carcass surveys
#'
#' @param unitCol name of the unit column in data_CO (required)
#'
#' @param sizeCol name of the carcass class column in data_CO (optional).
#'
#' @return numeric DWP array
#'
#' @export
#'
CO_DWP <- function(dwpsim, data_CO, unitCol, sizeCol = NULL){
  if (!"rdwp" %in% class(dwpsim))
    stop("dwpsim must be of class 'rdwp'")
  if (!unitCol %in% names(data_CO))
    stop("unitCol must be the name of a valid unit column in data_CO")
  if (length(unitCol) > 1)
    stop("unitCol must be the name of a unique, valid unit column in data_CO")

  if (is.null(sizeCol) || sizeCol == "placeholder"){
    if (is.list(dwpsim))
      stop("dwpsim should be an array rather than a list when no sizeCol is provided")
    if (!all(data_CO[ , unitCol] %in% row.names(dwpsim)))
      stop("some units in data_CO not represented in dwpsim")
    DWP <- dwpsim[match(data_CO[, unitCol], row.names(dwpsim)), ]
  } else {
    if (!sizeCol %in% names(data_CO))
      stop("sizeCol not in data_CO")
    if (!is.list(dwpsim) || !all(data_CO[ , sizeCol] %in% names(dwpsim)))
      stop("dwpsim must be a list to match sizes in data_CO[ , sizeCol]")
    DWP <- array(dim = c(nrow(data_CO), ncol(dwpsim[[1]])))
    for (ci in 1:nrow(data_CO)){
      DWP[ci, ] <- dwpsim[[data_CO[ci, sizeCol]]][data_CO[ci, unitCol], ]
    }
  }
  if (NCOL(DWP) == 1) DWP <- as.vector(DWP)
  return(DWP)
}
#' @title Calculate conditional probability of observation at a search
#'
#' @description Calculate the conditional probability of observing a carcass 
#'   at search oi as a function arrival interval (assuming carcass is not
#'   removed by scavengers before the time of the final search)
#'   
#' @param oi number of searches after arrival
#'
#' @param pk numeric array of searcher efficiency p and k parameters
#'  (p = pk[ , 1] and k = pk[ , 2])
#'
#' @param rng optional parameter giving the range of intervals to consider
#'
#' @return numeric array of probability of observing a carcass at oi for
#'   given that it arrived in intervals 1:oi if rng = NULL (or in intervals
#'   \code{rng}), assuming the carcass had not been previously discovered or
#'   removed by scavengers
#'
#' @export
#'
SEsi_left <- function (oi, pk, rng = NULL){
  # oi is the index for the search occasion (excluding t0)
  # pk is nsim x 2 array of simulated p and k parameters
  # rng is the intervals for which to calculate answer
  if (is.null(rng)) rng <- 1:oi
  if (is.null(dim(pk)) || nrow(pk) == 1) return(SEsi0(0:oi, pk))
  npk <- nrow(pk)
  nmiss <- oi - rng
  maxmiss <- max(nmiss)
  if (maxmiss == 0){
    pfind.si <- matrix(pk[, 1], ncol = 1)
  }
  else if (maxmiss == 1) {
    pfind.si <- cbind(pk[, 1], (1 - pk[, 1]) * pk[, 2] * pk[, 1])
  }
  else {
    powk <- array(rep(pk[, 2], maxmiss + 1), dim = c(npk, maxmiss + 1))
    powk[, 1] <- 1
    powk <- matrixStats::rowCumprods(powk)
    pfind.si <- pk[, 1] * powk * cbind(
      rep(1, npk),
      matrixStats::rowCumprods(1 - (pk[, 1] * powk[, 1:maxmiss]))
    )
  }
  return(pfind.si[ , oi - rng + 1])
}

#' @title Calculate conditional probability of observation after a series of 
#'   searches
#'
#' @description Calculate the conditional probability of observing a carcass 
#'   after i = 1:nsi searches (assuming carcass is not previous discovered by 
#'   searchers or removed by scavengers)
#'
#' @param nsi number of searches after arrival
#'
#' @param pk numeric array of searcher efficiency p and k parameters
#'  (p = pk[ , 1] and k = pk[ , 2])
#'
#' @return numeric nsi x dim(pk)[1] array of probabilities of observing a
#'  carcass after 1:nsi searches (assuming that the carcass had not been
#' previously discovered or removed by scavengers
#'
#' @export
#'
SEsi_right <- function(nsi, pk){
  # oi is the index for the search occasion (excluding t0)
  # pk is nsim x 2 array of simulated p and k parameters
  # rng is the intervals for which to calculate answer
  if (is.null(dim(pk)) || nrow(pk) == 1) return(t(SEsi0(0:nsi, pk)))
  npk <- nrow(pk)
  nmiss <- 1:nsi - 1
  maxmiss <- max(nmiss)
  if (maxmiss == 0){
    pfind.si <- matrix(pk[, 1], ncol = 1)
  }
  else if (maxmiss == 1) {
    pfind.si <- cbind(pk[, 1], (1 - pk[, 1]) * pk[, 2] * pk[, 1])
  }
  else {
    powk <- array(rep(pk[, 2], maxmiss + 1), dim = c(npk, maxmiss + 1))
    powk[, 1] <- 1
    powk <- matrixStats::rowCumprods(powk)
    pfind.si <- pk[, 1] * powk * cbind(
      rep(1, npk),
      matrixStats::rowCumprods(1 - (pk[, 1] * powk[, 1:maxmiss]))
    )
  }
  return(pfind.si)
}

#' @title Estimate generic g
#'
#' @description Generic g estimation by simulation from given SE model and CP 
#'   models under a specific search schedule.
#'
#' The g estimated by \code{estgGeneric} is a generic aggregate detection 
#'   probability and represents the probability of detecting a carcass that 
#'   arrives at a (uniform) random time during the period monitored, for each
#'   of the possible cell combinations, given the SE and CP models. This 
#'   is somethat different from the GenEst estimation of g when the purpose 
#'   is to estimate total mortality (M), in which case the detection 
#'   probability varies with carcass arrival interval and is difficult to 
#'   summarize statistically. The \code{estgGeneric} estimate is a useful 
#'   "big picture" summary of detection probability, but would be difficult
#'   to work with for estimating M with precision.
#'
#' @param nsim the number of simulation draws
#'
#' @param days Search schedule data as a vector of days searched
#'
#' @param model_SE Searcher Efficiency model (\code{pkm} object)
#'
#' @param model_CP Carcass Persistence model (\code{cpm} object)
#'
#' @return \code{gGeneric} object that is a list of [1] a list of g estimates,
#'    with one element in the list corresponding to each of the cells from the
#'    cross-model combination and [2] a table of predictors and cell names 
#'    associated with the gs
#'
#' @examples
#'   data(mock)
#'   model_SE <- pkm(formula_p = p ~ HabitatType, formula_k = k ~ 1,
#'                 data = mock$SE)
#'   model_CP <- cpm(formula_l = l ~ Visibility, formula_s = s ~ Visibility, 
#'                 data = mock$CP, left = "LastPresentDecimalDays", 
#'                 right = "FirstAbsentDecimalDays")
#'   avgSS <- averageSS(mock$SS)
#'   ghatsGeneric <- estgGeneric(days = avgSS, model_SE = model_SE,
#'    model_CP = model_CP)
#'
#' @export
#'
estgGeneric <- function(days, model_SE, model_CP, nsim = 1000){

  if (!is.vector(days) || !is.numeric(days))
    stop(" 'days' must be a numeric vector")
  if (!("pkm" %in% class(model_SE))) stop("Invalid pk model")
  if (anyNA(diag(model_SE$varbeta)) || sum(diag(model_SE$varbeta) < 0) > 0)
    stop("Cannot estimate variance for model_SE. Aborting estimation.")
  preds_SE <- model_SE$predictors
  preds_CP <- model_CP$predictors
  data_SE <- data.frame(model_SE$data[ , preds_SE], stringsAsFactors = FALSE)
  data_CP <- data.frame(model_CP$data[ , preds_CP], stringsAsFactors = FALSE)
  names(data_SE) <- preds_SE
  names(data_CP) <- preds_CP
  preds <- combinePredsAcrossModels(preds_CP, preds_SE, data_CP, data_SE)
  sim_SE <- rpk(n = nsim, model = model_SE)
  sim_CP <- rcp(n = nsim, model = model_CP, type = "ppersist")

  ncell <- nrow(preds)
  ghat <- vector("list", ncell)
  for (celli in 1:ncell){
    cell_SE <- preds$CellNames_SE[celli]
    cell_CP <- preds$CellNames_CP[celli]
    param_SE <- sim_SE[[cell_SE]]
    param_CP <- sim_CP[[cell_CP]]
    ghat[[celli]] <- calcg(days, param_SE, param_CP, model_CP$dist)
  }  
  names(ghat) <- preds$CellNames
  span <- max(days)
  I <- span/(length(days) - 1)
  out <- list("ghat" = ghat, "predictors" = preds, "SS" = list(span = span, I = I))
  class(out) <- c("gGeneric", "list")
  return(out)
}

#' @title Calculate detection probability for given SE and CP parameters and
#'  search schedule.
#'
#' @description Calculate detection probability (g) given SE and CP parameters
#'  and a search schedule.
#'
#' The g given by \code{calcg} is a generic aggregate detection
#'  probability and represents the probability of detecting a carcass that
#'  arrives at a (uniform) random time during the time spanned by the search
#'  schedule for the the given SE and CP parameters. This differs from the GenEst
#'  estimation of g when the purpose is to estimate total mortality (M), in
#'  which case the detection probability varies with carcass arrival interval
#'  and is difficult to summarize statistically. \code{calcg} provides a 
#'  useful "big picture" summary of detection probability, but would be 
#'  difficult to work with for estimating M with precision.
#'
#' @param days Search schedule (vector of days searched)
#'
#' @param param_SE numeric array of searcher efficiency parameters (p and k);
#'  must have the name number of rows as the \code{param_CP}.
#'
#' @param param_CP numeric array of carcass persistence parameters (a and b)
#'  must have the name number of rows as the \code{param_SE}.
#'
#' @param dist distribution for the CP model
#'
#' @return a vector of detection probabilities for each
#'
#' @export
#'
calcg <- function(days, param_SE, param_CP, dist){
  dist <- tolower(dist)
  samtype <- ifelse(length(unique(diff(days))) == 1, "Formula", "Custom")
  nsearch <- length(days) - 1
  # check format of param_SE and create formatted matrix pk
  # should be rectangle shape with columns p and k (named or not)
  if (any(c("array", "matrix", "data.frame") %in% class(param_SE))){
    if (ncol(param_SE) < 2){
      stop("param_SE must have columns for p and k")
    }
    if (!("p" %in% colnames(param_SE) & "k" %in% colnames(param_SE))){
      if (ncol(param_SE) > 2 | (ncol(param_SE) == 2 & is.null(colnames(param_SE))))
        stop("param_SE must have columns for p and k")
      # param_SE now either has two or more columns, which include p and k,
      # or param_SE is two-column structure with no column names
    }
    if (is.null(colnames(param_SE))){ # two-column structure
      colnames(param_SE) <- c("p", "k")
    }
    pk <- as.matrix(param_SE)
  } else { # should be a 2-vector with p and k (named or not)
    if (!is.numeric(param_SE) || !is.vector(param_SE))
      stop("param_SE must be numeric vector or array")
    if (length(param_SE) == 2){
      if (!is.null(names(param_SE)) & !identical(sort(names(param_SE)), c("k", "p")))
        stop("pamam_SE must have p and k values")
      if (!is.null(names(param_SE))){
        pknm <- names(param_SE)
      } else {
        pknm <- c("p", "k")
      }
      pk <- matrix(param_SE, nrow = 1)
      colnames(pk) <- pknm
    } else if (length(param_SE) == 1) {
      stop("param_SE must include values for both p and k")
    } else {
      stop("param_SE must be either 2-d array with two colunms ",
           "or a vector of length 2")
    }
  }
#  pk <- pk[, 1:2]
  n <- nrow(pk)
  # check format of param_CP and create formatted matrix cp
  if (is.vector(param_CP)){
    if (!is.numeric(param_CP)) stop("param_CP must be numeric")
    # must be of length 2 and n = 1
    # or exponential and length = n     (pdb)
    if (length(param_CP) == 2 & n == 1){
      if (!is.null(names(param_CP))){
        if (!identical(sort(names(param_CP)), c("pda", "pdb"))){
          stop("param_CP must contain values for pda and pdb ",
               "(either a vector of length 2 or a matrix with ",
               "columns for pda and pdb)"
          )
        } else {
          cpnm <- names(param_CP)
        }
      } else {
        cpnm <- c("pda", "pdb")
      }
      cp <- matrix(param_CP, nrow = 1)
      colnames(cp) <- cpnm
    } else if (dist == "exponential" & length(param_CP) == n & n > 1){
      # assumed that param_CP is pdb
      cp <- as.matrix(cbind(1/param_CP, param_CP))
      colnames(cp) <- c("pda", "pdb")
    } else if (dist == "exponential" & length(param_CP) == n & n == 1){
      # check whether param_CP is pda or pdb
      if (is.null(names(param_CP))){
        cp <- matrix(c(1/param_CP, param_CP), nrow = 1)
        colnames(cp) <- c("pda", "pdb")
      } else if (!names(param_CP) %in% c("pda", "pdb")){
          stop("param_CP with exponential distribution must be pda or pdb")
      } else {
        if (names(param_CP) == "pda"){
          cp <- matrix(c(param_CP, 1/param_CP), nrow = 1)
        } else {
          cp <- matrix(c(1/param_CP, param_CP, nrow = 1))
        }
        colnames(cp) <- c("pda", "pdb")
       }
      } else {
        stop("param_SE and param_CP have incompatible dimensions")
      }
  } else {
    if (!any(class(param_CP) %in% c("array", "matrix", "data.frame")))
      stop ("param_CP must be a 2-d array with columns for pda and pdb")
    cp <- as.matrix(param_CP)
    if (ncol(cp) == 1){
      if (dist != "exponential")
        stop ("param_CP must have columns for pda and pdb if dist != 'exponential'")
      if (is.null(colnames(cp))){
        cp <- cbind(1/cp, cp)
        colnames(cp) <- c("pda", "pdb")
      } else if (colnames(cp) == "pda"){
        cp <- cbind(cp, 1/cp)
        colnames(cp) <- c("pda", "pdb")
      } else if (colnames(cp) == "pdb"){
        cp <- cbind(1/cp, cp)
        colnames(cp) <- c("pda", "pdb")
      } else {
        stop("exponential param_CP must be pda or pdb (or unnamed)")
      }
    } else if (ncol(cp) == 2) {
      if (is.null(colnames(cp))){
        colnames(cp) <- c("pda", "pdb")
      } else if (!identical(sort(colnames(cp)), c("pda", "pdb"))){
        stop("param_CP must have columns pda and pdb")
      }
    } else if (ncol(cp) > 2){
      if (is.null(colnames(cp)) ||
          !"pda" %in% colnames(cp) | !"pdb" %in% colnames(cp)){
        stop("param_CP must have columns pda and pdb")
      }
      cp <- as.matrix(param_CP[, c("pda", "pdb")])
    }

  }
  # check for compatibility between cp and pk parameters
#  if (!identical(dim(pk), dim(cp)))
  if (nrow(pk) != nrow(cp))
    stop("param_SE and param_CP must be the same dimension")
 if (dist == "exponential"){
    pdb0 <- exp(mean(log(cp[ , "pdb"])))
    pda0 <- 1/pdb0
  } else {
    if (dist == "lognormal"){
      pdb0 <- mean(cp[ , "pdb"])
    } else {
      pdb0 <- exp(mean(log(cp[ , "pdb"])))
    }
    pda0 <- 1/mean(1/cp[ , "pda"])
  }
  if (length(days) == 2){ # then only one search, so it's easy!
    r_sim <- as.vector(GenEst::ppersist(dist = dist,
      t_arrive0 = days[1], t_arrive1 = days[2], t_search = days[2],
      pda = cp[ , "pda"], pdb = cp[ , "pdb"]))
    return(r_sim * pk[ , "p"])
  }


  ind1 <- rep(1:nsearch, times = nsearch:1)
  ind2 <- ind1 + 1
  ind3 <- unlist(lapply(1:nsearch, function(x) x:nsearch)) + 1
  schedule.index <- cbind(ind1, ind2, ind3)
  schedule <- cbind(days[ind1], days[ind2], days[ind3])
  nmiss <- schedule.index[,3] - schedule.index[,2]
  maxmiss <- max(nmiss)

  if (maxmiss == 0) {
    pfind.si <- pk[ , "p"]
  } else if (maxmiss == 1){
    pfind.si <- cbind(pk[ , "p"], (1 - pk[ , "p"]) * pk[ , "k"] * pk[ , "p"])
  } else {
    powk <- array(rep(pk[ , "k"], maxmiss + 1), dim = c(n, maxmiss + 1))
    powk[ , 1] <- 1
    powk <- matrixStats::rowCumprods(powk)
    val <- 1 - (pk[ , "p"] * powk[ , 1:maxmiss])
    if (is.null(dim(val))) val <- matrix(val, nrow = 1)
    pfind.si <- pk[ , "p"] * powk * cbind(rep(1, n), matrixStats::rowCumprods(val))
  }
  diffs <- cbind(schedule[,2] - schedule[,1], schedule[,3] - schedule[,2])
  intxsearch <- unique(diffs, MARGIN = 1)
  ppersu <- GenEst::ppersist(dist = dist, t_arrive0 = 0, t_arrive1 = intxsearch[ , 1],
              t_search = intxsearch[ , 1] + intxsearch[ , 2],
              pda = cp[ , "pda"], pdb = cp[ , "pdb"]
            )
  arrvec <- (schedule[ , 2] - schedule[ , 1]) / max(days)
  prob_obs <- numeric(n)
  if (maxmiss > 0){
    for (i in 1:dim(schedule)[1]){
      ind <- which(
               abs(intxsearch[,1] - (schedule[i,2] - schedule[i,1])) < 0.001 &
               abs(intxsearch[,2] - (schedule[i,3] - schedule[i,2])) < 0.001)
      prob_obs <- prob_obs +
                  pfind.si[ , nmiss[i] + 1] * ppersu[ind, ] * arrvec[i]
    }
  } else {
    for (i in 1:dim(schedule)[1]){
      ind <- which(
               abs(intxsearch[,1] - (schedule[i,2] - schedule[i,1])) < 0.001 &
               abs(intxsearch[,2] - (schedule[i,3] - schedule[i,2])) < 0.001)
     prob_obs <- prob_obs + pfind.si[nmiss[i] + 1] * ppersu[ind, ] * arrvec[i]
    }
  }
  names(prob_obs) <- NULL
  return(prob_obs)
}
#' @title Estimate generic detection probability for multiple carcass classes
#'
#' @description Generic g estimation for a combination of SE model and CP
#'   model under a given search schedule
#'
#' The g estimated by \code{estgGenericSize} is a generic aggregate detection
#'   probability and represents the probability of detecting a carcass that 
#'   arrives at a (uniform) random time during the period monitored, for each
#'   of the possible cell combinations, given the SE and CP models. This 
#'   is somethat different from the GenEst estimation of g when the purpose 
#'   is to estimate total mortality (M), in which case the detection 
#'   probability varies with carcass arrival interval and is difficult to 
#'   summarize statistically. The \code{estgGeneric} estimate is a useful 
#'   "big picture" summary of detection probability, but would be difficult
#'   to work with for estimating M with precision.
#'
#' @param nsim the number of simulation draws
#'
#' @param days Search schedule data as a vector of days searched 
#'
#' @param modelSetSize_SE Searcher Efficiency model set for multiple sizes
#'
#' @param modelSetSize_CP Carcass Persistence model set for multiple sizes
#'
#' @param modelSizeSelections_SE vector of SE models to use, one for each 
#'  size. Size names are required, and names must match those of
#'  modelSetSize_SE. E.g., 
#'  \code{c(lrg = "p ~ Visibility; k ~ 1", sml = "p ~ 1; k ~ 1")}.
#'  Model formulas are read as text and must have exact matches among models
#'  listed in modelSetSize_SE. For example, if one of the
#'  \code{modelSizeSelections_SE} elements is
#'  \code{lrg = "p ~ Visibility; k ~ 1"}, then \code{"p ~ Visibility; k ~ 1"}
#'  must be in \code{names(modelSizeSelections_SE)[["lrg"]]}.
#'
#' @param modelSizeSelections_CP vector of CP models to use, one for each size
#'
#' @return list of g estimates, with one element in the list corresponding
#'    to each of the cells from the cross-model combination
#'
#' @examples
#'   data(mock)
#'   pkmModsSize <- pkm(formula_p = p ~ HabitatType,
#'                    formula_k = k ~ HabitatType, data = mock$SE,
#'                    obsCol = c("Search1", "Search2", "Search3", "Search4"),
#'                    sizeCol = "Size", allCombos = TRUE)
#'   cpmModsSize <- cpm(formula_l = l ~ Visibility,
#'                    formula_s = s ~ Visibility, data = mock$CP,
#'                    left = "LastPresentDecimalDays",
#'                    right = "FirstAbsentDecimalDays",
#'                    dist = c("exponential", "lognormal"),
#'                    sizeCol = "Size", allCombos = TRUE)
#'
#'   pkMods <- c("S" = "p ~ 1; k ~ 1", "L" = "p ~ 1; k ~ 1",
#'              "M" = "p ~ 1; k ~ 1", "XL" = "p ~ 1; k ~ 1"
#'             )
#'   cpMods <- c("S" = "dist: exponential; l ~ 1; NULL", 
#'               "L" = "dist: exponential; l ~ 1; NULL",
#'               "M" = "dist: exponential; l ~ 1; NULL",
#'               "XL" = "dist: exponential; l ~ 1; NULL"
#'             )
#'   avgSS <- averageSS(mock$SS)
#'   gsGeneric <- estgGenericSize(nsim = 1000, days = avgSS,
#'                  modelSetSize_SE = pkmModsSize,
#'                  modelSetSize_CP = cpmModsSize,
#'                  modelSizeSelections_SE = pkMods,
#'                  modelSizeSelections_CP = cpMods
#'                )
#'
#' @export
#'
estgGenericSize <- function(days, modelSetSize_SE, modelSetSize_CP,
    modelSizeSelections_SE, modelSizeSelections_CP,
    nsim = 1000){

  if (!("pkmSetSize" %in% class(modelSetSize_SE))){
    stop("modelSetSize_SE must be a pkmSetSize object")
  }
  if (!("cpmSetSize" %in% class(modelSetSize_CP))){
    stop("modelSetSize_CP must be a cpmSetSize object")
  }
  sizeclasses_SE <- names(modelSetSize_SE)
  sizeclasses_CP <- names(modelSetSize_CP)
  if (!all(sizeclasses_SE %in% sizeclasses_CP) ||
      !all(sizeclasses_CP %in% sizeclasses_SE)){
    stop("Carcass classes don't match between SE and CP model sets")
  }
  sizeclasses <- sort(unique(c(sizeclasses_SE, sizeclasses_CP)))
  nsizeclass <- length(sizeclasses)
  # check whether k is included in every model. If not, error.
  for (sci in sizeclasses){
    if (modelSetSize_SE[[sci]][[modelSizeSelections_SE[sci]]]$pOnly){
     stop("k required for SE model for carcass class = ", sci)
    }
  }
  ghats <- list()
  for (sci in sizeclasses){
    if (any(unlist(lapply(modelSetSize_SE[[sci]], function(x) x$pOnly))))
      stop("No k included in SE model. Cannot estimate g")
    model_SE <- modelSetSize_SE[[sci]][[modelSizeSelections_SE[[sci]]]]
    model_CP <- modelSetSize_CP[[sci]][[modelSizeSelections_CP[[sci]]]]
    ghats[[sci]] <- estgGeneric(nsim = nsim, days = days,
      model_SE = model_SE, model_CP = model_CP)
  }

  class(ghats) <- c("gGenericSize", "list")
  return(ghats)
}

#' @title Tabulate an average search schedule from a multi-unit SS data table
#'
#' @description Given a multi-unit Search Schedule data table, produce an 
#'   average search schedule for use in generic detection probability 
#'   estimation. 
#'
#' @param data_SS a multi-unit SS data table, for which the average interval 
#'   will be tabulated. It is assumed that \code{data_SS} is properly 
#'   formatted, with a column of search dates and a column of 1s and 0s for 
#'   each unit indicating whether the unit was searched on the given date).
#'   Other columns are optional, but optional columns should not all contain
#'   at least on value that is not a 1 or 0.
#'
#' @param SSdate Column name for the date searched data (optional).
#'   if no \code{SSdate} is provided, \code{data_SS} will be parsed
#'   to extract the dates automatically. If there is more than one column with
#'   dates, then an error will be thrown and the user will be required to
#'   provide the name of the desired dates column.
#'
#' @return vector of the average search schedule
#'
#' @examples 
#'   data(mock)
#'   avgSS <- averageSS(mock$SS)
#'
#' @export
#'
averageSS <- function(data_SS, SSdate = NULL){
  SSdat <- prepSS(data_SS, SSdate = SSdate)
  schedules <- t(SSdat$searches_unit) * SSdat$days
  nintervals <- length(SSdat$days) - matrixStats::colCounts(schedules, value = 0)
  maxdays <- matrixStats::colMaxs(schedules)
  aveSS <- seq(0, max(maxdays), round(mean(maxdays/nintervals)))
  return(aveSS)
}
  
#' @title Summarize the gGeneric list to a simple table
#'
#' @description methods for \code{summary} applied to a \code{gGeneric} list
#'
#' @param object gGeneric output list (each element is a named vector of 
#'   gGeneric values for a cell in the model combinations)
#'
#' @param ... arguments to be passed down
#'
#' @param CL confidence level
#'
#' @return a summary table of g values (medians and confidence bounds) for 
#'   each cell combination within the gGeneric list
#'
#' @examples 
#'   data(mock)
#'   model_SE <- pkm(formula_p = p ~ HabitatType, formula_k = k ~ 1,
#'                 data = mock$SE)
#'   model_CP <- cpm(formula_l = l ~ Visibility, formula_s = s ~ Visibility, 
#'                 data = mock$CP, left = "LastPresentDecimalDays", 
#'                 right = "FirstAbsentDecimalDays")
#'   avgSS <- averageSS(mock$SS)
#'   ghatsGeneric <- estgGeneric(nsim = 1000, avgSS, model_SE, model_CP)
#'   summary(ghatsGeneric)
#'
#' @export
#'
summary.gGeneric <- function(object, ..., CL = 0.90){
  ghats <- object$ghat
  preds <- object$predictors
  cells <- names(ghats)
  ncell <- length(cells)
  predsByCell <- strsplit(cells, "\\.") 
  npred <- length(predsByCell[[1]])

  predsTab <- preds[ , -grep("CellNames", colnames(preds))]
  predsTab <- as.matrix(predsTab, ncol = npred, nrow = ncell)
  predNames <- colnames(preds)[-grep("CellNames", colnames(preds))]
  if (length(predNames) == 1 & predNames[1] == "group" & cells[1] == "all"){
    predNames <- "Group"
  }

  tableProbs <- c((1 - CL) / 2, 0.25, 0.5, 0.75, 1 - (1 - CL) / 2)

  colnames(predsTab) <- predNames
  gTab <- matrix(NA, ncell, 5)
  for (celli in 1:ncell){
    gspec <- ghats[[celli]]
    quants <- quantile(gspec, prob = tableProbs)
    gTab[celli, ] <- round(quants, 3)
  }

  out <- data.frame(predsTab, gTab)
  colnames(out)[npred + (1:5)] <- names(quants)
  return(out)
}

#' @title Summarize the gGenericSize list to a list of simple tables
#'
#' @description methods for \code{summary} applied to a \code{gGenericSize}
#'   list
#'
#' @param object gGenericSize output list (each element is a size-named 
#'   list of named vectors of gGeneric values for a cell in the model 
#'   combinations)
#'
#' @param ... arguments to be passed down
#'
#' @param CL confidence level
#'
#' @return a list of summary tables of g values (medians and confidence 
#'   bounds) for each cell combination within the gGeneric list
#'
#' @examples
#'   data(mock)
#'   pkmModsSize <- pkm(formula_p = p ~ HabitatType,
#'                    formula_k = k ~ HabitatType, data = mock$SE,
#'                    obsCol = c("Search1", "Search2", "Search3", "Search4"),
#'                    sizeCol = "Size", allCombos = TRUE)
#'   cpmModsSize <- cpm(formula_l = l ~ Visibility,
#'                    formula_s = s ~ Visibility, data = mock$CP,
#'                    left = "LastPresentDecimalDays",
#'                    right = "FirstAbsentDecimalDays",
#'                    dist = c("exponential", "lognormal"),
#'                    sizeCol = "Size", allCombos = TRUE)
#'   pkMods <- c("S" = "p ~ 1; k ~ 1", "L" = "p ~ 1; k ~ 1",
#'              "M" = "p ~ 1; k ~ 1", "XL" = "p ~ 1; k ~ 1"
#'             )
#'   cpMods <- c("S" = "dist: exponential; l ~ 1; NULL", 
#'               "L" = "dist: exponential; l ~ 1; NULL",
#'               "M" = "dist: exponential; l ~ 1; NULL",
#'               "XL" = "dist: exponential; l ~ 1; NULL"
#'             )
#'   avgSS <- averageSS(mock$SS)
#'   gsGeneric <- estgGenericSize(nsim = 1000, days = avgSS,
#'                  modelSetSize_SE = pkmModsSize,
#'                  modelSetSize_CP = cpmModsSize,
#'                  modelSizeSelections_SE = pkMods,
#'                  modelSizeSelections_CP = cpMods
#'                )
#'  summary(gsGeneric)
#'
#' @export
#'
summary.gGenericSize <- function(object, ..., CL = 0.90){

  sizeclasses <- names(object)
  out <- list()
  for (sci in sizeclasses){
    out[[sci]] <- summary(object[[sci]], ..., CL = CL)
  }
  return(out)
}

#' @title Create search schedule data into an prepSS object for convenient 
#'  splits analyses
#'
#' @description Since data_SS columns largely have a specific, required
#'   format, the \code{prepSS} function can often automatically decipher the
#'   data, but the user may specify explicit instructions for parsing the data
#'   for safety if desired. If the data are formatted properly, the automatic
#'   parsing is reliable in most cases. There are two exceptions. (1) If
#'   there is more than one column with possible dates (formatted as formal
#'   dates (as class \code{Date}, \code{POSIXlt} or \code{POSIXct}) or
#'   character strings or factors that can be unambiguously interpreted as
#'   dates (with assumed format "2018-05-15" or "2018/5/15"). In that case,
#'   the user must specify the desired dates as \code{dateColumn}. (2) If
#'   there is a covariate column consisting entirely of 0s and 1s. In that
#'   case, the user must specify the column(s) in \code{covars}.
#'
#' @param data_SS data frame or matrix with search schedule parameters,
#'  including columns for search dates, covariates (describing characteristics
#'  of the search intervals), and each unit (with 1s and 0s to indicate 
#'  whether the given unit was searched (= 1) or not (= 0) on the given date)
#'
#' @param SSdate name of the column with the search dates in it
#'  (optional). If no \code{SSdate} is given, \code{prepSS} will
#'  try to find the date column based on data formats. If there is exactly one
#'  column that can be interpreted as dates, that column will be taken as the
#'  dates searched. If more than one date column is found, \code{prepSS} exits
#'  with an error message.
#'
#' @param preds vector of character strings giving the names of columns to be
#'  interpreted as potential covariates (optional). Typically, it is not
#'  necessary for a user to provide a value for \code{preds}. It is used only
#'  to identify specific columns of 1s and 0s as covariates rather than as
#'  search schedules.
#'
#' @return \code{prepSS} object that can be conveniently used in the splitting
#'  functions.
#'
#' @examples
#'  data(mock)
#'  prepSS(mock$SS)
#'
#' @export
#'
prepSS <- function(data_SS, SSdate = NULL, preds = NULL){
  if ("prepSS" %in% class(data_SS)) return(data_SS)
  if (length(intersect(class(data_SS), c("data.frame", "matrix"))) == 0){
    stop("data_SS must be a data frame or matrix")
  } else if (is.null(colnames(data_SS))){
    stop("data_SS columns must be named")
  }
  i <- sapply(data_SS, is.factor)
  data_SS[i] <- lapply(data_SS[i], as.character)

  # if SSdate not provided, extract search dates (if possible)
  if (is.null(SSdate)){
    for (coli in colnames(data_SS)){
      tmp <- checkDate(data_SS[, coli])
      if (!is.null(tmp)){
        if (!is.null(SSdate)){
          stop(
            "more than 1 date column in data_SS, and ",
            "SSdate does not specify which to use."
          )
        } else {
          SSdate <- coli
          dates <- tmp
        }
      }
    }
    if (is.null(SSdate)){
      stop("no columns can be interpreted as dates")
    }
  } else {
    if (length(SSdate) > 1 || !is.character(SSdate)){
      stop("'SSdate' must be NULL or the name of a single column")
    }
    dates <- checkDate(data_SS[, SSdate])
    if (is.null(dates)){
      stop(paste(SSdate, "is not properly formatted as dates"))
    }
  }
  # extract (potential) units  (i.e. cols w/ 0-1 data)
  unitNames <- NULL
  preds <- SSdate
  for (coli in colnames(data_SS)){
    if (coli %in% preds) next
    if (!is.numeric(data_SS[, coli])){
      preds <- c(preds, coli)
      next
    }
    if (sum(!(data_SS[, coli] %in% 0:1)) > 0){
      preds <- c(preds, coli)
      next
    } else {
      unitNames <- c(unitNames, coli)
    }
  }
  if (grepl("-",paste(unitNames, collapse = ''))){
    stop("Unit names must not contain hyphens ( - )")
  }
  date0 <- min(dates)
  ans <- list()
  ans$date0 <- date0
  ans$days <- as.numeric(difftime(dates, date0, units = "days"))
  if (any(diff(ans$days) <= 0)){
    stop("search dates must be in increasing order")
  }
  ans[[SSdate]] <- dates
  for (i in 1:length(preds)){
    if (preds[i] == SSdate) next
    if (is.factor(data_SS[, preds[i]])){
      ans[[preds[i]]] <- as.character(data_SS[,preds[i]])
    } else {
      ans[[preds[i]]] <- data_SS[,preds[i]]
    }
  }
  ans$searches_unit <- t(as.matrix(data_SS[, unitNames]))
  ans$unit <- unitNames
  rownames(ans$searches_unit) <- ans$unit
  ans$SSdate <- SSdate
  class(ans) <- c("prepSS", "list")
  return(ans)
}