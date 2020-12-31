#' @title Estimate the number of fatalities in each search interval throughout
#'   the monitoring period.
#'
#' @description A carcass that is observed in a given search may have arrived 
#'   at any time prior to that search, so carcass discovery time is often not
#'   a reliable estimate of carcass arrival time. For each observed carcass, 
#'   \code{calcRate} takes into account the estimated probability of arrival
#'   in each possible arrival interval, adjusts by detection probability, and 
#'   sums to estimate the estimated number of carcass arrivals in every search
#'   interval.
#'
#' @param M Numeric array (ncarc x nsim) of estimated number of fatalities
#'   by observed carcass and simulation rep
#'
#' @param Aj Integer array (ncarc x nsim) of simulated arrival intervals for
#'   each observed carcass. Arrival intervals are given as integers j, 
#'   indicating that the given carcass (indexed by row) arrived in the jth
#'   search interval in the given simulation rep (indexed by column). Arrival 
#'   interval indices (j) are relative to indexed carcasses' search schedules.
#'
#' @param data_SS \code{\link{prepSS}} object that contains formatted data for
#'   calculating splits. Optional argument. Alternatively, user may provide 
#'   \code{days} and \code{searches_carcass}.
#'
#' @param days Vector of all dates that at least one unit was searched. Format
#'   is the number of days since the first search. For example, days = c(0, 7,
#'   14, 28, 35) for a simple 7-day search schedule in which searches were
#'   conducted every once per week on the same day for 5 weeks. Not all units
#'   need be searched on every search date.
#'
#' @param searches_carcass An ncarc x length(days) array of 0s and 1s to 
#'   indicate searches in which the indexed carcass could have been found. 
#'   For example, row i = \code{c(1, 0, 1, 0, 1)} indicates that the search 
#'   schedule for the location (unit) where carcass i was found would be 
#'   \code{days[c(1, 3, 5)]}.
#'
#' @return Numeric array (nsim x nsearch) of estimated fatalities in each
#'   search interval. NOTE: The search at time t = 0 does not correspond to an
#'   interval, and all carcasses found at that time are assumed to have 
#'   arrived prior to the monitoring period and are not included in mortality 
#'   estimates so \code{nsearch = length(days) - 1}.
#'
#' @export
#'
calcRate <- function(M, Aj, days = NULL, searches_carcass = NULL,
                     data_SS = NULL){
  if (!is.null(data_SS) && ("prepSS" %in% class(data_SS))){
    days <- data_SS$days
    unit <- rownames(Aj)
    x <- nrow(Aj)
    searches_carcass <- array(0, dim = c(x, length(days)))
    for (xi in 1:x){
       searches_carcass[xi, ] <- data_SS$searches_unit[unit[xi], ]
    }
    searches_carcass[, 1] <- 1
  } else if (!is.null(data_SS)){
    stop(
      "In arg list, data_SS must be an prepSS object. Alternatively, user ",
      "may provide 'days' and 'carcass_searches'."
    )
  }
  if (is.vector(M)){
    M <- matrix(M, nrow = 1)
    Aj <- matrix(Aj, nrow = 1)
    searches_carcass <- matrix(searches_carcass, nrow = 1)
  }
  calcRateC(M, Aj, days, searches_carcass)
}

#' @title Estimate the number of fatalities by time interval
#'
#' @description \code{calcTsplit()} is a lower-level function that requires 
#'   the output of \code{calcRate} as input. See \code{\link{calcSplits}} 
#'   for a more powerful, convenient, and flexible alternative.
#'
#' @param rate Array (nsim x nsearch) of arrival rates as number of fatalities
#'   per search interval. Typically, \code{rate} will be the return value of 
#'   the \code{calcRate} function. 
#'
#' @param days A vector of times representing search dates when at least one
#'   unit was searched. Times are formatted as number of days since the first
#'   search, e.g., c(0, 7, 14, 28, 35) would indicate a schedule in at least 
#'   one unit was searched every 7 days.
#'
#' @param tsplit A vector of times that splits the monitoring period into a
#'   set of time intervals for which \code{calcTsplit} will estimate the 
#'   number of fatalities. For example, if \code{tsplit = c(0, 14, 19, 35)}, 
#'   then \code{calcTsplit} estimates the number of fatalities occuring in 
#'   interval (0, 14], (14, 19], and (19, 35]. Times in \code{tsplit} must be 
#'   increasing and between 0 and max(days), inclusive.
#'
#' @return A numeric array with dimensions
#'  \code{dim = c(length(tsplit) - 1, nsim)} giving the estimated number of
#'  fatalities that occured in each time interval.
#'
#' @export
#'
calcTsplit <- function(rate, days, tsplit){
  if (is.vector(rate)){
    rate <- matrix(rate, nrow = 1)
  }
  calcTsplitC(rate, days, tsplit)
}

#' @title Estimate the number of fatalities by up to two splitting covariates
#'
#' @description Total mortality can be split into sub-categories, according to
#'   various splitting covariates such as species, visibility class, season, 
#'   site, unit, etc. Given the carcass search data, estimated mortalities, 
#'   and splitting covariates, \code{calcSplits()} gives the "splits" or 
#'   summaries the estimated mortalities by levels of the splitting 
#'   covariates. For example, user may specify \code{"season"} and 
#'   \code{"species"} as splitting variables to see estimated mortalities by
#'   season and species. Input would be arrays of estimated mortalities and
#'   arrival intervals when \code{ncarc} carcass have been discovered and 
#'   uncertainty in mortality estimates is captured via simulation with 
#'   \code{nsim} simulation draws.
#'
#' @details Arrival intervals (\code{Aj}) are given as integers, j, that
#'  indicate which search interval the given carcass (indexed by row) arrived
#'  in the given simulation draw (indexed by column). Arrival interval indices
#'  (j) are relative to indexed carcasses' search schedules.
#'
#'  No more than two splitting variables (\code{split_CO}, \code{split_SS}, 
#'  and \code{split_time}) in total may be used. \code{split_CO} variables
#'  describe qualitative characteristics of the observed carcasses or where
#'  they were found. Some examples include searcher (DHD, JPS, MMH), carcass
#'  size (S, M, L), species, age (fresh/dry or immature/mature), unit,
#'  visibility class (easy, moderate, difficult), etc.
#'
#' \code{split_SS} variables describe characteristics of the search intervals,
#'   such as season (spring, summer, fall, winter) or treatment
#'   (pre- or post-minimization). Each search interval is assigned a level of 
#'   the \code{split_SS} variable. For example, for a search schedule with
#'   5 searches (including a search at t = 0), and the \code{split_SS} 
#'   variable would have values for each of the 4 search intervals. The
#'   levels of the \code{split_SS} must be in contiguous blocks. For example,
#'   \code{season = c("S", "S", "F", "F")} would be acceptable, but
#'   \code{season = c("S", "F", "S", "F")} would not be.
#'
#' \code{split_time} variables are numeric vectors that split the monitoring
#'   period into distinct time intervals. For example,
#'   \code{split_time = c(0, 30, 60, 90, 120)} would split the 120 monitoring
#'   period into 30-day intervals, and \code{calcSplits()} would return 
#'   mortality estimates for each of the intervals.
#'
#' @param M \code{\link{estM}} object, containing numeric array (ncarc x nsim)
#'    of estimated mortalities and other pieces
#'   
#' @param split_CO Character vector of names of splitting covariates to be
#'   found in the \code{data_CO} data frame. No more than two \code{split_CO} 
#'   variables are allowed. Use \code{split_CO = NULL} if no CO splits are 
#'   desired. 
#'   
#' @param data_CO data frame that summarizes the carcass search data and must
#'   include columns specified by the \code{split_CO} arg. Each row includes
#'   search and discovery parameters associated with a single observed 
#'   carcass. Columns include carcass ID, carcass discovery date, unit, and 
#'   any number of covariates. \code{data_CO} is required if and only if 
#'   \code{split_CO} is non-NULL.
#'   
#' @param split_SS Character string giving the name of a splitting covariate 
#'   in the \code{data_SS} list, with \code{data_SS[[split_SS]]} describing
#'   characteristics of the search intervals (e.g., "season"). Note that
#'   \code{length(data_SS[[split_SS]]} must equal 
#'   \code{length(data_SS$days) - 1} because no inference is made about
#'   carcass arrivals prior to time t = 0, and the "interval" prior to t = 0 
#'   is not taken as a "search interval." If no \code{split_SS} split is 
#'   desired, use \code{split_SS = NULL}.
#'   
#' @param data_SS Search schedule data
#'   
#' @param split_time Numeric vector that defines time intervals for splits.
#'  Times must be numeric, strictly increasing, and span the monitoring period
#'  [0, \code{max(data_SS$days)}]. If no \code{split_time} is desired, use
#'  \code{split_time = NULL}. If \code{split_time} is NULL and \code{split_SS}
#'  is not NULL, \code{data_SS} is required.
#'   
#' @param ... arguments to be passed down
#'   
#' @return An object of class \code{splitFull} is returned. If one splitting
#'  covariate is given, then the output will be an array of estimated
#'  mortality in each level of the splitting covariate, with one row for each
#'  covariate level and one column for each simulation draw. If two splitting
#'  covariates are given, output will be a list of arrays. Each array gives
#'  the estimated mortalities for one level of the second splitting covariate
#'  and all levels of the first splitting covariate.
#'
#'  Objects of class \code{splitFull} have attributes \code{vars} (which gives
#'  the name of the splitting covariate(s)) and \code{type} (which specifies
#'  whether the covariate(s) are of type \code{split_CO}, \code{split_SS}, or
#'  \code{split_time}). A summary of a resulting \code{splitFull} object
#'  is returned from the S3 function \code{summary(splits, CL = 0.90, ...)},
#'  which gives the mean and a 5-number summary for each level of each
#'  covariate. The 5-number summary includes the alpha/2, 0.25, 0.5, 0.75,
#'  and 1 - alpha/2 quantiles, where alpha = 1 - CL. A graph summarizing the
#'  results can be drawn using \code{plot(splits, CL, ...)}, which gives
#'  a graphical representation of the \code{summary}.
#'
#' @examples
#'  \donttest{
#'   model_SE <- pkm(p ~ 1, k ~ 1, data = wind_RPbat$SE)
#'   model_CP <- cpm(l ~ 1, s ~ 1, data = wind_RPbat$CP, dist = "weibull",
#'     left = "LastPresent", right = "FirstAbsent")
#'   Mhat <- estM(nsim = 1000, data_CO = wind_RPbat$CO, 
#'     data_SS = wind_RPbat$SS, data_DWP = wind_RPbat$DWP, 
#'     model_SE = model_SE, model_CP = model_CP,
#'     unitCol = "Turbine", COdate = "DateFound")
#'
#'   M_spp <- calcSplits(M = Mhat, split_CO = "Species",
#'     data_CO = wind_RPbat$CO)
#'   summary(M_spp)
#'   plot(M_spp)
#'  }
#' @export
#'
calcSplits <- function(M, split_CO = NULL, data_CO = NULL,
                       split_SS = NULL, data_SS = NULL, split_time = NULL,
                       ...){
  ind <- sapply(data_CO, is.factor)
  data_CO[ind] <- lapply(data_CO[ind], as.character)
  ind <- sapply(data_SS, is.factor)
  data_SS[ind] <- lapply(data_SS[ind], as.character)

  ##### read data and check for errors
  if ("list" %in% class(M)){
    if (!all(c("Mhat", "Aj") %in% names(M))){
      stop("M must be a list with elements $M and $Aj")
    }
    Aj <- M$Aj
    M <- M$Mhat
    if (length(dim(M)) != 2){
      stop("M must be a list with elements $M and $Aj")
    }
    if (!is.numeric(M) || anyNA(M)){
      stop("M must be numeric and not contain missing values.")
    }
  } else {
    stop("M must be a list with elements $M and $Aj")
  }

  cleanInd <- which(Aj[, 1] == 0)
  if (length(intersect(split_CO, split_SS)) > 0){
    stop("CO and SS splitting variables must have distinct names")
  }
  if ((!is.null(split_SS) || !is.null(split_time)) && is.null(data_SS)){
    stop("data_SS must be provided if ",
      ifelse(is.null(split_SS), "split_time ", "split_SS "), "is")
  }
  if (!is.null(split_CO) || !is.null(split_SS) || !is.null(split_time)){
    if (is.null(data_CO)){
      stop("data_CO must be provided to perform non-null splits")
    } else if (length(cleanInd) > 0){
      if (dim(data_CO)[1] == dim(Aj)[1]){
        data_CO <- data_CO[-cleanInd,]
        Aj <- Aj[-cleanInd, ]
        M <- M[-cleanInd, ]
      } else {
        # better to check carcass IDs? or input data as a class?
        stop("mismatch between carcass numbers in data_CO and Aj")
      }
    }
  }
  unit <- rownames(Aj)
  # better than attaching units to rownames of Aj in estg would be to package
  # the unit and carcass ID to the list returned by estg and include these in
  # class estM

  if (!is.null(split_SS) || !is.null(split_time)){
    if (!("prepSS" %in% class(data_SS))) data_SS <- prepSS(data_SS)
  }
  ### declare traffic directing variables (updated after parsing inputs)
  # number of valid split variables:
  nvar <- 0
  # characterization of horizontal and vertical split variables as "SS"
  # or "CO":
  vartype <- NULL
  split_h <- NULL
  split_v <- NULL
  ### error-checking the split variables and interpreting inputs:
  if (!is.null(split_SS) && is.Date(data_SS[[split_SS]])){
    # split_SS is a temporal split
    split_time <- as.numeric(data_SS[[split_SS]]-data_SS$date0)
    split_SS <- NULL
  }
  if (!is.null(split_time)){
    if (!is.numeric(split_time) || !is.vector(split_time)){
      stop("split_time must be NULL or a numeric vector")
    }
    minspl <- min(split_time)
    if (minspl != 0){
      if (minspl < 0){
        stop("min(split_time) = ", minspl, " but must not be < 0" )
      }
      if (minspl > 0){
        warning(
          "min(split_time) = ", minspl, ". ",
          "Appending 0 to split_time, so first split is [0, ", minspl, "]."
        )
      }
    }
    maxspl <- max(split_time)
    maxdays <- max(data_SS$days)
    if (maxspl != maxdays){
      if (maxspl > maxdays){
        stop(
          "max(split_time) = ", maxspl, " ",
          "but must not be > max(data_SS$days) = ", maxdays
        )
      }
      if (maxspl > maxdays){
        warning(
          "max(split_time) = ", maxspl, " < max(data_SS$days) = ", maxdays, 
          ". ", "Appending max(days) to split_time so last split is ",
          "[", maxspl, ", ", max(data_SS$days),"]."
        )
      }
    }
    if (sum(diff(split_time) <= 0) > 0){
      stop("split_time must be strictly increasing")
    }
    # split_time has passed the initial error-checking, so assign attributes
    # to characterize the data
    split_h <- list()
    split_h[["name"]] <- "time"
    split_h[["vals"]] <- split_time
    if (minspl > min(data_SS$days)) split_h$vals <- c(0, split_h$vals)
    if (maxspl < maxdays) split_h$vals <- c(split_h$vals, maxdays)
    split_h$vals[length(split_h$vals)] <- data_SS$days[length(data_SS$days)]
    split_h[["level"]] <- unique(split_h$vals[-1])
    split_h[["nlev"]] <- length(split_h$level) - 1
    split_h[["type"]] <- "time"
  }
  if (!is.null(split_SS)){ # there is a split_SS variable
    if (!is.null(split_h)){
      stop(
        "Only one temporal split allowed. ",
        "Either split_SS or split_time must be NULL"
      )
    }
    if (!is.character(split_SS)){
      stop("split_SS must be NULL or the name of an element in data_SS")
    }
    if (length(split_SS) > 1){
      stop("At most 1 split_SS variable is allowed")
    }
    if (!(split_SS %in% names(data_SS))){
      stop(split_SS, " not found in ", data_SS)
    }
    if (is.numeric(data_SS[[split_SS]])){
      stop("split_SS column must be categorical or dates")
    }
    SSlevel <- data_SS[[split_SS]]
    if (length(unique(SSlevel)) == 1)
      stop("split_SS = ", split_SS, "has only one level...nothing to split on")
    if (min(diff(match(SSlevel, unique(SSlevel)))) < 0){
      stop(
        "split_SS levels must be in contiguous blocks in data_SS.\n",
        "For example, c('spring', 'spring', 'fall', 'fall') would be OK, ",
        "but c('spring', 'fall', 'fall', 'spring') is not."
      )
    }
    split_h <- list()
    split_h[["name"]] <- split_SS
    tmp <-
      cumsum(table(data_SS[[split_h$name]])[unique(data_SS[[split_h$name]])])
    split_h[["vals"]] <- c(0, data_SS$days[tmp])
    split_h[["level"]] <- unique(data_SS[[split_h$name]])
    split_h[["nlev"]] <- length(split_h$level)
    split_h[["type"]] <- "SS"
  }
  if (!is.null(split_CO)){
    if (length(split_CO) > 2 || (length(split_CO) > 1 & !is.null(split_h))){
      stop(
        "At most two split variables are allowed in total, i.e. ",
        "length(split_CO) + length(split_SS) + length(split_time) ",
        "must be <= 2"
      )
    }
    if (!is.character(split_CO)){
      stop("split_CO must be a name of a selected column in data_CO")
    }
    if (!all(split_CO %in% names(data_CO))){
      stop(split_CO[which(!(split_CO %in% names(data_CO)))], "not in data_SS")
    }
    if (length(split_CO) == 2 || !is.null(split_h)){
      split_v <- list()
      split_v[["name"]] <- ifelse(!is.null(split_h), split_CO[1], split_CO[2])
      split_v[["vals"]] <- data_CO[[split_v$name]]
      split_v[["level"]] <- gtools::mixedsort(unique(split_v$vals))
      split_v[["nlev"]] <- length(split_v$level)
      split_v[["type"]] <- "CO"
    }
    if (is.null(split_h)){
      split_h <- list()
      split_h[["name"]] <- split_CO[1]
      split_h[["vals"]] <- data_CO[[split_h$name]]
      split_h[["level"]] <- gtools::mixedsort(unique(split_h$vals))
      split_h[["nlev"]] <- length(split_h$level)
      split_h[["type"]] <- "CO"
    }
  }
  nvar <- sum(c(!is.null(split_h), !is.null(split_v)))

  # additional preprocessing
  if (is.vector(M)) M <- matrix(M, nrow = 1)
  x <- dim(M)[1] # total observed carcasses (assumes data_CO is error-checked)
  nsim <- dim(M)[2] # number of simulation draws (columns in M)
  if (!is.null(split_h) && (split_h$type %in% c("time", "SS"))){
      days <- data_SS$days
      searches_carcass <- array(0, dim = c(x, length(days)))
      for (xi in 1:x){
        searches_carcass[xi, ] <- data_SS$searches_unit[unit[xi], ]
      }
      searches_carcass[, 1] <- 1
  }
  # calculate splits
  splits <- list()
  if (nvar == 0){ # no splits...just calculate total M
    splits[["M"]] <- colSums(M)
    splits[["X"]] <- nrow(M)
  } else if (nvar == 1){ # just one split variable: split_h
    if (split_h$type == "CO"){
      splits[["M"]] <- array(dim = c(split_h$nlev, nsim))
      splits[["X"]] <- numeric(split_h$nlev)
      for (li in 1:split_h$nlev) {
        lind <- which(data_CO[, split_h$name] == split_h$level[li])
        if (length(lind) == 1){
          splits[["M"]][li, ] <- M[lind, ]
        } else {
          splits[["M"]][li, ] <- colSums(M[lind, ])
        }
        splits[["X"]][li] <- length(lind)
      }

    } else if (split_h$type %in% c("time", "SS")){
      days <- data_SS$days
      rate <- calcRate(M, Aj, days = days, 
                searches_carcass = searches_carcass)
      splits[["M"]] <- calcTsplit(rate, data_SS$days, split_h$vals)
      ratex <- calcRate(M^0, Aj = Aj, days = days, 
                 searches_carcass = searches_carcass)
      splits[["X"]] <- calcTsplit(ratex, data_SS$days, split_h$vals)
      splits[["X"]][which(splits[["M"]] == 0)] <- 0
      splits[["X"]] <- rowMeans(splits[["X"]])
    }
  } else if (nvar == 2){ # two split variables: split_h and split_v
    splits <- list()
    if (split_h$type == "CO"){
      for (vi in 1:split_v$nlev){
        splits[["M"]][[vi]] <- array(0, dim = c(split_h$nlev, nsim))
        splits[["X"]][[vi]] <- numeric(split_h$nlev)
        for (li in 1:split_h$nlev) {
          lind <- which(
            split_h$vals == split_h$level[li] &
            split_v$vals == split_v$level[vi]
          )
          if (length(lind) > 1){
            splits[["M"]][[vi]][li, ] <- colSums(M[lind, ])
          } else if (length(lind) == 1){
            splits[["M"]][[vi]][li, ] <- M[lind, ]
          }
          splits[["X"]][[vi]][li] <- length(lind)
        }
      }
    } else if (split_h$type %in% c("SS", "time")){
      for (vi in 1:split_v$nlev){
        lind <- which(split_v$vals == split_v$level[vi])
        rate <- calcRate(M = M[lind, ], Aj = Aj[lind, ], days = data_SS$days,
          searches_carcass = searches_carcass[lind,])
        splits[["M"]][[vi]] <- calcTsplit(rate, data_SS$days, split_h$vals)
        ratex <- calcRate(M[lind, ]^0, Aj = Aj[lind, ], days = data_SS$days,
                   searches_carcass = searches_carcass[lind,])
        splits[["X"]][[vi]] <- calcTsplit(ratex, data_SS$days, split_h$vals)
        splits[["X"]][[vi]][which(splits[["M"]][[vi]] == 0)] <- 0
        splits[["X"]][[vi]] <- rowMeans(splits[["X"]][[vi]])
      }
    }
  }
  #protection against unintended loss of attr's
#  splits <- sticky::sticky(splits)
# sticky() helps preserve object attributes when subsetted, but the function is
# being removed from CRAN. It is unlikely that the loss of it here will cause
# problems.

  attr(splits, "vars") <- c(split_h$name, split_v$name)
  attr(splits, "type") <- c(split_h$type, split_v$type)
  if (!is.null(split_h) && (split_h$type %in% c("time", "SS"))){
    attr(splits, "times") <- split_h$vals
  }
  if (nvar == 1){
    rownames(splits[["M"]]) <- split_h$level
    names(splits[["X"]]) <- split_h$level
  }
  if (nvar == 2){
    names(splits[["M"]]) <- split_v$level
    names(splits[["X"]]) <- split_v$level
    for (i in 1:length(splits[["M"]])){
      rownames(splits[["M"]][[i]]) <- split_h$level
      names(splits[["X"]][[i]]) <- split_h$level
    }
  }
  class(splits) <- "splitFull"
  return(splits)
}

#' @title Summarize results of mortality estimate splits
#'
#' @description Mortality estimates can be calculated for the various levels 
#'   of splitting covariates such as season, species, or visibility class
#'   using \code{\link{calcSplits}}, which gives full arrays of simulated M 
#'   estimates (i.e., for each level of each splitting covariate, each 
#'   discovered carcass, and each simulation draw). summary(splits, CL = 0.90,
#'   ...) gives summary statistics of the estimates.
#'
#' @param object A \code{splitFull} object (\code{\link{calcSplits}}) that 
#'   gives simulated mortality estimates for all combinations of levels of 1 
#'   or 2 splitting covariates.
#'   
#' @param CL desired confidence level for summary CIs (numeric scalar in 
#'  (0, 1))
#'   
#' @param ... to be passed down
#'   
#' @return an object of class \code{splitSummary}, which gives 5-number
#'   summaries for all combinations of levels among the splitting covariates 
#'   in the \code{splits}. The 5-number summaries include the mean and 
#'   alpha/2, 0.25, 0.5, 0.75, and 1 - alpha/2 quantiles of mortality 
#'   estimates, where alpha = 1 - CL. A graphical representation of the 
#'   results can be produced using \code{plot(splits, CL, ...)}. For splits
#'   along CO covariates, the levels are organized alphabetically (but with 
#'   numeric suffixes appearing in numeric order, e.g., "t1", "t2", "t10" 
#'   rather than "t1", "t10", "t2").
#'   
#' @export
#'
summary.splitFull <- function(object, CL = 0.90, ...){
  splits <- object
  alpha <- 1 - CL
  probs <- c(alpha/2, 0.25, 0.5, 0.75, 1 - alpha/2)
  if (length(attr(splits, "vars")) == 0){
    sumry <- c(quantile(splits$M, probs = probs))
    sumry <- pmax(sumry, splits$X)
    sumry <- c(X = splits$X, sumry)
    attr(sumry, "vars") <- NA
    attr(sumry, "type") <- NA
    attr(sumry, "times") <- 1
  } else if (length(attr(splits, "vars")) == 1){
    if (is.vector(splits$M)){
      splits$M <- matrix(splits$M, nrow = 1)
      splits$X <- matrix(splits$X, nrow = 1)
    }
    sumry <- matrixStats::rowQuantiles(splits$M, probs = probs)
    ind <- (sumry < splits$X)
    sumry <- (sumry * !ind) + (splits$X * ind)
    sumry <- cbind(X = splits$X, sumry)
  } else if (length(attr(splits, "vars")) == 2){
    if (is.vector(splits$M[[1]])){
      splits$M <- lapply(splits$M, function(x) matrix(x, nrow = 1))
      splits$X <- lapply(splits$X, function(x) matrix(x, nrow = 1))
    }
    sumry <- lapply(splits$M, function(x){
      cbind(matrixStats::rowQuantiles(x, probs = probs))
    })
    for (levi in 1:length(sumry)){
      ind <- (sumry[[levi]] < splits$X[[levi]])
      sumry[[levi]] <- (sumry[[levi]] * !ind) + (splits$X[[levi]] * ind)
      sumry[[levi]] <- cbind(X = splits$X[[levi]], sumry[[levi]])
    }
  } else {
    stop(
      "length(attr(splits, 'vars')) > 2.",
      "At most two split variables are allowed."
    )
  }
  # order the non-temporal dimensions "alphabetically" (done in calcSplits?)
  if (!is.null(attr(splits, "type"))){
    if (!is.list(splits$M)){
      if (attr(splits, "type") == "CO"){
        sumry <- sumry[gtools::mixedsort(rownames(sumry)), ]
      }
    } else {
      if (attr(splits, "type")[1] == "CO"){
        for (i in 1:length(splits)){
          sumry[[i]] <- sumry[[i]][gtools::mixedsort(rownames(sumry[[i]])), ]
        }
      }
      if (attr(splits, "type")[2] == "CO"){
        sumry <- sumry[gtools::mixedsort(names(sumry))]
       }
     }
     attr(sumry, "CL") <- CL
     attr(sumry, "vars") <- attr(splits, "vars")
     attr(sumry, "type") <- attr(splits, "type")
     attr(sumry, "times") <- attr(splits, "times")
   }
#  sumry <- sticky::sticky(sumry)
# sticky() helps preserve object attributes when subsetted, but the function is
# being removed from CRAN. It is unlikely that the loss of it here will cause
# problems because subsetting splitSummary objects by hand would be so unusual.

  class(sumry) <- "splitSummary"
  return(sumry)
}

#' @title Transpose a list of arrays
#'
#' @description A list of \code{n} arrays, each with dimension \code{m} x 
#'   \code{k} is redimensioned to a list of \code{m} arrays, each with 
#'   dimension \code{m} x \code{k}. NOTE: Attributes are not preserved.
#'
#' @param M a list of \code{n} \code{m} x \code{k} arrays
#'   
#' @return a list of \code{m} \code{n} x \code{k} arrays
#'
#' @export
#'
ltranspose <- function(M){
  if (!is.list(M))
    stop(substitute(M), " must be a list")
  if (!is.matrix(M[[1]]))
    stop("elements of ", substitute(M), " must be arrays.")
  adim <- dim(M[[1]])
  for (i in 1:length(M)){
    if (!isTRUE(all.equal(dim(M[[i]]), adim))){
      stop(
        "elements of ", substitute(M),
        " must be arrays with the same dimensions"
      )
    }
  }
  ans <- list()
  for (i in 1:adim[1]){
    ans[[i]] <- do.call("rbind", lapply(M, function(x) x[i, ]))
  }
  return(ans)
}

#' @title Transpose a \code{splitFull} array (preserving attributes)
#'
#' @description Transpose a \code{splitFull} array (preserving attributes)
#'   
#' @param splits a \code{splitFull} object, which is a list of \code{n}
#'   \code{m} x \code{k} arrays with attributes describing characteristics of 
#'   the splits
#'   
#' @return a list of \code{m} \code{n} x \code{k} arrays as a \code{splitFull}
#'  object
#'   
#' @export
#'   
transposeSplits <- function(splits){
  if (!(class(splits) %in% "splitFull")){
    stop(
      substitute(splits), " is not a splitFull object. ",
      "Only splitFull objects can be transposed using transposeSplits()"
    )
  }
  ans <- list()
  ans$M <- ltranspose(splits$M)
  tmp <- ltranspose(lapply(splits$X, 
           FUN = function(x) array(x, dim = c(length(x), 1))))
  ans$X <- lapply(tmp, as.vector)
  names(ans) <- c("M", "X")
  names(ans$M) <- names(ans$X) <- rownames(splits$M[[1]])
  class(ans) <- "splitFull"
  attr(ans, "vars") <- attr(splits, "vars")[2:1]
  attr(ans, "type") <- attr(splits, "type")[2:1]
  attr(ans, "times") <- attr(splits, "times")
  return(ans)
}
