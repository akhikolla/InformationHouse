#' @title Compute the logit or anti-logit
#' 
#' @param x A number. For  \code{logit}, a probability (between 0 and 1, 
#'   inclusive). For  \code{alogit}, any real number.
#'
#' @return \code{logit}: The logit of \code{x}.
#'
#' @examples
#'   logit(0.5)
#'
#' @export 
#'
logit <- function(x) {
  log(x / (1 - x))
}

#' @rdname logit
#' 
#' @return \code{alogit}:  The anti-logit of \code{x}.
#'
#' @examples
#'   alogit(0)
#'
#' @export 
#'
alogit <- function(x) {
  1 / (1 + exp(-x))
}

#' @title The CDF of the loglogistic distribution
#'
#' @param q a numeric vector of quantiles
#'
#' @param pda the \eqn{\alpha} parameter
#'
#' @param pdb the \eqn{\beta} parameter
#'
#' @details There are several common parameterizations of the loglogistic
#'  distribution. The one used here gives the following:
#'  \describe{
#'    \item{CDF}{\code{Pr(X <= x) = 1/(1 + (x/}\eqn{\beta})^-\eqn{\alpha}\code{)}}
#'    \item{PDF}{\code{Pr(X = x) = (}\eqn{\alpha}/\eqn{\beta}\code{) * (x/}\eqn{\beta})^(\eqn{\alpha}\code{ - 1)/(1 + (x/}\eqn{\beta})^\eqn{\alpha}\code{)^2}}
#'   }
#' @return \code{Pr(X <= q | pda, pdb)}
#'
#' @export
#'
pllogis <- function(q, pda, pdb){
  if (!is.numeric(q) || any(q < 0)) stop("q must be non-negative")
  if (!is.numeric(pda) || any(pda <= 0)) stop("pda must be positive")
  if (!is.numeric(pdb) || any(pdb <= 0)) stop("pdb must be positive")
  return(1/(1 + (q/pdb)^-pda))
}
#' @title Calculate day of study from calendar date
#'
#' @description Convert calendar date to integer day from a reference date
#'   (\code{ref}).
#'
#' @param date A date or vector of dates to convert to days.
#'
#' @param ref Reference date.
#'
#' @return Numeric value(s) of days from \code{ref}.
#'
#' @examples 
#'   x <- c("2018-01-01", "2018-02-01")
#'   dateToDay(x, x[1])
#'
#' @export 
#'
dateToDay <- function(date, ref = NULL){
  if (is.numeric(date)){
    return(date)
  }
  if (is.null(ref)){
    ref <- date
  }
  day <- as.numeric(difftime(date, ref, units = "days"))
  return(day)
}

#' Checks whether a vector of data can be interpreted as dates
#'
#' @description Checks whether the dates are in a standard format and 
#'  sensible. If so, function returns the dates converted to ISO 8601 
#'  yyyy-mm-dd format. Acceptable formats are yyyy-mm-dd, yyyy/mm/dd, 
#'  mm/dd/yyyy, and dd/mm/yyyy. If format is mm/dd/yyyy or dd/mm/yyyy, the 
#'  dates must be interpretable unambiguously. Also, dates must be later than 
#'  1900-01-01. This additional check provides some protection against common 
#'  data entry errors like entering a year as 0217 or 1017 instead of 2017.
#'
#' @param testdate Date(s) to check and format.
#'
#' @return dates formatted as yyyy-mm-dd (if possible) or NULL (if some value
#'  is not interpretable as a date after 1900-01-01).
#'
#' @examples
#'  checkDate("02/20/2018")
#'  checkDate("10/08/2018")
#'
#' @export
#'
checkDate <- function(testdate){
  beginningOfTime <- as.Date("1900-01-01")
  canDate <- try(as.Date(testdate), silent = TRUE)
  if (!("try-error" %in% class(canDate)) &&
      !anyNA(canDate) &&
      all(canDate > beginningOfTime)) return (canDate)

  formats <- list("%m/%d/%Y", "%d/%m/%Y", "%Y/%m/%d")
  canDate <- lapply(formats, 
               function(x) try(as.Date(testdate, tryFormats = x), 
                               silent = TRUE))
  canForm <- which(lapply(canDate, class) != "try-error")
  if (length(canForm) == 0) return (NULL)
  for (i in 1:length(canForm)){
    if (anyNA(canDate[[canForm[i]]])){
      canForm[i] <- NA
    } else if (any(canDate[[canForm[i]]] < beginningOfTime)) {
      canForm[i] <- NA
    }
  }
  if (all(is.na(canForm))) return(NULL)
  canForm <- canForm[!is.na(canForm)]
  if (length(canForm) > 1){
    warning("Cannot resolve ambiguous date format. Use yyyy-mm-dd for dates")
    return(NULL)
  }
  return(canDate[[canForm]])
}

#' @title Generic S3 function for summarizing AICc
#'
#' @description Extract AICc values from \code{pkm}, \code{pkmSet},
#'  \code{pkmSetSize}, \code{cpm}, \code{cpmSet}, and \code{cpmSetSize}.
#'
#' @param x Model or list of models to extract AICc values from.
#'
#' @param ... further arguments passed to or from other methods
#'
#' @return list of models sorted by AICc
#'
#' @export
aicc <- function(x, ... ){
  UseMethod("aicc", x)
}


#' @title Generic S3 function for printing corpus_frame
#'
#' @param x data frame to print
#'
#' @param ... other arguments
#'
#' @return prints data frame
#'
#' @export
print.corpus_frame <- function(x, ...){
  corpus::print.corpus_frame(x, rows = 80)
}

#' @title Auto-parsing to find the name of the unit column (\code{unitCol})
#'
#' If a unit column is not explicitly defined by user in the arg list to
#'  \code{estM} or \code{estg}, then \code{defineUnitCol} parses the CO, DWP,
#'  and SS files to extract the unit column if possible.
#'
#' Criteria that a column must meet to be a unit column are that it is found
#'  in both \code{data_CO} and \code{data_DWP}, all units in \code{data_CO} must
#'  also be included among units in \code{data_DWP}, all units in both
#'  \code{data_CO} and \code{data_DWP} must be included among the column names
#'  in \code{data_SS}. If \code{data_DWP = NULL}, then the unit column must be
#'  included in \code{data_CO} and all its units must be included among the
#'  column names of \code{data_SS}.
#'
#'
#' @param data_CO carcass observation data (data frame)
#'
#' @param data_SS search schedule data (data frame)
#'
#' @param data_DWP density-weighted proportion data (data frame)
#'
#' @return name of unit column (\code{unitCol}), if a unique unit column can be
#'  identified. If no unit column is present or there is more than one unit
#'  column, \code{defineUnitCol} stops with an error.
#'
#' @export

defineUnitCol <- function(data_CO, data_SS = NULL, data_DWP = NULL){
  if (!"prepSS" %in% class(data_SS)){
    ind <- sapply(data_SS, is.factor)
    data_SS[ind] <- lapply(data_SS[ind], as.character)
  }
  ind <- sapply(data_CO, is.factor)
  data_CO[ind] <- lapply(data_CO[ind], as.character)
  ind <- sapply(data_DWP, is.factor)
  data_DWP[ind] <- lapply(data_DWP[ind], as.character)

  if (is.null(data_CO)) stop("data_CO empty")
  if (is.null(data_SS) && is.null(data_DWP))
    stop("cannot find an unambiguous unit column in data_CO")
  if (is.null(data_SS)){ # unique column with matching names and levels in CO and DWP?
    unitCol <- intersect(colnames(data_CO), colnames(data_DWP))
    if (length(unitCol) == 0){
      stop(
        "No columns in data_CO and data_DWP share a common name to use as a ",
        "unit column. Cannot estimate M"
      )
    }
    if (length(unitCol) == 1 && any(!(data_CO[ , unitCol] %in% data_DWP[ , unitCol]))){
      ind <- which(!(data_CO[ , unitCol] %in% data_DWP[ , unitCol]))
      if (length(ind) > 1)
        stop("Units ", paste0(data_CO[ind, unitCol], collapse = ", "), " are ",
          "represented in data_CO but not in data_DWP. Cannot estimate M.")
        stop("Unit ", data_CO[ind, unitCol], " is ",
          "represented in data_CO but not in data_DWP. Cannot estimate M.")
    }
    if (length(unitCol) > 1){
      badind <- NULL
      for (ui in unitCol){
        if (length(data_DWP[ , ui]) != length(unique(data_DWP[, ui]))){
          badind <- c(badind, ui)
          next
        }
        if (any(! data_CO[, ui] %in% data_DWP[ ,ui])){
          badind <- c(badind, ui)
          next
        }
      }
      if (length(badind) > 0){
        unitCol <- unitCol[!unitCol %in% badind]
      }
      if (length(unitCol) == 0){
        stop("No shared column in DWP and CO files meets the criteria for ",
             "a unit column: either no common column name, or some CO ",
             "'units' not represented in corresponding DWP 'units', or ",
             "'units' listed in DWP are not unique.")
      }
    }
    return(unitCol)
  }
  if (is.null(data_DWP)){
    unitCol <- colnames(data_CO)
  } else {
    unitCol <- intersect(colnames(data_CO), colnames(data_DWP))
    if (length(unitCol) == 0){
      stop(
        "No columns in data_CO and data_DWP share a common name to use as a ",
        "unit column. Cannot estimate M"
      )
    }
    if (length(unitCol) == 1 & any(!(data_CO[ , unitCol] %in% data_DWP[ , unitCol]))){
      ind <- which(!(data_CO[ , unitCol] %in% data_DWP[ , unitCol]))
      if (length(ind) > 1)
        stop("Units ", paste0(data_CO[ind, unitCol], collapse = ", "), " are ",
          "represented in data_CO but not in data_DWP. Cannot estimate M.")
      stop("Unit ", data_CO[ind, unitCol], " is ",
        "represented in data_CO but not in data_DWP. Cannot estimate M.")
    }
  }
  if ("prepSS" %in% class(data_SS)){
    SSuname <- data_SS$unit
  } else {
    SSuname <- names(data_SS)
  }
  if (length(unitCol) == 1){
    ind <- !(data_CO[ , unitCol] %in% names(data_SS))
    if (any(ind))
      stop("Carcasses were found (CO) at units (",
        paste0(data_CO[ind, unitCol], collapse = "), "), " which are not ",
        "included in the search schedule (SS). Cannot estimate g or M.")
  }
  if (length(unitCol) > 1){
    bad <- NULL
    for (ni in 1:length(unitCol)){
      if (any(!(data_CO[ , unitCol[ni]] %in% SSuname))){
        bad <- c(bad, ni)
        next
      }
    }
    if (length(bad) != length(unitCol) - 1){
      stop(
        "No unitCol provided, and data_CO and data_DWP do not have a column ",
        "that can unambiguously serve as unitCol. Cannot estimate M."
      )
    } else {
      unitCol <- unitCol[-bad]
    }
  }
  unitCol
}
