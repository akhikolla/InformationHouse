#' @title Fit density-weighted proportion (DWP) models.
#'
#' @description Carcass density is modeled as a function of distance from
#'  turbine. Format and usage parallel that of common \code{R} functions
#'  \code{lm}, \code{glm}, and \code{gam} and the GenEst functions \code{pkm}
#'  and \code{cpm}.
#'
#' @details The fraction of carcasses falling in the area searched at a turbine
#'  may be a function of carcass class (e.g., large or small) and/or direction
#'  from the turbine. Data may be provided for fitting a distance model(s) or,
#'  alternatively, simulated turbine-wise DWP data from custom-fitted models may
#'  be provided. If pre-fit, pre-simulated data are used, then \code{glm} returns
#'  a \code{dwpm} object with \code{type = data}.
#'
#'  To fit a model, \code{data_DWP} should be a data frame with a row for each carcass
#'  and columns giving (at a minimum) unique carcass IDs, turbine ID, distance
#'  from turbine, and fraction of area searched at the given distance at the
#'  given turbine. Optional columns may include carcass class, covariates that may
#'  influence detection probability (e.g., visibility class), and direction.
#'  If covariates are to be included in the model, then the fraction of area
#'  column would give the fraction of the area in the given covariate level at
#'  that distance. Alternatively, prefab data may be provided in a dataframe,
#'  with structure depending on data type. The simplest case would be that
#'  point estimates only are provided. In that case, if there are no distinctions
#'  among carcass classes (e.g., size), then \code{data_DWP} should be a dataframe
#'  with one column giving the unit (e.g., turbine) and one column with the DWP
#'  at each unit; if distinctions are made among carcass classes, then \code{data_DWP} would
#'  be a data frame with a unit column and a DWP column for each carcass class. If
#'  the DWP estimates incorporate uncertainties, then \code{data_DWP} should be
#'  an array with \code{n_unit * nsim} rows and with colunms for units and DWPs for
#'  each carcass class.
#'
#' @param data_DWP data frame with structure depending on model
#'  type. In general, \code{data_DWP} would be a data frame if a model is to be
#'  fit or if point estimates only are provided as pre-simulated DWP data, and,
#'  if pre-simulated data with variation are provided, then a 2-d array (if one
#'  carcass class) or a list of 2-d arrays (if more than one carcass class). See
#'  "Details" for details.
#'
#' @param type model type may be \code{rings}, \code{glm}, \code{TWL}, or
#'  \code{data}. Currently, only the \code{data} type is supported.
#'
#' @param unitCol name of the column with the units, which must be non-numeric
#'
#' @param dwpCols name(s) of the columns with the DWP data
#'
#' @return an object of an object of class \code{dwpm}, which is a list with
#'  model \code{type} (currently only \code{type = data} is supported) and
#'  \code{model}, which gives the simulated DWP values as an array (if there's
#'  only a single carcass class) or a list of arrays (if there are more than one
#'  carcass classes).
#'
#' @export
#'
dwpm <- function(data_DWP, type = "data", unitCol = NULL, dwpCols = NULL){
  if (is.null(data_DWP)){
    message("data_DWP missing. Assuming DWP = 1.")
    out <- 1
    attr(out, "type") <- "data"
    class(out) <- "dwpm"
    return(out)
  }
  if (!type %in% "data") stop("\"", type, "\" DWP model type not supported.")
  if (!"data.frame" %in% class(data_DWP)) stop("data_DWP must be a data frame")
  if (type %in% "data"){
    numericColumns <- which(unlist(lapply(data_DWP,
        function(x) is.numeric(x) & !any(is.na(x)))))
    # check unitCol
    if (!is.null(unitCol)){ # check format of user-provided unitCol
      if (length(unitCol) > 1)
        stop("unitCol must be the name of a single column in data_DWP")
      if (!unitCol %in% names(data_DWP))
        stop("unitCol (", unitCol, ") not in data_DWP")
      if (is.numeric(unitCol))
        stop("unitCol must be non-numeric")
      unittab <- table(data_DWP[ , unitCol])
      if (length(unique(unittab)) > 1)
        stop("Each unit in unitCol must have the same number of reps")
      nsim <- unittab[1]
    } else { # identify the unitCol when user has not provided ore
      if (length(names(data_DWP)) > length(numericColumns) + 1){
        stop("more than one potential unitCol in data_DWP = ",
          deparse(substitute(data_DWP)), ". A unique unitCol must be provided.")
      } else if (length(names(data_DWP)) == length(numericColumns)) {
         stop("A non-numeric unit column (unitCol) must be present in data_DWP")
      } else {
        unitCol <- names(data_DWP)[-numericColumns]
      }
      unittab <- table(data_DWP[ , unitCol])
      if (length(unique(unittab)) > 1)
        stop("unitCol = NULL but no suitable unit column found in data_DWP ",
            "(all units must have the same number of reps)")
      nsim <- unittab[1]
    }

    unitNames <- unique(data_DWP[ , unitCol])

    # check dwpCols
    if (!is.null(dwpCols)){ # quick check format of user-provided unitCol
      if (any(!dwpCols %in% names(data_DWP))){
        stop("some dwpCols not in data_DWP")
      }
    } else { # identify dwpCols if not provided
      for (i in numericColumns){
        if (sum(data_DWP[ , i] <= 0 | data_DWP[ , i] > 1) == 0)
          dwpCols <- c(dwpCols, names(data_DWP[i]))
      }
      if (length(dwpCols) == 0)
        stop(data_DWP, " must contain a DWP column of values in (0, 1]")
    }

    # construct output structure: matrix if no carcass class distinctions; list otherwise
    if (length(dwpCols) == 1){ # no carcass classes --> matrix
      dwpmat <- array(dim = c(length(unitNames), nsim))
      row.names(dwpmat) <- unitNames
      if (nsim == 1){
        dwpmat[ , 1] <- data_DWP[ , dwpCols]
      } else {
        for (ui in unitNames){
          dwpmat[ui , ] <- data_DWP[data_DWP[ , unitCol] == ui, dwpCols]
        }
      }
      out <- dwpmat
    } else { # more than one carcass class --> list of dwpmat's
      out <- list()
      dwpmat <- array(dim = c(length(unitNames), nsim))
      row.names(dwpmat) <- unitNames
      for (di in dwpCols){
        if (nsim == 1){
          dwpmat[ , 1] <- data_DWP[ , di]
        } else {
          for (ui in unitNames){
            dwpmat[ui , ] <- data_DWP[data_DWP[ , unitCol] == ui, di]
          }
        }
        out[[di]] <- dwpmat
      }
    }
  }
  attr(out, "type") <- type
  class(out) <- "dwpm"
  return(out)
}


#' @title Simulate parameters from a fitted dwp model
#'
#' @description Simulate parameters from a \code{\link{dwpm}} model object
#'
#' @details If the model type = \code{data}, then the number of simulated columns
#'  must be either \code{>=n} (in which case the first n colunms are taken as the
#'  simulated DWP) or 1 (in which case, DWP is assumed constant).
#'
#' @param n the number of simulation draws
#'
#' @param model A \code{\link{dwpm}} object (which is returned from
#'   \code{dwpm()})
#'
#' @return array of \code{n} simulated \code{dwp} values for each unit.
#'  Dimensions = c(n, number of units).
#'
#' @export
#'
rdwp <- function(n, model){
  if (!"dwpm" %in% class(model)) stop("model not of class dwpm")
  if ("data" %in% attr(model, "type")){
    attr(model, "type") <- NULL
    attr(model, "class") <- NULL
    if (!is.list(model)){
      if (NCOL(model) == 1){
        output <- model
      } else if (ncol(model) >= n) {
        output <- model[ , 1:n]
      } else {
        stop("n > number of simulated DWPs")
      }
    } else {
      if (NCOL(model[[1]]) == 1){
        output <- model
      } else if (ncol(model[[1]]) >= n) {
        output <- lapply(model, "[", 1:nrow(model[[1]]), 1:n)
      } else {
        stop("n > number of simulated DWPs")
      }
    }
  }
  class(output) <- "rdwp"
  return(output)
}

