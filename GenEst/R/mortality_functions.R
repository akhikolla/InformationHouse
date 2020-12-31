#' @title Estimate mortality
#'
#' @description Given given fitted Searcher Efficiency and Carcass 
#'   Persistence models; Search Schedule, Density Weighted Proportion,
#'   and Carcass Observation data; and information about the fraction of the
#'   the facility that was surveyed. 
#'
#' @param data_CO Carcass Observation data
#'
#' @param data_SS Search Schedule data
#'
#' @param data_DWP Survey unit (rows) by carcass class (columns) density weighted
#'   proportion table
#'
#' @param frac fraction carcasses on ground that was surveyed but not accounted
#'  for in DWP
#'
#' @param COdate Column name for the date found data
#'
#' @param model_SE Searcher Efficiency model (or list of models if there are
#'   multiple carcass classes)
#'
#' @param model_CP Carcass Persistence model (or list of models if there are
#'   multiple carcass classes)
#'
#' @param model_DWP fitted dwp model (optional)
#'
#' @param unitCol Column name for the unit indicator (optional)
#'
#' @param SSdate Column name for the date searched data
#'
#' @param sizeCol Name of colum in \code{data_CO} where the carcass classes
#'  are recorded. Optional. If none provided, it is assumed there is no
#'  distinctions among carcass classes.
#'
#' @param IDcol column with unique carcass (CO) identifier
#'
#' @param DWPCol Column name for the DWP values in the DWP table when no
#'   carcass class is used and there is more than one column in \code{data_DWP}
#'   that could be interpreted as DWP.
#'
#' @param nsim the number of simulation draws
#'
#' @param max_intervals maximum number of arrival intervals to consider
#'  for each carcass
#'
#' @return list of Mhat, Aj, ghat, DWP (by carcass), and Xtot = total number of
#'  carcasses observe
#'
#' @examples 
#'  \dontrun{
#'  data(mock)
#'  model_SE <- pkm(formula_p = p ~ HabitatType, formula_k = k ~ 1,
#'               data = mock$SE
#'              )
#'  model_CP <- cpm(formula_l = l ~ Visibility, formula_s = s ~ Visibility, 
#'                data = mock$CP, dist = "weibull",
#'                left = "LastPresentDecimalDays", 
#'                right = "FirstAbsentDecimalDays"
#'              )
#'  eM <- estM(nsim = 1000, data_CO = mock$CO, data_SS = mock$SS, 
#'          data_DWP = mock$DWP, frac = 1, model_SE = model_SE, 
#'          model_CP = model_CP, COdate = "DateFound",
#'          DWPCol = "S", sizeCol = NULL
#'        )
#'  }
#'
#' @export 
#'
estM <- function(data_CO, data_SS, data_DWP = NULL, frac = 1,
                 COdate = "DateFound", model_SE, model_CP, model_DWP = NULL,
                 unitCol = NULL, SSdate = NULL, sizeCol = NULL, IDcol = NULL,
                 DWPCol = NULL, nsim = 1000, max_intervals = 8){

  i <- sapply(data_CO, is.factor)
  data_CO[i] <- lapply(data_CO[i], as.character)
  i <- sapply(data_SS, is.factor)
  data_SS[i] <- lapply(data_SS[i], as.character)
  if (!is.null(data_DWP) && "data.frame" %in% class(data_DWP)){
    i <- sapply(data_DWP, is.factor)
    data_DWP[i] <- lapply(data_DWP[i], as.character)
  }
  if (!(COdate %in% colnames(data_CO))){
    stop("COdate not found in data_CO")
  }
  data_CO[ , COdate] <- checkDate(data_CO[ , COdate])
  if (is.null(data_CO[ , COdate]))
    stop("dates_CO not unambiguously intepretable as dates")

  if (is.null(unitCol))
    unitCol <- defineUnitCol(data_CO = data_CO, data_DWP = data_DWP)
  # if no sizeCol is provided, then the later analysis is done without
  #   making distinctions between carcass classes; no error-checking here
  # if sizeCol is provided, it must be present in CO. Its levels must
  #  also all be present in DWP, but the check is done in the DWPbyCarcass 
  #  function, which allow DWPbyCarcass to more readily be used as a 
  #  standalone function if user wishes.
  if (!is.null(sizeCol)){
    if (!(sizeCol %in% colnames(data_CO))){
      stop("carcass class column not in carcass data.")
    } else if (!all(data_CO[, sizeCol] %in% names(data_DWP))){
      stop("a carcass class in data_CO is missing from data_DWP")
    } else {
      sizeclasses <- sort(as.character(unique(data_CO[ , sizeCol])))
      nsizeclass <- length(sizeclasses)
    }
  }
  # error-checking for match b/t DWP and CO data is done in DWPbyCarcass
  if (!is.null(data_DWP) && !is.null(model_DWP))
    stop("provide either data_DWP or model_DWP, not both")
  if (is.null(data_DWP) && is.null(model_DWP)){
    model_DWP <- dwpm(data_DWP = NULL)
  }
  if (!is.null(data_DWP)){
      model_DWP <- dwpm(data_DWP = data_DWP, type = "data", unitCol = unitCol,
          dwpCols = DWPCol)
  }
  if (!"dwpm" %in% class(model_DWP)) stop("fitted DWP model required")

  est <- estg(data_CO = data_CO, COdate = COdate,
      data_SS = data_SS, SSdate = SSdate,
      model_SE = model_SE, model_CP = model_CP, model_DWP = model_DWP,
      unitCol = unitCol, sizeCol = sizeCol, IDcol = IDcol,
      nsim = nsim, max_intervals = max_intervals)
  gDf <- est$ghat * est$DWP * frac

  c_out <- which(rowSums(gDf) == 0)
  if (length(c_out) == 0){
    n <- length(gDf)
    Mhat <- (rbinom(n, round(1/gDf), gDf) - (round(1/gDf)*gDf - 1))/gDf
  } else {
    Mhat <- array(0, dim = c(dim(data_CO)[1], nsim))
    gDf <- gDf[-c_out, ]
    n <- length(gDf)
    Mhat[-c_out, ] <- (rbinom(n, round(1/gDf), gDf) - (round(1/gDf)*gDf - 1))/gDf
  }
  row.names(Mhat) <- row.names(est$ghat)
  out <- list(Mhat = Mhat, Aj = est$Aj, ghat = est$ghat, DWP = est$DWP,
    Xtot = nrow(data_CO) - length(c_out))
  class(out) <- c("estM", "list")
  return(out)
}

#' @title Summarize total mortality estimation
#'
#' @description \code{summary} defined for class \code{estM} objects
#'
#' @param object \code{estM} object
#'
#' @param ... arguments to pass down
#'
#' @param CL confidence level
#'
#' @export
#'
summary.estM <- function(object, ..., CL = 0.90){
  alpha <- 1 - CL
  Mtot <- colSums(object$Mhat) 
  Mtot[Mtot < object$Xtot] <- object$Xtot
  out <- round(quantile(Mtot, probs = c(0.5, alpha/2, 1 - alpha/2)), 2)
  names(out) <- c("median", paste0(100*c(alpha/2, 1- alpha/2), "%"))
  return(out)
}
