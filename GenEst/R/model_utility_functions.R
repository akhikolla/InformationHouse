#' @title Combine predictors
#'
#' @description Create a table of combinations of factor levels for the
#'  the predictors in a searcher efficiency or carcass persistence analysis.
#'  This is a utility function called by \code{pkm0} and \code{cpm0}.
#'
#' @param preds Vector of character strings with the names of predictor
#'  variables included in the model.
#'
#' @param data Searcher efficiency or carcass persistence data frame with
#'  columns for each predictor and rows corresponding to carcasses in the field
#'  trials.
#'
#' @return Data frame with columns for each predictor in \code{preds} and rows
#'  for each factor level combination among the predictors. In addition there is
#'  a column with \code{CellNames}, which are the combinations of predictor
#'  levels separated by periods ( . ).
#'
#' @export
#'
combinePreds <- function(preds, data){
  ind <- sapply(data, is.factor)
  data[ind] <- lapply(data[ind], as.character)

  npred <- length(preds)
  if (npred == 0){
    return(data.frame(group = "all", CellNames = "all"))
  } else {
    if(any(is.na(match(preds, names(data))))) {
        stop("At least one predictor missing from data.")
    }
    predLevels <- list()
    for (predi in preds) predLevels[[predi]] <- levels(as.factor(data[[predi]]))
    output <- expand.grid(predLevels, stringsAsFactors = FALSE)
  }
  names(output) <- preds
  output$CellNames <- apply(output, 1, paste0, collapse = ".")
  return(output)
}

#' @title Combine predictors across models
#'
#' @description Create a table of factor combinations of predictors in given
#'  searcher efficiency and carcass persistence models. This is a utility
#'  function called by estgGeneric and governs the cells for which detection
#'  probabilities are calculated.
#'
#' @param preds_CP Character vector with names of carcass persistence predictors.
#'
#' @param preds_SE Character vector with names of searcher efficiency predictors.
#'
#' @param data_CP data frame with columns for each predictor and rows
#'  corresponding to carcasses in the field trials.
#'
#' @param data_SE data frame with columns for each predictor and rows
#'  corresponding to carcasses in the field trials.
#'
#' @return Data frame with columns for each predictor in \code{preds} and rows
#'  for each factor level combination among the predictors. In addition there
#'  are column with \code{CellNames}, \code{CellNames_CP}, and
#'  \code{CellNames_SE}, which are the combinations of predictor levels for all
#'  predictors, CP predictors, and SE predictors (respectively), separated by
#'  periods ( . ).
#'
#' @export
#'
combinePredsAcrossModels <- function(preds_CP, preds_SE, data_CP, data_SE){
  if (length(data_CP) > 0){
    ind <- sapply(data_CP, is.factor)
    data_CP[ind] <- lapply(data_CP[ind], as.character)
  }
  if (length(data_SE) > 0){
    ind <- sapply(data_SE, is.factor)
    data_SE[ind] <- lapply(data_SE[ind], as.character)
  }
  preds <- unique(c(preds_SE, preds_CP))
  npred <- length(preds)
  if (npred == 0){
    preds <- data.frame(group = "all", CellNames = "all",
      CellNames_SE = "all", CellNames_CP = "all")
    return(preds)
  }
  if(any(is.na(match(preds_SE, names(data_SE))))) {
    stop("At least one Searcher Efficiency predictor missing from data.")
  }
  if(any(is.na(match(preds_CP, names(data_CP))))) {
    stop("At least one Carcass Persistence predictor missing from data.")
  }

  predLevels <- list()
  for(predi in preds){
    predLevels_SE <- levels(as.factor(data_SE[[predi]]))
    predLevels_CP <- levels(as.factor(data_CP[[predi]]))
     if (length(predLevels_SE) > 0 & length(predLevels_CP) > 0){
      if (!(all(predLevels_SE %in% predLevels_CP) &
            all(predLevels_CP %in% predLevels_SE))){
        stop("Identical factor has different levels in SE and CP data.")
      }
    }
    predLevels[[predi]] <- unique(c(predLevels_SE, predLevels_CP))
  }
  output <- expand.grid(predLevels, stringsAsFactors = FALSE)
  names(output) <- preds
  output$CellNames <- apply(output, 1, paste0, collapse = ".")
  output_SE <- data.frame(output[ , preds_SE])
  CellNames_SE <- apply(output_SE, 1, paste0, collapse = ".")
  output_CP <- data.frame(output[ , preds_CP])
  CellNames_CP <- apply(output_CP, 1, paste0, collapse = ".")
  if (all(CellNames_SE == "")) CellNames_SE <- "all"
  if (all(CellNames_CP == "")) CellNames_CP <- "all"
  output$CellNames_SE <- CellNames_SE
  output$CellNames_CP <- CellNames_CP
  return(output)
}

#' @title Return the model with the greatest log-likelihood
#'
#' @description  Compares all fitted models in a list and returns the model
#'  with the greatest log-likelihood
#'
#' @param modelSet a list of fitted models with a \code{loglik} element.
#'  Models may be \code{pkm}, \code{cpm}, \code{survreg} objects or any
#'  objects with a \code{loglik} component.
#'
#' @return The model object with the greatest log-likelihood among
#'  the models in \code{modelSet}
#'
#' @export
#'
refMod <- function(modelSet){
  llvec <- sapply(modelSet, "[[", "loglik")
  out <- modelSet[[which(llvec == max(llvec))]]
  return(out)
}

#' @title model utility functions (not exported)
#' @description model utility functions that are not exported
#' @param specific specific model compared against the full set
#' @param modelSet full model set to compare to the specific
#' @name model_utility_functions
NULL

#' @title Check for model components
#'
#' @description Check if all component terms and interactions are
#'   included in a formula. Terms are automatically alphabatized within
#'   interactions.
#'
#' @param formula A formula object
#'
#' @return a logical regarding complete set of terms and interactions
#'
#' @export
#'
checkComponents <- function(formula){

  termOrders <- attr(terms(formula), "order")
  termNames <- attr(terms(formula), "term.labels")
  nterm <- length(termNames)

  if(nterm == 0){
    return(TRUE)
  }
  if(max(termOrders) == 1){
    return(TRUE)
  }

  termNamesAlph <- rep(NA, nterm)
  for (termi in 1:nterm){
    temp <- strsplit(termNames[termi], ":")[[1]]
    temp <- sort(temp)
    termNamesAlph[termi] <- paste(temp, collapse = ":")
  }
  termOrdersSort <- termOrders[order(termOrders, termNamesAlph)]
  termNamesSort <- termNamesAlph[order(termOrders, termNamesAlph)]

  ixns <- which(termOrdersSort > 1)
  nixn <- length(ixns)
  ixnOrders <- termOrdersSort[ixns]
  allThere <- rep(FALSE, nixn)

  for (ixni in 1:nixn){
    ixnName <- termNamesSort[ixns[ixni]]
    ixnParts <- strsplit(ixnName, ":")[[1]]
    ixnSubOrder <- ixnOrders[ixni] - 1
    toCheck <- NULL

    for (orderi in 1:ixnSubOrder){
      combPartsOptions <- combn(ixnParts, orderi)
      addToCheck <- apply(combPartsOptions, 2, paste, collapse = ":")
      toCheck <- c(toCheck, addToCheck)

    }

    allThere[ixni] <- all(toCheck %in% termNamesAlph)
  }
  output <- all(allThere)
  return(output)
}


#' @rdname model_utility_functions
matchCells <- function(specific, modelSet){
  cells_spec <- specific$cells
  cells_set <- modelSetCells(modelSet)
  ncell_set <- nrow(cells_set)
  predictors_spec <- specific$predictors
  predictors_set <- modelSetPredictors(modelSet)
  predictorsMatch_spec <- which(predictors_spec %in% predictors_set)
  predictorsMatch_set <- which(predictors_set %in% predictors_spec)
  if (!all(predictors_spec %in% predictors_set)){
    stop("Complete model set does not include specific model's terms.")
  }

  if (length(predictorsMatch_spec) == 0){
    out <- as.character(rep(cells_spec$CellNames, ncell_set))
  } else{

    matchedPredictors_set <- sort(colnames(cells_set)[predictorsMatch_set])
    matchedPredictors_spec <- sort(colnames(cells_spec)[predictorsMatch_spec])
    cellSubSet_spec <- data.frame(cells_spec[ , matchedPredictors_spec])
    cellSubSet_set <- data.frame(cells_set[ , matchedPredictors_set])
    pasteCellNames_spec <- apply(cellSubSet_spec, 1, paste, collapse = ".")
    pasteCellNames_set <- apply(cellSubSet_set, 1, paste, collapse = ".")

    out <- rep(NA, ncell_set)
    for (celli in 1:ncell_set){
      cellMatch <- which(pasteCellNames_spec == pasteCellNames_set[celli])
      out[celli] <- cells_spec$CellNames[cellMatch]
    }
  }
  return(out)
}
#' @title Trim a Model-Set-Size Complex to a Single Model Per Size
#'
#' @description Select a single model from each carcass class (based on the model
#'   names).
#'
#' @param modSetSize modSetSize complex (cpm or pkm)
#'
#' @param mods named (according to carcass classes) vector of model names to use
#'
#' @return modSetSize reduced to a single model per carcass class
#'
#' @export
#'
trimSetSize <- function(modSetSize, mods){
  sizeclasses <- names(modSetSize)
  nsizeclass <- length(sizeclasses)
  models <- list()
  for (sz in sizeclasses){
    if (! (mods[[sz]] %in% names(modSetSize[[sz]]))){
      stop("requested model for ", sz, " (",  mods[[sz]] , ") not found")
    }
    models[[sz]] <- modSetSize[[sz]][[mods[[sz]]]]
  }
  return(models)
}

#' @rdname model_utility_functions
modelSetModelPredictors <- function(modelSet){
  nmod <- length(modelSet)
  out <- vector("list", length = nmod)
  for (modi in 1:nmod){
     if (!(grepl("Failed model fit", modelSet[[modi]][1]))){
       if ( length(modelSet[[modi]]$predictors) > 0){
         out[[modi]] <- modelSet[[modi]]$predictors
       }
     }
  }
  names(out) <- names(modelSet)
  return(out)
}

#' @rdname model_utility_functions
modelSetPredictors <- function(modelSet){
  unique(unlist(modelSetModelPredictors(modelSet)))
}

#' @rdname model_utility_functions
modelSetModelCells <- function(modelSet){
  nmod <- length(modelSet)
  out <- vector("list", length = nmod)
  for (modi in 1:nmod){
    if (!(grepl("Failed model fit", modelSet[[modi]][1]))){
      out[[modi]] <- modelSet[[modi]]$cells
    }
  }
  names(out) <- names(modelSet)
  return(out)
}

#' @rdname model_utility_functions
modelSetCells <- function(modelSet){
  modelCells <- modelSetModelCells(modelSet)
  modelPreds <- modelSetPredictors(modelSet)

  if (is.null(modelPreds)){
    out <- data.frame("group" = "all", "CellNames" = "all")
    return(out)
  }

  nmod <- length(modelSet)
  npred <- length(modelPreds)
  predLevels <- vector("list", npred)

  for (modi in 1:nmod){
    modiCells <- modelCells[[modi]]
    for (predi in 1:npred){
      if (modelPreds[predi] %in% colnames(modiCells)){
        existing <- predLevels[[predi]]
        new <- as.character(modiCells[ , modelPreds[predi]])
        predLevels[[predi]] <- unique(c(existing, new))
      }
    }
  }
  out <- expand.grid(predLevels, stringsAsFactors = FALSE)
  colnames(out) <- modelPreds
  out$CellNames <- apply(out, 1, paste, collapse = ".")
  return(out)
}

