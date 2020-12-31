#' Modify the function name prefix
#'
#' Add some text to the function name so it fits in with the error message.
#'
#' @param functionName The function name.
#'
#' @keywords internal
modifyFunctionPrefix <- function(functionName){
    return( paste0("in ", functionName, " :\n    ") )
}



#' Modify the string containing the value
#'
#' Create a string with the value to be used in an error message.
#'
#' @param value The value of the offending parameter.
#'
#' @keywords internal
modifyValueStr <- function(value){
    return( paste0(" (value=", value, ")") )
}



#' Create the error message for a certain parameter that is missing
#'
#' For a particular parameter, create an error message specifying that
#' the parameter is missing.
#'
#' @param paramName Name of parameter, e.g. \code{lambda}
#'
#' @param algoName Name of algorithm, e.g. \code{AFF change detector}
#'
#' @param functionName The name of the original function which is being
#'                      called to check the parameters. Will help user
#'                      with debugging.
#'
#' @keywords internal
getMissingErrorMessage <- function(paramName, algoName=NULL, functionName=""){
    algoString <- ""
    if (is.null(algoName)){
        #do nothing
    } else {
        algoString <- paste0("of ", algoName, " ")
    }
    if (functionName != ""){
        functionName <- modifyFunctionPrefix(functionName)
    }
    missingMessage <- paste0(functionName, paramName, " is missing - initialisation ", algoString, "failed.")
    return(missingMessage)
}



#' Create the error message for a certain parameter that is not finite
#'
#' For a particular parameter, create an error message specifying that
#' the parameter is not finite (i.e. not numeric and not NaN)
#'
#' @param paramName Name of parameter, e.g. \code{lambda}
#'
#' @param value The value of the non-finite parameter
#'
#' @param algoName Name of algorithm, e.g. \code{AFF change detector}
#'
#' @param functionName The name of the original function which is being
#'                      called to check the parameters. Will help user
#'                      with debugging.
#'
#' @keywords internal
getFiniteErrorMessage <- function(paramName, value=NULL, algoName=NULL, functionName=""){
    algoString <- ""
    if (is.null(algoName)){
        #do nothing
    } else {
        algoString <- paste0("of ", algoName, " ")
    }
    if (functionName != ""){
        functionName <- modifyFunctionPrefix(functionName)
    }
    valueStr <- ""
    if ( !(is.null(value)) ){
        valueStr <- modifyValueStr(value)
    }
    finiteMessage <- paste0(functionName, paramName, valueStr, " is not a finite numeric value - initialisation ", algoString, "failed.")
    return(finiteMessage)
}



#' Create the error message for a certain parameter outside its bounds
#'
#' For a particular parameter, create an error message specifying that
#' the parameter is outside its bounds.
#'
#' @param paramName Name of parameter, e.g. \code{lambda}
#'
#' @param value The value of the non-finite parameter
#'
#' @param lowerBound Lower bound of range, e.g. \code{0}.
#'
#' @param upperBound Upper bound of range, e.g. \code{1}.
#'
#' @param algoName Name of algorithm, e.g. \code{AFF change detector}
#'
#' @param strictly Use "[" or "(" \code{strictly=FALSE} means "[" and  
#'        \code{strictly=TRUE} means "("
#'
#' @param functionName The name of the original function which is being
#'                      called to check the parameters. Will help user
#'                      with debugging. Default is \code{""}, in which
#'                      case no function name will be displayed.
#'
#' @keywords internal
getInBoundsErrorMessage <- function(paramName, value=NULL, lowerBound, upperBound, algoName=NULL, strictly=FALSE, functionName=""){
    algoString <- ""
    if (is.null(algoName)){
        #do nothing
    } else {
        algoString <- paste0("of ", algoName, " ")
    }
    openB <- "["
    closedB <- "]"
    if (strictly){
        openB <- "("
        closedB <- ")"
    }
    if (functionName != ""){
        functionName <- modifyFunctionPrefix(functionName)
    }
    valueStr <- ""
    if ( !(is.null(value)) ){
        valueStr <- modifyValueStr(value)
    }
    inBoundsMessage <- paste0(functionName, paramName, valueStr, " is not in range ", openB, lowerBound, ", ", upperBound, closedB, " - initialisation ", algoString, "failed.")
    return(inBoundsMessage)
}


#' Check if an object exists
#'
#' For a particular parameter, check if the object with that name exists
#'
#' @param object The object.
#'
#' @keywords internal
doesNotExist <- function(object){
    #From Brian Ripley:
    #https://stat.ethz.ch/pipermail/r-help/2002-January/018277.html
    if ( exists(as.character(substitute(object))) ){
#        cat("EXISTS value: ", object, "\n")
        return (FALSE)
    }
#        cat("DOES NOT exist value: ", object, "\n")
    return(TRUE)
}



#' Create the error message for a certain parameter that does not exist
#'
#' For a particular parameter, create an error message specifying that
#' the parameter does not exist (e.g. \code{Na}, which is not the same
#' as \code{NA}).
#'
#' @param paramName Name of parameter, e.g. \code{lambda}
#'
#' @param value The name of the object that does not exist.
#'
#' @param algoName Name of algorithm, e.g. \code{AFF change detector}
#'
#' @param functionName The name of the original function which is being
#'                      called to check the parameters. Will help user
#'                      with debugging.
#'
#' @keywords internal
getExistsErrorMessage <- function(paramName, value=NULL, algoName=NULL, functionName=""){
    algoString <- ""
    if (is.null(algoName)){
        #do nothing
    } else {
        algoString <- paste0("of ", algoName, " ")
    }
    if (functionName != ""){
        functionName <- modifyFunctionPrefix(functionName)
    }
    valueStr <- ""
    if ( !(is.null(value)) ){
#        valueStr <- modifyValueStr( as.character(substitute(value)) )
        valueStr <- paste0(" ", value)
    }
    finiteMessage <- paste0(functionName, "value used for ", paramName, " does not exist - initialisation ", algoString, "failed.")
    return(finiteMessage)
}




#' Create the error message for a certain parameter not above a threshold
#'
#' For a particular parameter, create an error message specifying that
#' the parameter is not above a minimum threshold
#'
#' @param paramName Name of parameter, e.g. \code{lambda}
#'
#' @param value The value of the non-finite parameter
#'
#' @param lowerBound Lower bound of range, e.g. \code{0}.
#'
#' @param algoName Name of algorithm, e.g. \code{AFF change detector}
#'
#' @param functionName The name of the original function which is being
#'                      called to check the parameters. Will help user
#'                      with debugging. Default is \code{""}, in which
#'                      case no function name will be displayed.
#'
#' @keywords internal
getIsAboveErrorMessage <- function(paramName, value=NULL, lowerBound, algoName=NULL, functionName=""){
    algoString <- ""
    if (is.null(algoName)){
        #do nothing
    } else {
        algoString <- paste0("of ", algoName, " ")
    }
    if (functionName != ""){
        functionName <- modifyFunctionPrefix(functionName)
    }
    valueStr <- ""
    if ( !(is.null(value)) ){
        valueStr <- modifyValueStr(value)
    }
    isAboveMessage <- paste0(functionName, paramName, value, " is not greater than ", lowerBound, " - initialisation ", algoString, "failed.")
    return(isAboveMessage)
}




#' Check the arguments for FFF (no change detection)
#'
#' A function which will do all the checks for FFF, check the parameters
#' are in the proper ranges, etc. Does not return any values.
#'
#' @param lambda The value of the FFF. Must be finite, in range [0,1], etc.
#'
#' @param functionName The name of the original function which is being
#'                      called to check the parameters. Will help user
#'                      with debugging. Default is \code{""}, in which
#'                      case no function name will be displayed.
#'
#' @keywords internal
checkFFFargs <- function(lambda, functionName=""){
    LAMBDAMIN <- 0
    LAMBDAMAX <- 1
#    if ( doesNotExist(lambda) ){
#        stop(getExistsErrorMessage(paramName="lambda", algoName="AFF", functionName=functionName), call.=FALSE)
#    }

    #check it is numeric, else warning
    #is.finite better than is.numeric, because is.numeric return TRUE for NaN
    if (is.finite(lambda)){ 
        if ( isInBounds(lambda, LAMBDAMIN, LAMBDAMAX) ){
            #nothing
        }
        else{
            stop(getInBoundsErrorMessage(paramName="lambda", value=lambda, lowerBound=LAMBDAMIN, upperBound=LAMBDAMAX, algoName="FFF", functionName=functionName), call.=FALSE)
        }
    } else {
        #if is NOT finite
        stop(getFiniteErrorMessage(paramName="lambda", value=lambda, algoName="FFF", functionName=functionName), call.=FALSE)
    }
    
    #no need to return anything.
}





#' Check the arguments for AFF (no change detection)
#'
#' A function which will do all the checks for AFF, check the parameters
#' are in the proper ranges, etc. Does not return any values.
#'
#' @param eta The value of the step size in the gradient descent. 
#'            Must be finite, in range [0,1], etc.
#'
#' @param functionName The name of the original function which is being
#'                      called to check the parameters. Will help user
#'                      with debugging. Default is \code{""}, in which
#'                      case no function name will be displayed.
#'
#' @keywords internal
checkAFFargs <- function(eta, functionName=""){
    ETAMIN <- 0
    ETAMAX <- 1

#    if ( doesNotExist(eta) ){
#        stop(getExistsErrorMessage(paramName="eta", algoName="AFF", functionName=functionName), call.=FALSE)
#    }

    #check it is numeric, else warning
    #is.finite better than is.numeric, because is.numeric return TRUE for NaN
    if (is.finite(eta)){ 
        if ( isInBounds(eta, ETAMIN, ETAMAX) ){
            #nothing!
        }
        else{
            stop(getInBoundsErrorMessage(paramName="eta", value=eta, lowerBound=ETAMIN, upperBound=ETAMAX, algoName="AFF", functionName=functionName), call.=FALSE)
        }
    } else {
        #eta is NOT finite
        stop(getFiniteErrorMessage(paramName="eta", value=eta, algoName="AFF", functionName=functionName), call.=FALSE)
    }
    
    #no need to return anything.
}



#' Check the arguments for FFF change detection initialisation
#'
#' A function which will do all the checks for FFFcd, check the parameters
#' are in the proper ranges, etc. Does not return any values.
#'
#' @param lambda The value of the FFF. Must be finite, in range [0,1], etc.
#'
#' @param alpha The value of the threshold. Must be finite, in (0,1), etc.
#'
#' @param BL The value of the burn-in length. Must be finite, above 0, etc.
#'
#' @param functionName The name of the original function which is being
#'                      called to check the parameters. Will help user
#'                      with debugging. Default is \code{""}, in which
#'                      case no function name will be displayed.
#'
#' @return The truncated burn-in length (since it must be an integer).
#'
#' @keywords internal
checkFFFMeanCDargs <- function(alpha, lambda, BL, functionName=""){
    LAMBDAMIN <- 0
    LAMBDAMAX <- 1
    BLMIN <- 0
    ALPHAMIN <- 0
    ALPHAMAX <- 1


    #lambda cannot be missing
    if (missing(lambda)){
        stop(getMissingErrorMessage(paramName="lambda", algoName="FFF change detector", functionName=functionName), call.=FALSE)
    }
#    if ( doesNotExist(lambda) ){
#        stop(getExistsErrorMessage(paramName="lambda", algoName="FFF change detector", functionName=functionName), call.=FALSE)
#    }
    #check it is numeric, else warning
    #is.finite better than is.numeric, because is.numeric return TRUE for NaN
    if (! (is.finite(lambda))){
        stop(getFiniteErrorMessage(paramName="lambda", algoName="FFF change detector", functionName=functionName), call.=FALSE)
    }
    if ( !(isInBounds(lambda, LAMBDAMIN, LAMBDAMAX)) ){
        stop(getInBoundsErrorMessage(paramName="lambda", lowerBound=LAMBDAMIN, upperBound=LAMBDAMAX, algoName="FFF change detector", functionName=functionName), call.=FALSE)
    }

    #alpha cannot be missing
    if (missing(alpha)){
        stop(getMissingErrorMessage(paramName="alpha", algoName="FFF change detector", functionName=functionName), call.=FALSE)
    }
#    if ( doesNotExist(alpha) ){
#        stop(getExistsErrorMessage(paramName="alpha", algoName="FFF change detector", functionName=functionName), call.=FALSE)
#    }
    #check it is numeric, else warning
    #is.finite better than is.numeric, because is.numeric return TRUE for NaN
    if (! (is.finite(alpha))){
        stop(getFiniteErrorMessage(paramName="alpha", algoName="FFF change detector", functionName=functionName), call.=FALSE)
    }
    #alpha must be strictly in range (0, 1)
    if ( !(isInBoundsStrictly(alpha, ALPHAMIN, ALPHAMAX)) ){
        stop(getInBoundsErrorMessage(paramName="alpha", lowerBound=ALPHAMIN, upperBound=ALPHAMAX, algoName="FFF change detector", strictly=TRUE, functionName=functionName), call.=FALSE)
    }
    
    #BL must be finite
#    if ( doesNotExist(BL) ){
#        stop(getExistsErrorMessage(paramName="BL", algoName="FFF change detector", functionName=functionName), call.=FALSE)
#    }
    if (! (is.finite(BL)) ){
        stop(getFiniteErrorMessage(paramName="BL", algoName="FFF change detector", functionName=functionName), call.=FALSE)
    }
    #BL must be above min, and truncated
    if ( !(isAboveBound(BL, BLMIN)) ){
        stop(getIsAboveErrorMessage(paramName="BL", lowerBound=BLMIN, algoName="FFF change detector", functionName=functionName), call.=FALSE)
    }
    #BL must be an integer - so force truncation
    BL <- trunc(BL)

    return(BL)

    #no need to return anything.
}



#' Check the arguments for AFF mean change detector
#'
#' A function which will do all the checks for AFF, check the parameters
#' are in the proper ranges, etc. Does not return any values.
#'
#' @param alpha The value of the threshold. Must be finite, in (0,1), etc.
#'
#' @param eta The value of the step size in the gradient descent. 
#'            Must be finite, in range [0,1], etc.
#'
#' @param BL The value of the burn-in length. Must be finite, above 0, etc.
#'
#' @param functionName The name of the original function which is being
#'                      called to check the parameters. Will help user
#'                      with debugging. Default is \code{""}, in which
#'                      case no function name will be displayed.
#'
#' @keywords internal
checkAFFMeanCDargs <- function(alpha, eta, BL, functionName=""){
    BLMIN <- 2
    BL_ZERO <- 0
    ALPHAMIN <- 0
    ALPHAMAX <- 1
    ETAMIN <- 0
    ETAMAX <- 1

    #alpha cannot be missing
    if (missing(alpha)){
        stop(getMissingErrorMessage(paramName="alpha", algoName="AFF change detector", functionName=functionName), call.=FALSE)
    }
#    if ( doesNotExist(alpha) ){
#        stop(getExistsErrorMessage(paramName="alpha", algoName="AFF change detector", functionName=functionName), call.=FALSE)
#    }
    #check it is numeric, else warning
    #is.finite better than is.numeric, because is.numeric return TRUE for NaN
    if (! (is.finite(alpha))){
        stop(getFiniteErrorMessage(paramName="alpha", algoName="AFF change detector", functionName=functionName), call.=FALSE)
    }
    #alpha must be strictly in range (0, 1)
    if ( !(isInBoundsStrictly(alpha, ALPHAMIN, ALPHAMAX)) ){
        stop(getInBoundsErrorMessage(paramName="alpha", lowerBound=ALPHAMIN, upperBound=ALPHAMAX, algoName="AFF change detector", strictly=TRUE, functionName=functionName), call.=FALSE)
    }

#    if ( doesNotExist(eta) ){
#        stop(getExistsErrorMessage(paramName="eta", algoName="AFF change detector", functionName=functionName), call.=FALSE)
#    }
    #check it is numeric, else warning
    #is.finite better than is.numeric, because is.numeric return TRUE for NaN
    if (is.finite(eta)){ 
        if ( isInBounds(eta, ETAMIN, ETAMAX) ){
            #nothing!
        }
        else{
            etaRangeMessage <- paste("eta is not in range [", ETAMIN, ", ", ETAMAX, "] - initialisation failed.", sep="")
        stop(getInBoundsErrorMessage(paramName="eta", lowerBound=ETAMIN, upperBound=ETAMAX, algoName="AFF change detector", functionName=functionName), call.=FALSE)
        }
    } else {
        stop(getFiniteErrorMessage(paramName="eta", algoName="AFF change detector", functionName=functionName), call.=FALSE)
    }

    #BL must be finite
    if (! (is.finite(BL)) ){
        stop(getFiniteErrorMessage(paramName="BL", algoName="AFF change detector", functionName=functionName), call.=FALSE)
    }
    #BL must be above min, and truncated
    if ( !(isAboveBound(BL, BLMIN)) && !(BL == BL_ZERO) ){
        stop(getIsAboveErrorMessage(paramName="BL", lowerBound=BLMIN, algoName="AFF change detector", functionName=functionName), call.=FALSE)
    }
    #BL must be an integer - so force truncation
    BL <- trunc(BL)
    return(BL)
}






#' Check the arguments for CUSUM mean change detector
#'
#' A function which will do all the checks for AFF, check the parameters
#' are in the proper ranges, etc. Does not return any values.
#'
#' @param k One of the CUSUM control parameters. Must be finite, etc.
#'
#' @param h One of the CUSUM control parameters. Must be finite, etc.
#'
#' @param BL The value of the burn-in length. Must be finite, above 0, etc.
#'
#' @param functionName The name of the original function which is being
#'                      called to check the parameters. Will help user
#'                      with debugging. Default is \code{""}, in which
#'                      case no function name will be displayed.
#'
#' @keywords internal
checkCUSUMMeanCDargs <- function(k, h, BL, functionName=""){
    BLMIN <- 0
    #check it is numeric, else warning
    #is.finite better than is.numeric, because is.numeric return TRUE for NaN
    if (! (is.finite(k))){
        stop(getFiniteErrorMessage(paramName="k", algoName="CUSUM change detector", functionName=functionName), call.=FALSE)
    }
    if (! (is.finite(h))){
        stop(getFiniteErrorMessage(paramName="h", algoName="CUSUM change detector", functionName=functionName), call.=FALSE)
    }
    #BL must be finite
    if (! (is.finite(BL)) ){
        stop(getFiniteErrorMessage(paramName="BL", algoName="CUSUM change detector", functionName=functionName), call.=FALSE)
    }
    #BL must be above min, and truncated
    if ( !(isAboveBound(BL, BLMIN)) ){
        stop(getIsAboveErrorMessage(paramName="BL", lowerBound=BLMIN, algoName="CUSUM change detector",functionName=functionName), call.=FALSE)
    }
    #BL must be an integer - so force truncation
    BL <- trunc(BL)
    return(BL)
}




#' Check the arguments for EWMA mean change detector
#'
#' A function which will do all the checks for AFF, check the parameters
#' are in the proper ranges, etc. Does not return any values.
#'
#' @param r One of the EWMA control parameters. Must be finite, etc.
#'
#' @param L One of the EWMA control parameters. Must be finite, etc.
#'
#' @param BL The value of the burn-in length. Must be finite, above 0, etc.
#'
#' @param functionName The name of the original function which is being
#'                      called to check the parameters. Will help user
#'                      with debugging. Default is \code{""}, in which
#'                      case no function name will be displayed.
#'
#' @keywords internal
checkEWMAMeanCDargs <- function(r, L, BL, functionName=""){
    BLMIN <- 0
    RMIN <- 0
    RMAX <- 1

    #check it is numeric, else warning
    #is.finite better than is.numeric, because is.numeric return TRUE for NaN
    if (is.finite(r)){
        if ( isInBounds(r, RMIN, RMAX) ){
            #nothing
        }
        else{
            stop(getInBoundsErrorMessage(paramName="r", value=r, lowerBound=RMIN, upperBound=RMAX, algoName="EWMA", functionName=functionName), call.=FALSE)
        }
    } else{
        stop(getFiniteErrorMessage(paramName="r", algoName="EWMA change detector", functionName=functionName), call.=FALSE)
    }

    if (! (is.finite(L))){
        stop(getFiniteErrorMessage(paramName="L", algoName="EWMA change detector", functionName=functionName), call.=FALSE)
    }
    #BL must be finite
    if (! (is.finite(BL)) ){
        stop(getFiniteErrorMessage(paramName="BL", algoName="EWMA change detector", functionName=functionName), call.=FALSE)
    }
    #BL must be above min, and truncated
    if ( !(isAboveBound(BL, BLMIN)) ){
        stop(getIsAboveErrorMessage(paramName="BL", lowerBound=BLMIN, algoName="EWMA change detector",functionName=functionName), call.=FALSE)
    }
    #BL must be an integer - so force truncation
    BL <- trunc(BL)
    return(BL)
}



#' Check the prechange mean, sigma and variance are appropriately defined
#'
#' If \code{usePrechange} was set to \code{TRUE}, then \code{prechangeMean}
#' must not be \code{NULL}, and at least one of \code{prechangeSigma} and
#' \code{prechangeVar} must not be \code{NULL}.
#'
#' @param prechangeMean The prechange mean. Should not be \code{NULL}.
#'
#' @param prechangeSigma The prechange standard deviation. 
#'                       Either it or \code{prechangeVar} must not 
#'                       be \code{NULL}.
#'
#' @param prechangeVar The prechange variance.
#'                     Either it or \code{prechangeSigma} must not 
#'                     be \code{NULL}.
#'
#' @param functionName The name of the original function which is being
#'                      called to check the parameters. Will help user
#'                      with debugging. Default is \code{""}, in which
#'                      case no function name will be displayed.
#'
#' @keywords internal
checkPrechange <- function(prechangeMean, prechangeSigma, prechangeVar, functionName=""){

    if (is.null(prechangeMean)){
        message <- paste0(functionName, " - usePrechange is TRUE, but prechangeMean is NULL - initialisation failed.")
        stop(message)
    }

    if ( !(is.null(prechangeMean)) & !(is.finite(prechangeMean)) ){
        stop(getFiniteErrorMessage(paramName="prechangeMean", functionName=functionName))
    }

    if ( (is.null(prechangeVar)) & (is.null(prechangeSigma)) ){
        message <- paste0(functionName, " - usePrechange is TRUE, but both prechangeSigma and prechangeVar are NULL - initialisation failed.")
        stop(message)
    }

    if ( !(is.null(prechangeSigma)) & !(is.finite(prechangeMean)) ){
        stop(getFiniteErrorMessage(paramName="prechangeSigma", functionName=functionName))
    }
    if ( !(is.null(prechangeVar)) & !(is.finite(prechangeMean)) ){
        stop(getFiniteErrorMessage(paramName="prechangeVar", functionName=functionName))
    }
}


#' Check the multiple and single flags, as well as usePrechange and skipCheck
#'
#' Make sure the booleans \code{multiple} and \code{single} are not both 
#' \code{NULL}. Also 
#'
#' @param multiple A boolean
#'
#' @param single A boolean
#'
#' @param usePrechange A boolean
#'
#' @param skipCheck A boolean
#'
#' @param functionName The name of the original function which is being
#'                      called to check the parameters. Will help user
#'                      with debugging. Default is \code{""}, in which
#'                      case no function name will be displayed.
#'
#' @keywords internal
checkBooleans <- function(multiple, single, usePrechange, skipCheck, functionName=""){

    if ( !(is.logical(multiple)) ) {
        message <- paste0(functionName, " - 'multiple' is not a boolean - initialisation failed.")
        stop(message)
    }

    if ( !(is.logical(single)) ) {
        message <- paste0(functionName, " - 'single' is not a boolean - initialisation failed.")
        stop(message)
    }


    #check if usePrechange is not null
    if (is.null(usePrechange)==FALSE){
        #then it must be a boolean
        if ( !(is.logical(usePrechange)) ) {
            #if not, stop
            message <- paste0(functionName, " - 'usePrechange' is NULL and also not a booleans - initialisation failed.")
            stop(message)
        }
    }

    if ( !(is.logical(skipCheck)) ) {
        message <- paste0(functionName, " - 'skipCheck' is not a boolean - initialisation failed.")
        stop(message)
    }
}



#' Check the stream consists of finite numeric values
#'
#' Check the stream does not have any NULLs, NAs, infs, etc
#'
#' @param x The stream
#'
#' @param functionName The name of the original function which is being
#'                      called to check the parameters. Will help user
#'                      with debugging. Default is \code{""}, in which
#'                      case no function name will be displayed.
#'
#' @keywords internal
checkStream <- function(x, functionName=""){
    if (is.null(x)){
        message <- paste0(functionName, " the stream 'x' is NULL - initialisation failed.")
        stop(message)
    }

    if (any(is.na(x))){
        message <- paste0(functionName, " the stream 'x' contains NA values - initialisation failed.")
        stop(message)
    }

    if (any(is.nan(x))){
        message <- paste0(functionName, " the stream 'x' contains NaN values - initialisation failed.")
        stop(message)
    }

    if (any(is.infinite(x))){
        message <- paste0(functionName, " the stream 'x' contains non-finite values - initialisation failed.")
        stop(message)
    }
}


#' Check skipChecks
#'
#' Check if the skipChecks parameter is a boolean (it should be).
#'
#' @param skipChecks A boolean
#'
#' @param functionName The name of the original function which is being
#'                      called to check the parameters. Will help user
#'                      with debugging. Default is \code{""}, in which
#'                      case no function name will be displayed.
#'
#' @keywords internal
checkStream <- function(x, functionName=""){
    if (is.null(x)){
        message <- paste0(functionName, " the stream 'x' is NULL - initialisation failed.")
        stop(message)
    }

    if (any(is.na(x))){
        message <- paste0(functionName, " the stream 'x' contains NA values - initialisation failed.")
        stop(message)
    }

    if (any(is.nan(x))){
        message <- paste0(functionName, " the stream 'x' contains NaN values - initialisation failed.")
        stop(message)
    }

    if (any(is.infinite(x))){
        message <- paste0(functionName, " the stream 'x' contains non-finite values - initialisation failed.")
        stop(message)
    }
}

