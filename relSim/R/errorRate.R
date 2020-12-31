#' Returns the false positive or false negative rates for a set of IBS and/or
#' KI thresholds
#' 
#' This function is used to calcalate the various tables in the work of Ge et
#' al. and Balding et al. Specifically it can be used to calculate the false
#' positive rate for unrelated pairs being identified as full-sibs or
#' parent-child pairs under differing levels of IBS or KI (or both) thresholds.
#' It can also be used to calculate the false negative rates for full-sib, or
#' parent-child, pairs being identified as unrelated, again with differing
#' levels of IBS, KI or both.
#' 
#' 
#' @param simResults A data.frame with three columns labelled sib, pc and ibs.
#' This will usually be obtained from a call to \code{sim} or
#' \code{readResults}.
#' @param bIBS If \code{TRUE} then IBS thresholds are used to generate the
#' error rates. If both \code{bIBS} and \code{bKI} are \code{TRUE} then both
#' criteria are used.
#' @param bKI If \code{TRUE} then KI thresholds are used to generate the error
#' rates. If both \code{bIBS} and \code{bKI} are \code{TRUE} then both criteria
#' are used.
#' @param rel The relationship used in the simulation. Must be one of 'UN',
#' 'FS' or 'PC'.
#' @param IBSthresh A vector of IBS values that can be used to classify the
#' results as being related (or not).
#' @param KIthresh A vector of KI threshold values that can be used to classify
#' the results as being related (or not).
#' @param nLoci The number of loci being used in the multiplex. This dictates
#' the upper bound on the IBS values.
#' @return A vector (or a two-column matrix) of false negative or false
#' positive rates. If the relationship is 'UN' then false positive rates are
#' returned for parent-child and full-sibs, with parent-child being in column 1
#' and full-sibs in column 2. If the relationship is 'PC' then the false
#' negative rate is returned for parent-child pairs, and if it is 'FS' then the
#' false negative rate for full-sibs.
#' @author James M. Curran
#' @seealso sim, readResults
#' @examples
#' 
#' ## not run
#' \dontrun{data(fbiCaucs)
#' unrel = sim(10000)
#' errorRate(unrel)
#' }
#' 
#' @export errorRate
errorRate = function(simResults, bIBS = TRUE, bKI = FALSE,
                     rel = "UN", IBSthresh = 14:17,
                     KIthresh = c(1e3,1e4,1e5,1e6), nLoci = 13){

    rel = toupper(rel)
    if(!grepl("(UN|FS|PC)", rel)){
        stop("rel must be one of 'UN', 'FS', or 'PC'")
    }

    if(bIBS & !bKI){
        o = tabulate(simResults$ibs+1, nbins = 2*nLoci + 1)
        p = cumsum(o)/sum(o)

        if(rel == "UN"){
            return(1-p[IBSthresh])
        }else{
            return(p[IBSthresh])
        }
    }else if(!bIBS & bKI){
        if(rel == "UN"){
            results = matrix(rep(0, 2*length(KIthresh)), ncol = 2)
            N = nrow(simResults)

            for(i in 1:length(KIthresh)){
                ki = KIthresh[i]
                results[i,1] = sum(simResults$pc >= ki)/N
                results[i,2] = sum(simResults$sib >= ki)/N
            }

            return(results)
        }else if(grepl("(FS|PC)",rel)){
            results = rep(0, length(KIthresh))
            N = nrow(simResults)

            for(i in 1:length(KIthresh)){
                ki = KIthresh[i]
                if(rel == "FS"){
                    results[i] = sum(simResults$sib < ki)/N
                }else{
                    results[i] = sum(simResults$pc < ki)/N
                }
            }

            return(results)
        }
    }else if(bIBS & bKI){ ## both
        ## in this situation IBSthresh and KIBSthresh will have
        ## repeated elements
        cat("Both\n")

        if(length(IBSthresh) != length(KIthresh)){
            stop("IBSthresh and KIthresh must be of equal length when using both criteria")
        }

        if(rel == "UN"){
            results = matrix(rep(0, 2*length(KIthresh)), ncol = 2)
            N = nrow(simResults)

            for(i in 1:length(KIthresh)){
                ki = KIthresh[i]
                ibs = IBSthresh[i]
                results[i,1] = sum(simResults$ibs >= ibs
                                   & simResults$pc >= ki)/N
                results[i,2] = sum(simResults$ibs >= ibs
                                   & simResults$sib >= ki)/N
            }

            return(results)
        }else if(grepl("(FS|PC)",rel)){
            results = rep(0, length(KIthresh))
            N = nrow(simResults)

            for(i in 1:length(KIthresh)){
                ki = KIthresh[i]
                ibs = IBSthresh[i]
                if(rel == "FS"){
                    results[i] = sum(simResults$ibs < ibs
                                     | simResults$sib < ki)/N
                }else{
                    results[i] = sum(simResults$ibs < ibs
                                     | simResults$pc < ki)/N
                }
            }

            return(results)
        }
    }
}


