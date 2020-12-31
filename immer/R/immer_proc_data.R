## File Name: immer_proc_data.R
## File Version: 0.16

immer_proc_data <- function( dat, pid=NULL, rater=NULL, weights=NULL, maxK=NULL)
{
    #-- indicator variable for providing rating data
    is_rating_data <- ! is.null(rater)
    #-- include person IDs if not provided
    N1 <- nrow(dat)
    I <- ncol(dat)
    if (is.null(pid)){
        pid <- 1:N1
    }
    if ( is.null(rater) ){
        rater <- rep(0,N1)
    }
    #-- apply sirt::rm_proc_data() function for processing rating data
    res <- sirt::rm_proc_data( dat=dat, pid=pid, rater=rater, rater_item_int=FALSE, reference_rater=NULL )
    N <- res$N
    res$pid <- res$person.index$pid
    #-- maxK and K
    if ( is.null(maxK) ){
        maxK <- apply( dat, 2, max, na.rm=TRUE)
    }
    dataproc.vars <- res$dataproc.vars
    res$maxK <- maxK[ dataproc.vars$item.index ]
    res$K <- max( maxK )

    #-- create weights if not provided
    if ( ! is.null(weights) ){
        a1 <- stats::aggregate( weights, list(pid), mean )
        weights <- a1[,2]
    }
    if ( is.null(weights) ){
        weights <- rep(1,N)
    }
    res$weights <- weights
    res$W <- sum(weights)
    res$ND <- nrow(res$dat)
    res$I <- I
    res$dataproc.vars$pseudoitems <- colnames(res$dat2)
    res$is_rating_data <- is_rating_data
    if ( ! is_rating_data ){
        res$ND <- NA
        res$RR <- NA
        res$rater <- NULL
        res$dat2.NA <- dat
    }
    #-- pseudoitems_design
    res$pseudoitems_design <- list( dataproc.vars=res$dataproc.vars,
                                person.index=res$person.index, rater.index=res$rater.index,
                                maxK=res$maxK, K=res$K, is_rating_data=res$is_rating_data )

    #-- create item table
    pseudoitems_design <- res$pseudoitems_design
    dfr <- as.data.frame(pseudoitems_design$dataproc.vars)
    itemtable <- data.frame( item=as.factor(dfr$item), rater=as.factor(dfr$rater) )
    itemtable$maxK <- res$maxK
    itemtable$item_name <- colnames(res$dat)[ res$dataproc.vars$item.index[ dfr$item ] ]
    itemtable$rater_name <- paste0( res$rater.index$rater[ dfr$rater ] )
    itemtable$pseudoitem_name <- colnames(res$dat2)
    res$itemtable <- itemtable
    #--- output
    return(res)
}
