## File Name: immer_jml_proc_shortcut.R
## File Version: 0.14

immer_jml_proc_shortcut <- function(dat, pid, shortcut, weights)
{
    N <- nrow(dat)
    if (is.null(pid)){
        pid <- 1:N
    }
    if (is.null(weights)){
        weights <- rep(1,N)
    }
    dfr <- data.frame( orig=1:N, pid=pid, weights=weights)

    #-- shortcut
    if (shortcut){
        res <- TAM::tam_NA_pattern(x=dat)
        mp.index <- res$mp.index
        dfr$mp_index <- mp.index
        dfr$score <- rowSums( dat, na.rm=TRUE )
        dfr$index <- dfr$mp_index*10^round( log10( max(dfr$score) ) + 2 ) + dfr$score
        dfr <- dfr[ order(dfr$index), ]
        dfr$update <- c(1,diff(dfr$index) > 0 )
        a1 <- rowsum(weights, dfr$index)
        dfr$update_weights <- a1[ match( dfr$index, rownames(a1) ),1] * dfr$update
        dat <- dat[ dfr$orig, ]
        pid <- pid[ dfr$orig ]
        weights <- weights[ dfr$orig ]
        rownames(dat) <- NULL
    }
    #-- no shortcut
    if (! shortcut){
        dfr$update <- 1
        dfr$update_weights <- 1
    }
    #--- output
    res <- list( dat=dat, pid=pid, shortcut=shortcut, N=N, shortcut_index=dfr, weights=weights)
    return(res)
}
