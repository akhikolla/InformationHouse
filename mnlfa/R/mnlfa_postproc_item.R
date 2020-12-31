## File Name: mnlfa_postproc_item.R
## File Version: 0.05


mnlfa_postproc_item <- function(parm_list, items, item_type, parms_estimated,
    parms_regularized, parms_regular_types, parms_regular_lam)
{
    item <- NULL
    I <- length(items)
    for (ii in 1:I){
        parm_ii <- parm_list[[ii]]
        parm_ii <- c( parm_ii$b, parm_ii$a )
        item_ii <- paste(items[ii])
        np_ii <- length(parm_ii)
        dfr <- data.frame( itemid=ii, intid=1:np_ii, item=rep(item_ii,np_ii))
        dfr$item_type <- item_type[ii]
        dfr$parm <- names(parm_ii)
        dfr$est <- as.vector(parm_ii)
        dfr$estim <- parms_estimated[[ii]]
        dfr$regul <- parms_regularized[[ii]]
        dfr$reg_type <- parms_regular_types[[ii]]
        dfr$reg_lam <- parms_regular_lam[[ii]]
        rownames(dfr) <- NULL
        item <- rbind( item, dfr )
    }
    return(item)
}
