## File Name: MHprop_refresh.R
## File Version: 0.14


MHprop_refresh <- function( MHprop )
{
    vars <- MHprop$VARS_refreshing
    V <- length(vars)
    for (vv in 1:V){
        var.vv <- vars[vv]
        ri <- MHprop$refresh_iter[[ var.vv ]]
        accept <- MHprop$accept[[ var.vv ]] / MHprop$refresh_iter[[ var.vv ]]
        SD <- MHprop$SD[[ var.vv ]]
        SDchange <- MHprop$refresh_SDchange[[ var.vv ]]
        MHprop$SD[[ var.vv ]] <- MHprop_refresh_parstype( accept=accept, SD=SD, MHprop=MHprop, SDchange=SDchange )
        MHprop$accept[[ var.vv ]] <- 0*MHprop$accept[[var.vv]]
        MHprop$refresh_count[[var.vv]] <- 0
    }
    return(MHprop)
}
