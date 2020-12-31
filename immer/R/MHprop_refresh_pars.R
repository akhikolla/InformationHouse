## File Name: MHprop_refresh_pars.R
## File Version: 0.04


MHprop_refresh_pars <- function( acc, SD.pp, MHprop, SDchange )
{
    target <- mean( MHprop$accept_bounds )
    if (MHprop$refresh_formula){
        SD.pp <- ifelse( acc < MHprop$accept_bounds[1], SD.pp / ( 2 - acc / target ), SD.pp )
        SD.pp <- ifelse( acc > MHprop$accept_bounds[2], SD.pp * ( 2 - (1-acc)/(1-target) ), SD.pp )
    } else {
        SD.pp <- ifelse( acc < MHprop$accept_bounds[1], SD.pp - SDchange, SD.pp )
        SD.pp <- ifelse( acc > MHprop$accept_bounds[2], SD.pp + SDchange, SD.pp )
        SD.pp <- ifelse( SD.pp < SDchange, SDchange, SD.pp )
    }
    return(SD.pp)
}
