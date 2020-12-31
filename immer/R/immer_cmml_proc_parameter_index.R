## File Name: immer_cmml_proc_parameter_index.R
## File Version: 0.04

immer_cmml_proc_parameter_index <- function( W_par, par_index )
{
    NL <- max(par_index)
    if (NL==0){
        NL <- 1
    }
    if ( is.null(W_par)){
        W_par <- matrix( 0, nrow=NL, ncol=1)
    }
    return(W_par)
}
