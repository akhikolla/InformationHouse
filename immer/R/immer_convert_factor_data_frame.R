## File Name: immer_convert_factor_data_frame.R
## File Version: 0.01

immer_convert_factor_data_frame <- function(x, variable )
{
    if ( variable %in% colnames(x) ){
        x[, variable ] <- as.factor( x[,variable ] )
    }
    #-- output
    return(x)
}
