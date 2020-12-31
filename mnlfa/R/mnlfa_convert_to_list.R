## File Name: mnlfa_convert_to_list.R
## File Version: 0.05


mnlfa_convert_to_list <- function(x, I, names_list=NULL, create_list=FALSE)
{
    if ( ( ! is.list(x) ) | create_list ){
        x0 <- x
        x <- list()
        for (ii in 1:I){
            x[[ii]] <- x0
        }
    }
    if ( ! is.null(names_list) ){
        names(x) <- names_list
    }
    #-- output
    return(x)
}
