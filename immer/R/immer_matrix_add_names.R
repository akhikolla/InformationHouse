## File Name: immer_matrix_add_names.R
## File Version: 0.01

immer_matrix_add_names <- function(x, row=NULL, col=NULL)
{
    if ( ! is.null(row) ){
        rownames(x) <- row
    }
    if ( ! is.null(col) ){
        colnames(x) <- col
    }
    return(x)
}
