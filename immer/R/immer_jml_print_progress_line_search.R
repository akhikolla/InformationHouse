## File Name: immer_jml_print_progress_line_search.R
## File Version: 0.06

immer_jml_print_progress_line_search <- function( verbose, deviance, digits_deviance=6)
{
    if (verbose){
        cat("  - Line search | Deviance=", round( deviance, 6 ), "\n")
        utils::flush.console()
    }
}
