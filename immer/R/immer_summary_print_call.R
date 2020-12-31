## File Name: immer_summary_print_call.R
## File Version: 0.05

immer_summary_print_call <- function(CALL)
{
    cat("Call:\n", paste(deparse(CALL), sep="\n", collapse="\n"),
                "\n\n", sep="")
}
