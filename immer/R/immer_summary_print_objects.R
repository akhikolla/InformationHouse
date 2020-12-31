## File Name: immer_summary_print_objects.R
## File Version: 0.02


immer_summary_print_objects <- function(obji, from=NULL, to=NULL,
        digits=3, rownames_null=TRUE)
{
    sirt::sirt_summary_print_objects(obji=obji, from=from, to=to, digits=digits,
            rownames_null=rownames_null)
}
