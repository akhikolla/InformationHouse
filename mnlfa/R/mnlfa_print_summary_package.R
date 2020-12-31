## File Name: mnlfa_print_summary_package.R
## File Version: 0.05


mnlfa_print_summary_package <- function (pack="mnlfa")
{
    d1 <- utils::packageDescription(pack)
    cat(paste(d1$Package, " ", d1$Version, " (", d1$Date, ")", sep=""), "\n")
}
