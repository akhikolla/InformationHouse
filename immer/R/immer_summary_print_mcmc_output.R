## File Name: immer_summary_print_mcmc_output.R
## File Version: 0.04

immer_summary_print_mcmc_output <- function(object, digits)
{
    summ <- object$summary.mcmcobj
    obji <- summ[, c("parameter", "Mean", "MAP", "SD", "Q2.5", "Q97.5" ) ]
    # obji$MCMC_SE <- summ$Time.series.SE
    vars <- c("Rhat","effSize" )
    for (vv in vars ){
        obji[,vv] <- summ[,vv]
    }
    obji$effSize <- round(obji$effSize, 0 )
    sirt::sirt_summary_print_objects(obji=obji, from=2, digits=digits, rownames_null=TRUE)
}
