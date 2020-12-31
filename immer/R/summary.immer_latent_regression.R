## File Name: summary.immer_latent_regression.R
## File Version: 0.09


#*******************************************************
# Summary for immer object
summary.immer_latent_regression <- function( object, digits=3, file=NULL, ... )
{
    #-- open sink
    immer_osink( file=file )

    #-- package description
    cat("-----------------------------------------------------------------\n")
    immer_summary_print_package_rsession(pack="immer")

    cat( object$description, "\n\n")

    #-- call
    immer_summary_print_call(CALL=object$CALL)

    #-- computation time
    immer_summary_print_computation_time(object=object )

    cat( "Number of iterations=", object$iter, "\n" )
    cat("\n")

    cat( "Deviance=", round( object$ic$dev, 2 ), " | " )
    cat( "Log Likelihood=", round( -object$ic$dev/2, 2 ), "\n" )
    cat( "Number of persons=", object$ic$n, "\n" )
    cat( "Number of groups=", object$ic$G, "\n" )

    cat( "Number of estimated parameters=", object$ic$np, "\n" )
    cat("\n")

    #-- print information criteria
    immer_summary_print_ic(object=object)

    cat("-----------------------------------------------------------------\n")
    cat( "Regression coefficients\n" )
    immer_summary_print_objects(obji=object$beta_stat, from=2, digits=digits, rownames_null=TRUE)

    cat("-----------------------------------------------------------------\n")
    cat( "Standard deviations\n" )
    immer_summary_print_objects(obji=object$gamma_stat, from=2, digits=digits, rownames_null=TRUE)

    #-- close sink
    immer_csink( file=file )
}
#*******************************************************
