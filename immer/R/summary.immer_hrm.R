## File Name: summary.immer_hrm.R
## File Version: 0.35


#*******************************************************
# Summary for immer object
summary.immer_hrm <- function( object, digits=3, file=NULL, ... )
{
    #-- open sink
    immer_osink( file=file )

    #-- package description
    cat("-----------------------------------------------------------------\n")
    immer_summary_print_package_rsession(pack="immer")

    #-- analysis description
    cat( object$description, "\n\n")

    #-- call
    immer_summary_print_call(CALL=object$CALL)

    #-- computation time
    immer_summary_print_computation_time( object=object )

    cat( "Number of iterations ","=", object$iter, "\n" )
    cat( "Number of burnin iterations ","=", object$burnin, "\n" )
    cat( "Number of saved iterations ","=", object$N.save, "\n\n" )

    cat("-----------------------------------------------------------------\n")

    cat( "Deviance ","="," ", round( object$ic$dev, 2 ), " | " )
    cat( "Log Likelihood ","="," ", round( -object$ic$dev/2, 2 ), "\n" )
    cat( "Number of person-rater-interactions ","="," ", object$ic$ND, "\n" )
    cat( "Number of persons   ","="," ",  object$ic$N, "\n" )
    cat( "Number of items   ","="," ",  object$ic$I, "\n" )
    cat( "Number of raters   ","="," ",  object$ic$R, "\n\n" )

    cat( "Number of estimated parameters=", object$ic$np, "\n" )
    cat( "                      # mu=", object$ic$Npars["mu"], "\n" )
    cat( "                      # sigma=", object$ic$Npars["sigma"], "\n" )
    cat( "                      # a =", object$ic$Npars["a"], "\n" )
    cat( "                      # b =", object$ic$Npars["b"], "\n" )
    cat( "                      # phi=", object$ic$Npars["phi"], "\n" )
    cat( "                      # psi=", object$ic$Npars["psi"], "\n\n" )

    #-- print information criteria
    immer_summary_print_ic(object=object)

    cat("-----------------------------------------------------------------\n")
    cat( "Trait Distribution\n" )
    cat( "Mean=", round( object$mu, digits), " SD=", round( object$sigma, digits), "\n\n")
    cat( "EAP Reliability=")
    cat(round( object$EAP.rel, digits ) )
    cat( "\n")

    cat("-----------------------------------------------------------------\n")
    cat("Item Parameters \n")
    immer_summary_print_objects(obji=object$item, from=2, digits=digits, rownames_null=TRUE)

    cat("-----------------------------------------------------------------\n")
    cat("Rater Parameters \n")
    immer_summary_print_objects(obji=object$rater_pars, from=4, digits=digits, rownames_null=TRUE)

    cat("-----------------------------------------------------------------\n")
    cat("MCMC Diagnostics \n")
    immer_summary_print_mcmc_output(object=object, digits=digits)

    #-- close sink
    immer_csink( file=file )
}
#*******************************************************
