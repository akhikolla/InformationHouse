## File Name: summary.immer_jml.R
## File Version: 0.20


#*******************************************************
# Summary for immer object
summary.immer_jml <- function( object, digits=3, file=NULL, ... )
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
    cat( "Iteration with minimal deviance=", object$iter_opt, "\n" )
    cat("\n")

    cat("Estimation method=", object$est_method, "\n")
    cat(" Epilson parameter=", object$eps, "\n")
    if ( object$est_method=="jml_bc" ){
        cat(" Bias correction factor=", round(object$bc_adj_fac, digits), "\n")
    }
    cat("\n")

    cat( "Deviance=", round( object$ic$dev, 2 ), " | " )
    cat( "Log Likelihood=", round( -object$ic$dev/2, 2 ), "\n" )
    cat( "Number of person-rater-interactions=", object$ic$ND, "\n" )
    cat( "Number of persons=", object$ic$n, "\n" )
    cat( "Number of items=", object$ic$I, "\n" )
    cat( "Number of raters=", object$ic$R, "\n\n" )

    cat( "Number of estimated parameters=", object$ic$np, "\n" )
    cat( "                       # theta=", object$ic$ntheta, "\n" )
    cat( "                       # xsi=", object$ic$nxsi, "\n" )
    cat("\n")

    cat( "irtmodel=", object$irtmodel, "\n" )
    cat("\n")


    #-- print information criteria
    immer_summary_print_ic(object=object)

    cat("-----------------------------------------------------------------\n")
    cat( "Trait Distribution\n" )
    cat("  M=", round( object$person_desc$mean, 3), "\n",
        " SD manifest=", round( object$person_desc$sd_obs, 3), "\n",
        " SD latent=", round( object$person_desc$sd_lat, 3), "\n",
        " MLE reliability=", round( object$person_desc$mle_rel, 3), "\n",
        " Min=", round( object$person_desc$min, 3), "\n",
        " Max=", round( object$person_desc$max, 3), "\n"
        )

    cat("-----------------------------------------------------------------\n")
    cat("Item Parameters \n")
    immer_summary_print_objects(obji=object$item, from=2, digits=digits, rownames_null=TRUE)

    cat("-----------------------------------------------------------------\n")
    cat("Basis Item Parameters \n")
    immer_summary_print_objects(obji=object$xsi_dfr, from=2, digits=digits, rownames_null=TRUE)

    #-- close sink
    immer_csink( file=file )
}
#*******************************************************
