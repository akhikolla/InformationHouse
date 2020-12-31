## File Name: summary.immer_ccml.R
## File Version: 0.08


#*******************************************************
# Summary for immer object
summary.immer_ccml <- function( object, digits=3, file=NULL, ... )
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

    cat( "Number of iterations ","=", object$nlminb_result$iterations, "\n" )
    cat( "Convergence code ","=", object$nlminb_result$convergence, "\n" )
    cat("\n")

    #cat("Estimation method=", object$est_method, "\n")
    # cat(" Epilson parameter=", object$eps, "\n")

    cat( "Objective function value ","=", round( object$objective, 2 ), "\n" )
    cat( "Number of person-rater-interactions ","=", object$ic$ND, "\n" )
    cat( "Number of persons   ","=", object$ic$N, "\n" )
    cat( "Number of items   ","=", object$ic$I, "\n" )
    cat( "Number of raters   ","=", object$ic$R, "\n\n" )

    cat( "Number of estimated item parameters ","=", object$ic$np, "\n" )
    cat("\n")

    #-- print information criteria
    immer_ccml_summary_print_ic(object=object)

    cat("-----------------------------------------------------------------\n")
    cat("Item Parameters \n")
    immer_summary_print_objects(obji=object$item, from=2, digits=digits, rownames_null=TRUE)

    cat("-----------------------------------------------------------------\n")
    cat("Basis Item Parameters \n")
    immer_summary_print_objects(obji=object$xsi_out, from=1, digits=digits, rownames_null=FALSE)

    #-- close sink
    immer_csink( file=file )
}
#*******************************************************
