## File Name: summary.mnlfa.R
## File Version: 0.17


summary.mnlfa <- function( object, file=NULL, ... )
{
    CDM::osink( file=file, suffix=paste0( "__SUMMARY.Rout") )

    cat("-----------------------------------------------------------------------------\n")

    #-- print package
    # mnlfa_print_summary_package(pack="mnlfa")
    cat("\n")

    #-- print call
    CDM::cdm_print_summary_call(object=object)

    #-- print computation time
    CDM::cdm_print_summary_computation_time(object=object)

    cat("Moderated nonlinear factor analysis\n")

    CDM::cat_paste("\n-----------------------------------------------------------------------------\n")
    CDM::cat_paste("Number of iterations", xx(), object$iter, "\n" )
    if ( ! object$converged ){ cat("Maximum number of iterations was reached.\n") }

    CDM::cat_paste( "\nDeviance", xx(), round( object$deviance, 2 ), " | " )
    CDM::cat_paste( "Log Likelihood", xx(), round( -object$deviance/2, 2 ), "\n" )
    CDM::cat_paste( "Penalty", xx(), round( object$regular_penalty, 2 ), "\n" )

    CDM::cat_paste( "Number of persons", xx(), object$ic$n, "\n" )
    CDM::cat_paste( "Number of estimated parameters", xx(), object$ic$np, "\n" )
    CDM::cat_paste( "Number of regularized parameters", xx(), object$ic$numb_reg_pars, "\n\n" )

    #-- information criteria
    CDM::cdm_print_summary_information_criteria(object=object)

    cat("-----------------------------------------------------------------------------\n")
    cat("Item Parameters \n")
    obji <- object$item
    CDM::cdm_print_summary_data_frame(obji, digits=3)

    cat("-----------------------------------------------------------------------------\n")
    cat("Trait Distribution Parameters \n\n")
    obji <- object$parm_trait
    cat("Mean parameters\n")
    print( round(obji$mu, digits=4))
    cat("\nLog standard deviation parameters\n")
    print( round(obji$sigma, digits=4))

    CDM::csink( file=file )
}
#*******************************************************
