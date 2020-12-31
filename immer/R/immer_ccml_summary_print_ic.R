## File Name: immer_ccml_summary_print_ic.R
## File Version: 0.07

immer_ccml_summary_print_ic <- function(object)
{
    cat( "CLAIC=", round( object$ic$CLAIC, 2 ), " | penalty=",
            round( object$ic$CLAIC - object$ic$dev,2 ),
            "   | CLAIC=-2*CL + 2*tr(J*H^-1)  \n" )
    cat( "CLBIC=", round( object$ic$CLBIC, 2 ), " | penalty=",
                round( object$ic$CLBIC - object$ic$dev,2 ),
            "   | CLBIC=-2*CL + log(n)*tr(J*H^-1)  \n" )
    cat("\n")
}
