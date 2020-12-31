## File Name: immer_summary_print_ic.R
## File Version: 0.07

immer_summary_print_ic <- function(object)
{
    cat( "AIC=", round( object$ic$AIC, 2 ), " | penalty=",
            round( object$ic$AIC - object$ic$dev,2 ),
            "   | AIC=-2*LL + 2*p  \n" )
    cat( "AIC3=", round( object$ic$AIC3, 2 ), " | penalty=",
                round( object$ic$AIC3 - object$ic$dev,2 ),
                "   | AIC3=-2*LL + 3*p  \n" )
    cat( "AICc=", round( object$ic$AICc, 2 )," | penalty=",
            round( object$ic$AICc - object$ic$dev,2 ) )
        cat("    | AICc=-2*LL + 2*p + 2*p*(p+1)/(n-p-1)  (bias corrected AIC)\n" )
    cat( "BIC=", round( object$ic$BIC, 2 ), " | penalty=",
                round( object$ic$BIC - object$ic$dev,2 ),
            "   | BIC=-2*LL + log(n)*p  \n" )
    cat( "aBIC=", round( object$ic$aBIC, 2 ), " | penalty=",
                round( object$ic$aBIC - object$ic$dev,2 ),
            "   | aBIC=-2*LL + log((n-2)/24)*p  (adjusted BIC) \n" )
    cat( "CAIC=", round( object$ic$CAIC, 2 )," | penalty=",
                    round( object$ic$CAIC - object$ic$dev,2 ) )
        cat("   | CAIC=-2*LL + [log(n)+1]*p  (consistent AIC)\n\n" )
}
