## File Name: summary.immer_agree2.R
## File Version: 0.11


#*******************************************************
# Summary for immer object
summary.immer_agree2 <- function( object, digits=3,... )
{
    cat("-----------------------------------------------------------------\n")
    d1 <- utils::packageDescription("immer")
    cat( paste( d1$Package, " ", d1$Version, " (", d1$Date, ")", sep=""), "\n\n" )

    # cat( object$description, "\n\n")

    cat("Call:\n", paste(deparse(object$CALL), sep="\n", collapse="\n"),
                "\n\n", sep="")


    cat("-----------------------------------------------------------------\n")

    cat( "Number of persons (sum of weights) ","=", object$nobs, "\n" )
    cat( "Number of categories ","=", object$ncat, "\n" )

    cat("-----------------------------------------------------------------\n")
    cat("Contingency Table \n\n")

    obji <- round( object$agree_table    , digits=digits )
    print( obji )

    cat("-----------------------------------------------------------------\n")
    cat("Marginal Frequencies \n\n")
    obji <- round( object$marg    , digits=digits )
    print( obji )

    cat("-----------------------------------------------------------------\n")
    cat("Raw Agreement Statistics \n")
    obji <- object$agree_raw
    obji <- round( obji, digits=digits )
    print(obji)

    cat("-----------------------------------------------------------------\n")
    cat("Agreement Statistics \n\n")

    cat( "Raw agreement=", round(object$Pa, digits=digits ), "\n" )
    cat( "Scott's Pi=",
            round(object$agree_stats["pi"], digits=digits ), "\n" )
    cat( "Cohen's Kappa ","=",
            round(object$agree_stats["kappa"], digits=digits ), "\n" )
    cat( "Aicken's Alpha ","=",
            round(object$agree_stats["Aicken"], digits=digits ), "\n" )

    cat( "Gwet's AC1 ","=",
            round(object$agree_stats["AC1"], digits=digits ), "\n" )

    cat("-----------------------------------------------------------------\n")
    cat("Conditional Probabilities for Hard-to-classify Persons (Aicken, 1990)\n\n")
    obji <- round( object$PH    , digits=digits )
    print( obji )
}
#*******************************************************
