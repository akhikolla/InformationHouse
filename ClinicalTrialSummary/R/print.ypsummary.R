print.ypsummary <- function(x, ...) {
    digit <- paste("%.", max(3, getOption("digits") - 3), "f", sep = "")
    tau <- x$tau
    #--- The average hazard ratio ---#
    x_ahrypt <- x$values$ahrypt
    print_ahrypt <- lapply(X = x_ahrypt, FUN = fun_digit, digit = digit)
    print_ahrypt$pvalue <- fun_less(print_ahrypt$pvalue)

    #--- The weighted average hazard ratio ---#
    x_ahrf <- x$values$ahrf
    print_ahrf <- lapply(X = x_ahrf, FUN = fun_digit, digit = digit)
    print_ahrf$pvalue <- fun_less(print_ahrf$pvalue)

    #--- The restricted superiority probability ratio ---#
    x_rspypt <- x$values$rspypt
    print_rspypt <- lapply(X = x_rspypt, FUN = fun_digit, digit = digit)
    print_rspypt$pvalue <- fun_less(print_rspypt$pvalue)

    #--- The restricted mean survival difference ---#
    x_mdypt <- x$values$mdypt
    print_mdypt <- lapply(X = x_mdypt, FUN = fun_digit, digit = digit)
    print_mdypt$pvalue <- fun_less(print_mdypt$pvalue)

    #--- The ratio of restricted mean times lost ---#
    x_mrypt <- x$values$mrypt
    print_mrypt <- lapply(X = x_mrypt, FUN = fun_digit, digit = digit)
    print_mrypt$pvalue <- fun_less(print_mrypt$pvalue)

    cat("\n=====================================================\n")
    cat(paste("    Summary Measures (tau = ", tau, ")", sep = ""))
    cat("\n=====================================================\n")
    cat("    AHR      Estimate (CI):", print_ahrypt$estimate, paste("(", print_ahrypt$lower,
                                                                    ", ", print_ahrypt$upper, ")", sep = ""), "\n")
    cat("             z-value (p-value):", print_ahrypt$z, paste("(", print_ahrypt$pvalue,
                                                                 ")", sep = ""))
    cat("\n -----------------------------------------------------\n")
    cat("    WAHR     Estimate (CI):", print_ahrf$estimate, paste("(", print_ahrf$lower,
                                                                  ", ", print_ahrf$upper, ")", sep = ""), "\n")
    cat("             z-value (p-value):", print_ahrf$z, paste("(", print_ahrf$pvalue,
                                                               ")", sep = ""))
    cat("\n -----------------------------------------------------\n")
    cat("    RSPR     Estimate (CI):", print_rspypt$estimate, paste("(", print_rspypt$lower,
                                                                    ", ", print_rspypt$upper, ")", sep = ""), "\n")
    cat("             z-value (p-value):", print_rspypt$z, paste("(", print_rspypt$pvalue,
                                                                 ")", sep = ""))
    cat("\n -----------------------------------------------------\n")
    cat("    RMSD     Estimate (CI):", print_mdypt$estimate, paste("(", print_mdypt$lower,
                                                                   ", ", print_mdypt$upper, ")", sep = ""), "\n")
    cat("             z-value (p-value):", print_mdypt$z, paste("(", print_mdypt$pvalue,
                                                                ")", sep = ""))
    cat("\n -----------------------------------------------------\n")
    cat("    RRMTL    Estimate (CI):", print_mrypt$estimate, paste("(", print_mrypt$lower,
                                                                   ", ", print_mrypt$upper, ")", sep = ""), "\n")
    cat("             z-value (p-value):", print_mrypt$z, paste("(",print_mrypt$pvalue,
                                                                ")", sep = ""))
    cat("\n=====================================================\n","\n")
}

