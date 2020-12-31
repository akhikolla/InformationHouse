#' @export
#' @method print emfrail_summary
#' @keywords internal
print.emfrail_summary <- function(x,
                                  digits = max(getOption('digits')-4, 3),
                                  signif.stars =
                                    getOption("show.signif.stars"),
                                  ...) {

  # cat("Summary of emfrail fit\n")
  cat("Call: \n")
  dput(attr(x, "call"))
  cat("\n")

  if(!identical(attr(x, "print_opts")$coef, FALSE)) {
    if(!is.null(x$coefmat)) {
      cat("Regression coefficients:\n")
      printCoefmat(x$coefmat,
                   digits = digits,
                   signif.stars = signif.stars, ...)
    }
  }

  if(!identical(attr(x, "print_opts")$dist, FALSE)) {
    cat("Estimated distribution:", x$est_dist$dist, "/ left truncation:", x$est_dist$left_truncation,"\n")
    if(x$est_dist$dist == "pvf") {
      cat("PVF m =", x$est_dist$pvfm," ")
      if(x$est_dist$pvfm == -0.5) cat("(Inverse Gaussian)")
      cat("\n")
    }
    cat("\n")
  }

  if(!identical(attr(x, "print_opts")$fit, FALSE)) {
    cat("Fit summary:\n")
    if(!is.null(x$ca_test))
      cat("Commenges-Andersen test for heterogeneity: p-val ", format(x$ca_test[3], digits = 3), "\n")
    cat("no-frailty Log-likelihood:", round(x$loglik[1], digits = 3), "\n")
    cat("Log-likelihood:", round(x$loglik[2], digits = 3), "\n")
    if(x$loglik[4] < 0.5)
      cat("LRT: 1/2 * pchisq(", format(x$loglik[3], digits = 3),"), p-val ",
        format(x$loglik[4], digits = 3), "\n", sep = "") else
          cat("LRT: 1/2 * pchisq(", format(x$loglik[3], digits = 3),"), p-val>",
              format(x$loglik[4], digits = 3), "\n", sep = "")

    cat("\n")
  }


  if(!identical(attr(x, "print_opts")$frailty, FALSE)) {

    cat("Frailty summary:\n")

    fsum <- data.frame()

    # gamma and pvf have this
    if(!is.null(x$fr_var)) {
      fsum <- rbind(fsum,  c(round(x$fr_var[1], digits = 3), round(x$fr_var[3], digits = 3), round(x$fr_var[4], digits = 3)))
      dimnames(fsum)[[1]] <- "Var[Z]"
    }

      # cat("frailty variance = ",
      #     round(x$fr_var[1], digits = 3),
      #     # " (",
      #     # round(x$fr_var[2], digits = 2),
      #     # ") / 95% CI: [",
      #     " / 95% CI: [",
      #     round(x$fr_var[3], digits = 3),
      #     ", ",
      #     round(x$fr_var[4], digits = 3),
      #     "]\n", sep = "")

    if(!identical(attr(x, "print_opts")$verbose_frailty, FALSE)) {

    if(!is.null(x$gamma_pars)) {

      tmp <- as.data.frame(with(x$gamma_pars, rbind(c(round(tau[[1]], digits = 3),  round(tau[[3]], digits = 3), round(tau[[4]], digits = 3)),
                                                    c(round(kappa[[1]], digits = 3),  round(kappa[[3]], digits = 3), round(kappa[[4]], digits = 3)),
                                                    c(round(e_log_z[[1]], digits = 3), round(e_log_z[[3]], digits = 3),  round(e_log_z[[4]], digits = 3)),
                                                    c(round(var_log_z[[1]], digits = 3), round(var_log_z[[3]], digits = 3), round(var_log_z[[4]], digits = 3))))
      )

      dimnames(tmp)[[2]] <- dimnames(fsum)[[2]]
      dimnames(tmp)[[1]] <- c("Kendall's tau", "Median concordance", "E[logZ]", "Var[logZ]")
      fsum <- rbind(fsum, tmp)



      }

      # browser()
    if(!is.null(x$stable_pars)) {


      tmp <- as.data.frame(with(x$stable_pars,
                                rbind(c(round(tau[[1]], digits = 3),  round(tau[[3]], digits = 3), round(tau[[4]], digits = 3)),
                                c(round(kappa[[1]], digits = 3),  round(kappa[[3]], digits = 3), round(kappa[[4]], digits = 3)),
                                c(round(e_log_z[[1]], digits = 3), round(e_log_z[[3]], digits = 3),  round(e_log_z[[4]], digits = 3)),
                                c(round(var_log_z[[1]], digits = 3), round(var_log_z[[3]], digits = 3), round(var_log_z[[4]], digits = 3)),
                                c(round(attenuation[[1]], digits = 3), round(attenuation[[3]], digits = 3),round(attenuation[[4]], digits = 3))

      )))


      dimnames(tmp)[[2]] <- c("A","B","C")
      dimnames(tmp)[[1]] <- c("Kendall's tau", "Median concordance", "E[logZ]", "Var[logZ]", "Attenuation")
      fsum <- rbind(fsum, tmp)

    }


    if(!is.null(x$pvf_pars)) {
        if(x$est_dist$pvfm < 0) {
          tmp <- as.data.frame(with(x$pvf_pars, rbind(c(round(tau[[1]], digits = 3),  round(tau[[3]], digits = 3), round(tau[[4]], digits = 3)),
                                                      c(round(kappa[[1]], digits = 3),  round(kappa[[3]], digits = 3), round(kappa[[4]], digits = 3)))

          ))

          dimnames(tmp)[[2]] <- dimnames(fsum)[[2]]
          dimnames(tmp)[[1]] <- c("Kendall's tau", "Median concordance")
          fsum <- rbind(fsum, tmp)
        }

          if(x$est_dist$pvfm > 0)
          cat("Estimated mass at 0:",
              x$pvf_pars[[1]],
              "\n")

      }

    tmp <- as.data.frame(t(c(round(x$theta[1], digits = 3),  round(x$theta[3], digits = 3), round(x$theta[4], digits = 3))))

    dimnames(tmp)[[2]] <- dimnames(fsum)[[2]]
      dimnames(tmp)[[1]] <- c("theta")
      fsum <- rbind(fsum, tmp)


      dimnames(fsum)[[2]] <- c("estimate", "lower 95%", "upper 95%")
      print(fsum)

    if(!is.null(x$cens_test)) {
      cat("\n")
      cat("Score test for dependent censoring: p-val", format(x$cens_test[2], digits = 3))
      cat("\n")
    }

    # if(!is.null(x$pvf_pars))
    #   cat("Estimated mass at 0:",
    #       x$pvf_pars[[1]],
    #       "\n")

    if(isTRUE(x$lik_ci)) cat("Confidence intervals based on the likelihood function") else
      cat("Confidence intervals based on the delta method")

    cat("\n")
    }

  }



}
