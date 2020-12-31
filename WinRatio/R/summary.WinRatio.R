summary.WinRatio <- function(object, ..., digits = 2){
  
  cat("Number of subjects: N =", object$n, "\n")
  cat(paste0("Number of subjects in group '", object$group1, "': ", object$n1, "\n")) 
  cat(paste0("Number of subjects in group '", object$group0, "': ", object$n0, "\n"))
  cat(paste0("Number of paired comparison: ", object$n1 * object$n0, "\n"))
  cat(paste0("Active group = group '", object$group1, "'\n"))
  
  cat(" ","\n")
  mat <- matrix(data = NA, nrow = length(object$call$outcomes)+1, ncol = 2)
  rownames(mat) <- c(paste0("Outcome ", 1:length(object$call$outcomes)), "Total")
  colnames(mat) <- c("Numbers of 'winners'", "Numbers of 'losers'")
  mat[length(object$call$outcomes)+1,] <- c(object$total.wins, object$total.loss)
  outcomes <- NULL
  for (k in 1:length(object$call$outcomes)){
    if (object$call$outcomes[[k]][2] == "s"){
      outcomes <- c(outcomes, paste0("Outcome ", k, " = ", object$call$outcomes[[k]][1], 
                                     ", ",object$call$outcomes[[k]][3], " (survival event)"))
    } else if (object$call$outcomes[[k]][2] == "r"){
      outcomes <- c(outcomes,paste0("Outcome ", k, " = ", 
                                    paste0(utils::head(object$call$outcomes[[k]][[1]], 1), " ... ", utils::tail(object$call$outcomes[[k]][[1]], 1)), ", ",
                                    paste0(utils::head(object$call$outcomes[[k]][[3]], 1), " ... ", utils::tail(object$call$outcomes[[k]][[3]], 1)),
                                    " (repeated events)"))
    } else if (object$call$outcomes[[k]][2] == "c"){    
      outcomes <- c(outcomes, paste0("Outcome ", k, " = ", object$call$outcomes[[k]][1],
                                     " (continuous event, direction = '", object$call$outcomes[[k]][3], "')"))
    }
    mat[k,] <- c(object$wins[k], object$loss[k])
  }
  paste(object$call$outcomes[[k]][1])
  cat(paste0(outcomes, collapse = "\n"))
  
  cat("\n\n")
  
  print(mat)
  
  signif <- dplyr::case_when(object$p.value < 0.001 ~ "***",
                      object$p.value < 0.01 ~ "**",
                      object$p.value < 0.05 ~ "*",
                      T ~ "")
  pval <- dplyr::if_else(object$p.value < 0.0001, 
                  paste0("p-value < 0.0001", signif), 
                  paste0("p-value = ", sprintf(paste0("%.", 4, "f"), object$p.value), signif))
  
  cat("\n")
  cat(paste0("Total number of ties: ", object$total.ties))
  
  cat("\n\n")
  cat(paste0("Win Ratio (CI 95%) = ",sprintf(paste0("%.", digits, "f"), object$wr), " (", 
             sprintf(paste0("%.", digits, "f"), object$wr.lower), " - ",
             sprintf(paste0("%.", digits, "f"), object$wr.upper), "), ", pval,
             "\n\nSignif. codes: '***' p < 0.001, '**' p < 0.01, '*' p < 0.05."))
  

}
