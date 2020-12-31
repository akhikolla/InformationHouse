#HELPER: MULTIVARIATE NORMAL DISTRIBUTION
#SIMPLIFIED VERSION OF mvrnorm ADAPTED FROM MASS 
### LAST UPDATE: NA 
#' @importFrom stats rnorm
simp.mvrnorm <- function(n=1, mu, Sigma){
  p <- length(mu)
  eS <- eigen(Sigma, symmetric = TRUE)
  ev <- eS$values
  X <- matrix(rnorm(p * n), n)
  X <- drop(mu)+eS$vectors%*%diag(sqrt(pmax(ev,0)),p)%*%t(X)
  if (n == 1)
    drop(X)
  else t(X)
}

#PROGRESS BAR WITH MESSAGE
### LAST UPDATE: 11/08/2020; Le Bao
#' @importFrom utils flush.console
progBar <- function(iter, total, progress.bar="message",nugget, decay, partial.sill, 
                    interac = interactive()) {
  if (!interac) {
    if (total < 10 | iter %% round(total/10) == 0) {
      if (progress.bar=="message"){
        pb.noniter.ar(iter=iter,total=total,nugget=nugget,
                      decay=decay,partial.sill=partial.sill)
      } else if (progress.bar==TRUE) {pb.noniter(iter=iter,total=total)}
    }
  } else {
    if (iter < 5 | iter %% round(total/100,0) == 0 | round(total/100, 0) == 0 |
        iter > total-5) {
      if (progress.bar=="message"){
        pb.ar(iter=iter,total=total,nugget=nugget,
              decay=decay,partial.sill=partial.sill)
      } else if (progress.bar==TRUE){
        pb(iter=iter,total=total)
      }  
    }
  }
}

  
pb.ar <-function(iter,total,nugget,decay,partial.sill){
  p <- floor(iter/total*100)
  if (p < 10){
    mess <- sprintf("Nugget=%.3f; Decay=%.3f; Partial sill=%.3f............................%i%%", 
                    nugget, decay, partial.sill, p)
  }
  if (p >= 10){
    mess <- sprintf("Nugget=%.3f; Decay=%.3f; Partial sill=%.3f...........................%i%%", 
                    nugget, decay, partial.sill, p)
  }
  if (p == 100){
    mess <- sprintf("Nugget=%.3f; Decay=%.3f; Partial sill=%.3f..........................%i%%", 
                    nugget, decay, partial.sill, p)
  }
  nmess <- nchar(mess)
  if (iter > 1) cat(.makeMessage(paste(rep("\b", nmess), collapse = "")))
  cat(.makeMessage(mess))
  if (iter == total) { cat("\n") }
  flush.console()
}

#SIMPLE PROGRESS BAR
### LAST UPDATE: 11/08/2020; Le Bao

#' @importFrom utils setTxtProgressBar txtProgressBar
pb <- function(iter,total){
  pb.set <- txtProgressBar(min=0, max = total, style = 3)
  # update progress bar
  setTxtProgressBar(pb.set, iter)
}

pb.noniter.ar <- function(iter,total,nugget,decay,partial.sill){
  p <- round(iter/total*100)
  if (p < 100){
  cat(sprintf("Nugget=%.3f; Decay=%.3f; Partial sill=%.3f...........................%i%%", 
          nugget, decay, partial.sill, p), "\n")
  }
  if (p == 100){
    cat(sprintf("Nugget=%.3f; Decay=%.3f; Partial sill=%.3f..........................%i%%", 
                    nugget, decay, partial.sill, p), "\n")
  }
}
  
pb.noniter <- function(iter,total){
  p <- round(iter/total*100)
  if (p < 100){
    cat(sprintf("%i%%...", p))
  }
  if (p == 100){
    cat(sprintf("%i%%", p), "\n")
  }
}

# HELPER: ACCEPTANCE RATES CHECKER
### LAST UPDATE: 11/08/2020; Le Bao
ar.check <- function(ar.rate){
  cat("\n")
  if(ar.rate$beta.rate<.2) 
    message("*Coefficient acceptance rate is low, which can indicate slow mixing. You may need to run many iterations or consider adjusting 'b.tune'.")
  if(ar.rate$beta.rate>.7) 
    message("*Coefficient acceptance rate is high. If coefficients are nonconvergent, consider adjusting 'b.tune'.")
  if(ar.rate$tau2.rate<.2) 
    message("*Nugget acceptance rate is low, which can indicate slow mixing. You may need to run many iterations or consider adjusting 'nugget.tune'.")
  if(ar.rate$tau2.rate>.7) 
    message("*Nugget acceptance rate is high. If tau2 is nonconvergent, consider adjusting 'nugget.tune'.")
  if(ar.rate$phi.rate<.2) 
    message("*Decay term acceptance rate is low, which can indicate slow mixing. You may need to run many iterations or consider adjusting 'range.tol'.")
  if(ar.rate$phi.rate>.7) 
    message("*Decay term acceptance rate is high. If phi is nonconvergent, consider adjusting 'range.tol'.")
  if(ar.rate$sigma2.rate<.2) 
    message("*Partial sill acceptance rate is low, which can indicate slow mixing. You may need to run many iterations or consider adjusting 'psill.tune'.")
  if(ar.rate$sigma2.rate>.7) 
    message("*Partial sill acceptance rate is high. If sigma2 is nonconvergent, consider adjusting 'psill.tune'.")
}