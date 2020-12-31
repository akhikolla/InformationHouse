HWAIC <-
function(x, y, tracing = 0, tol=0.000001) 
{
  z <- x + y  # x males y females
  pam <- af(x)
  paf <- af(y)
  pa <- af(z)
  fm <- HWf(x)
  ff <- HWf(y)
  f <- HWf(z)
  
  l11 <- loglik.M11(pa, z)
  
  l12 <- loglik.M12(x, y, tol=tol)
  
  l13 <- loglik.M12(y, x, tol=tol) # swap males and females to fit the model
  
  l14 <- loglik.M14(pa, z, f)
  
  l15 <- loglik.M15(pa, f, x, y, tracing = tracing)
  
  l21 <- loglik.M21(x, y, pam, paf)
  
  l22 <- loglik.M22(pam, paf, x, y, ff)    
  
  l23 <- loglik.M23(pam, paf, x, y, fm) 
  
  l24 <- loglik.M24(x, y, tracing = tracing)
  
  l25 <- loglik.M25(x, y, pam, paf, fm, ff)
  
  # loglikelihoods
  
  ll <- c(l11[1], l12[1], l13[1], l14[1], l15[1], 
          l21[1], l22[1], l23[1], l24[1], l25[1])
  
  # number of parameters of each scenario
  
  np <- c(1, 2, 2, 2, 3, 
          2, 3, 3, 3, 4)
  
  aic <- 2 * np - 2 * ll
  
  names(aic) <- c("M_11", "M_12", "M_13", "M_14", "M_15", 
                  "M_21", "M_22", "M_23", "M_24", "M_25")
  
  indexbest <- which.min(aic)
  cat("Best fitting",names(aic[indexbest]),aic[indexbest],"\n")
  return(aic)
}
