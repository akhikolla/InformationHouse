HWPosterior <-
function(X, verbose = TRUE, prior.af = c(0.5, 0.5), 
                            prior.gf = c(0.333, 0.333, 0.333), x.linked = FALSE, 
                            precision = 0.05) 
{
  if(x.linked) { # X-chromosomal
    lab <- names(X)
    if (is.null(lab)) {
      cat("No genotype labels given, default order c(A,B,AA,AB,BB) assumed.\n")
      X <- genlabels(X)
      lab <- names(X)
    }
    if (length(X) != 5) {
      cat("Improper number of genotype counts.\n")
      stop()
    }
    if (!all(lab %in% c("A", "B", "AA", "AB", "BB"))) 
      stop("Unknown genotypes occurred. Supply counts as a named vector c(A,AA,AB,B,BB)")
    n <- sum(X)
    nfAA <- X[lab == "AA"]
    nfAB <- X[lab == "AB"]
    nfBB <- X[lab == "BB"]
    nmA <- X[lab == "A"]
    nmB <- X[lab == "B"]
    x.m <- c(nmA, nmB)
    x.f <- c(nfAA, nfAB, nfBB)
    p.h0 <- H0(x.f, x.m, alpha = prior.af)
    p.h1 <- H1(x.f, x.m, alpha = prior.gf)
    p.h2 <- H2(x.f, x.m, alpha.f = prior.af, alpha.m = prior.af)
    p.h3 <- H3(x.f, x.m, alpha.f = prior.gf, alpha.m = prior.af)
    posterior.prob <- c(p.h0, p.h1, p.h2, p.h3)
    posterior.prob <- posterior.prob/sum(posterior.prob)
    lBF <- log10(3 * posterior.prob/(1 - posterior.prob))
    r.posterior.prob <- round(posterior.prob, 4)
    Res <- cbind(posterior.prob, lBF)
    colnames(Res) <- c("Posterior_Prob", "log10(Bayes Factor)")
    rownames(Res) <- c("M0 (HWE):", "M1 (f!=0):", "M2 (d!=1):", 
                       "M3 (f!=0 & d!=1:)")
    if (verbose) {
      cat("Bayesian test for Hardy-Weinberg equilibrium of X-chromosomal variants.\n\n")
      print(round(Res, digits = 4))
    }
    out <- Res 
  } else { # autosomal
    lab <- names(X)
    if (is.null(lab)) {
      cat("No genotype labels given, default order c(AA,AB,BB,AA,AB,BB) assumed, for males and females respectively\n")
      X <- c("AA","AB","BB","AA","AB","BB")
      lab <- names(X)
    }
    if (length(X) != 6) {
      cat("Improper number of genotype counts.\n")
      stop()
    }
    posterior.prob <- c(M11p(X,alpha=prior.af),
              M12p(X,alpha=prior.gf,r=precision),
              M13p(X,alpha=prior.gf,r=precision),
              M14p(X,alpha=prior.gf),
              M15p(X,alpha.f=prior.gf,r=precision), # alpha.m default
              M21p(X,alpha.f=prior.af,alpha.m=prior.af),
              M22p(X,alpha.f=prior.gf,alpha.m=prior.af),
              M23p(X,alpha.f=prior.af,alpha.m=prior.gf),
              M24p(X,alpha.f=prior.gf,alpha.m=prior.af,r=precision),
              M25p(X,alpha.f=prior.gf,alpha.m=prior.gf))
    posterior.prob <- posterior.prob/sum(posterior.prob)
    names(posterior.prob) <- c("M_11", "M_12", "M_13", "M_14", "M_15", "M_21", 
                    "M_22", "M_23", "M_24", "M_25")
    if(verbose) {
      print(round(posterior.prob,digits=4))
      indexbest <- which.max(posterior.prob)
      cat("Best fitting", names(posterior.prob[indexbest]), posterior.prob[indexbest], 
          "\n")
    }
    out <- posterior.prob
  }
}
