dlevene.bi <-
function (x) # Levene's density for bi-allelic autosomal
{
  nAA <- x["AA"]
  nAB <- x["AB"]
  nBB <- x["BB"]
  n <- sum(x)
  nA <- 2 * nAA + nAB
  nB <- 2 * n - nA
    numer <- lfactorial(n) + lfactorial(nA) + lfactorial(2 * 
                                                         n - nA) + nAB * log(2)
  denom <- lfactorial(2 * n) + lfactorial(nAB) + lfactorial(0.5 * 
                                                              (nA - nAB)) + lfactorial(0.5 * (nB - nAB))
  prob <- exp(numer - denom)
  return(prob)
}
