ddir <-
function (x, alpha, log = T) {  # log = F
    x <- as.vector(x)
    alpha <- as.vector(alpha)
    if (!identical(length(x), length(alpha))) 
        stop("x and alpha differ in length.")
    dens <- lgamma(sum(alpha)) +  sum((alpha-1)*log(x)) -
		sum(lgamma(alpha)) 
    if (log == FALSE) 
        dens <- exp(dens)
    return(dens)
}
