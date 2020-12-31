fun_covij <- function(x, b, length_b, h, arglist) {

  ej <- x[1]
  ek <- x[2]

  b_ej <- b + h*canon(ej, len = length_b)
  b_ek <- b + h*canon(ek, len = length_b)
  b_ej_ek <- b + h*canon(ek, len = length_b) + h*canon(ej, len = length_b)

  pen_lik_b <- log_penlikelihood(b, arglist)

  if (ej == ek) {
    pen_lik_b_ek <- pen_lik_b_ej <- log_penlikelihood(b_ej, arglist)
  } else {
    pen_lik_b_ej <- log_penlikelihood(b_ej, arglist)
    pen_lik_b_ek <- log_penlikelihood(b_ek, arglist)
  }

  pen_lik_b_ej_ek <- log_penlikelihood(b_ej_ek, arglist)
  value <- - (pen_lik_b_ej_ek - pen_lik_b_ej - pen_lik_b_ek + pen_lik_b)/h^2

  return(c(ej, ek, value))

}
