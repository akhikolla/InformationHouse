create_share <- function(fun_grad) {
  env <- new.env()
  env$xold <- NULL
  env$fun <- function(x) {
    if (is.null(env$xold) || any(env$xold != x)) {
      out <- fun_grad(x)
      env$xold <- x
      env$fun_val <- out[[1]]
      env$grad_val <- out[[2]]
    }
    env$fun_val
  }
  env$grad <- function(x = NULL) {
    if (is.null(env$xold) || any(env$xold != x)) {
      out <- fun_grad(x)
      env$xold <- x
      env$fun_val <- out[[1]]
      env$grad_val <- out[[2]]
    }
    env$grad_val
  }
  env
}

optim_fad <- function(par, fngr, ...) {
  env <- create_share(fngr)
  optim(par=par, fn=env$fun, gr=env$grad, ...)
}


