projcsp <- function(xstar, yhat) {
    theta = solve(t(xstar) %*% xstar) %*% t(xstar) %*% yhat
    new = xstar %*% theta
    list(betahat = theta, yhat = new)
}
