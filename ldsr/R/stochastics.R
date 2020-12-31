#' One LDS replicate
#'
#' Generate a single stochastic time series from an LDS model
#' @param rep.num The ID number of the replicate
#' @param theta A list of parameters: A, B, C, D, Q, R, x0, v0
#' @inheritParams LDS_EM
#' @param mu Mean of the log-transformed streamflow process
#' @param years The years of the study horizon
#' @param exp.trans Whether exponential transformation back to the streamflow space is required.
#' If TRUE, both Y and Q are returned, otherwise only Y.
#' @return A data.table. The first column is the years of the study horizon, as supplied by `year`.
#' Subsequent columns are `simX`, `simY`, and `simQ` which are the simulated catchment state (X),
#' log-transformed and centralized flow (Y) and flow (Q). The last column is the replicate ID number.
#' @examples
#' # Learn theta
#' one_LDS_rep(1, theta, t(NPpc), t(NPpc), 1200:2012, mu = mean(log(NPannual$Qa)))
#' @export
one_LDS_rep <- function(rep.num, theta, u = NULL, v = NULL, years, mu = 0, exp.trans = TRUE) {

    n <- length(years)
    sim.X <- matrix(0, 1, n + 1)
    sim.Y <- matrix(0, 1, n)
    sim.X[, 1] <- stats::rnorm(1, 0, sqrt(theta$V1))

    q <- t(stats::rnorm(n, 0, sqrt(theta$Q)))
    r <- t(stats::rnorm(n, 0, sqrt(theta$R)))

    if (is.null(u)) {
        for (t in 1:n) {
            sim.X[, t + 1] <- theta$A %*% sim.X[,t] + q[,t]
            sim.Y[, t    ] <- theta$C %*% sim.X[,t] + r[,t]
        }
    } else {
        for (t in 1:n) {
            sim.X[, t + 1] <- theta$A %*% sim.X[,t] + theta$B %*% u[, t] + q[t]
            sim.Y[, t    ] <- theta$C %*% sim.X[,t] + theta$D %*% v[, t] + r[t]
        }
    }

    if (exp.trans) sim.Q <- exp(sim.Y + mu) else sim.Q <- sim.Y + mu
    out <- data.table(year = years,
                      simX = as.vector(sim.X[, 1:n]),
                      simY = as.vector(sim.Y),
                      simQ = as.vector(sim.Q),
                      rep  = rep(rep.num, n))
    out
}

#' Multiple LDS replicates
#'
#' Generate multiple stochastic time series from an LDS model. This is a convenient wrapper for [one_LDS_rep].
#' @inheritParams one_LDS_rep
#' @param num.reps The number of stochastic replicates#'
#' @return Same as [one_LDS_rep], but the data.table consists of multiple replicates.
#' @examples
#' LDS_rep(theta, t(NPpc), t(NPpc), 1200:2012, num.reps = 10, mu = mean(log(NPannual$Qa)))
#' @export
LDS_rep <- function(theta, u = NULL, v = NULL, years, num.reps = 100, mu = 0, exp.trans = TRUE) {

    rbindlist(lapply(seq_len(num.reps),
                     function(i) one_LDS_rep(i, theta, u, v, years, mu, exp.trans)))

}
