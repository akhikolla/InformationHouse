###############################################################################
#
#    grpSLOPE: Group SLOPE (Group Sorted L1 Penalized Estimation)
#    Copyright (C) 2016 Alexej Gossmann
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
###############################################################################

# lambdas of Theorem 2.5 and equation (2.16) in Brzyski et. al. (2016)
lambdaChiOrtho <- function(fdr, n.group, group.sizes, wt, method) {
  lambda.max <- rep(NA, n.group)
  lambda.min <- rep(NA, n.group)

  for (i in 1:n.group) {
    qchi.seq <- rep(NA, n.group)
    for (j in 1:n.group) {
      qchi.seq[j] <- sqrt(qchisq(1 - fdr*i/n.group, df=group.sizes[j])) / wt[j]
    }
    lambda.max[i] <- max(qchi.seq)
    lambda.min[i] <- min(qchi.seq)
  }

  # stop here if method is "max"
  if (method=="max") return(lambda.max)

  cdfMean <- function(x) {
    pchi.seq <- rep(NA, n.group)
    for (i in 1:n.group) {
      pchi.seq[i] <- pchisq((wt[i]*x)^2, df=group.sizes[i])
    }
    return(mean(pchi.seq))
  }

  lambda.mean <- rep(NA, n.group)
  for (k in 1:n.group) {
    if (lambda.min[k] == lambda.max[k]) {
      lambda.mean[k] <- lambda.max[k]
    } else {
      # compute inverse of cdfMean at 1-fdr*k/n.group
      cdfMean.inv <- uniroot(function(y) (cdfMean(y) - (1-fdr*k/n.group)),
                             lower = lambda.min[k], upper = lambda.max[k], extendInt="yes")
      lambda.mean[k] <- cdfMean.inv$root
    }
  }

  return(lambda.mean)
}

# Procedure 6 in Brzyski et. al. (2016)
lambdaChiEqual <- function(fdr, n.obs, n.group, m, w) {
  lambda.chi    <- rep(NA, n.group)
  lambda.chi[1] <- sqrt(qchisq(1 - fdr / n.group, df=m)) / w

  for (i in 2:n.group) {
    # prevent division by 0 or sqrt of a negative number later on
    if ( (n.obs - m*(i-1) - 1) <= 0 ) {
      stop("Corrected lambdas cannot be computed unless groups sizes are small enough compared to sample size.")
    }

    s <- (n.obs - m*(i-1)) / n.obs + (w^2 * sum(lambda.chi[1:(i-1)]^2)) / (n.obs - m*(i-1) - 1)
    s <- sqrt(s)
    lambda.tmp <- (s/w) * sqrt(qchisq(1 - fdr * i / n.group, df=m))
    if (lambda.tmp <= lambda.chi[i-1]) {
      lambda.chi[i] <- lambda.tmp
    } else {
      lambda.chi[i:n.group] <- lambda.chi[i-1]
      break
    }
  }

  return(lambda.chi)
}

# Procedure 1 in Brzyski et. al. (2016)
lambdaChiMean <- function(fdr, n.obs, n.group, group.sizes, wt) {
  lambda.chi.mean <- rep(NA, n.group)

  cdfMean <- function(x) {
    pchi.seq <- rep(NA, n.group)
    for (i in 1:n.group) {
      pchi.seq[i] <- pchisq((wt[i]*x)^2, df=group.sizes[i])
    }
    return(mean(pchi.seq))
  }

  # get upper and lower bounds for lambda.chi.mean[1]
  qchi.seq <- rep(NA, n.group)
  for (j in 1:n.group) {
    qchi.seq[j] <- sqrt(qchisq(1 - fdr/n.group, df=group.sizes[j])) / wt[j]
  }
  upperchi <- max(qchi.seq)
  lowerchi <- min(qchi.seq)

  if (upperchi == lowerchi) {
    lambda.chi.mean[1] <- upperchi
  } else {
    lambda.chi.mean[1] <- uniroot(function(y) (cdfMean(y) - (1-fdr/n.group)),
                                  lower = lowerchi, upper = upperchi, extendInt="yes")$root
  }

  # get lambda.chi.mean[2:n.group]
  for (i in 2:n.group) {
    s <- rep(NA, n.group)
    for (j in 1:n.group) {
      # prevent division by 0 or sqrt of a negative number later on
      if ((n.obs - group.sizes[j]*(i-1) - 1) <= 0) {
        stop("Corrected lambdas cannot be computed unless groups sizes are small enough compared to sample size.")
      }

      s[j] <- (n.obs - group.sizes[j]*(i-1)) / n.obs + 
        (wt[j]^2 * sum(lambda.chi.mean[1:(i-1)]^2)) / (n.obs - group.sizes[j]*(i-1) - 1)
      s[j] <- sqrt(s[j])
    }

    cdfMean <- function(x) {
      pchi.seq <- rep(NA, n.group)
      for (j in 1:n.group) {
        pchi.seq[j] <- pchisq((wt[j]/s[j] * x)^2, df=group.sizes[j])
      }
      return(mean(pchi.seq))
    }

    cdfMean.inv <- uniroot(function(y) (cdfMean(y) - (1-fdr*i/n.group)),
                           lower = 0, upper = upperchi, extendInt="upX")$root

    if (cdfMean.inv <= lambda.chi.mean[i-1]) {
      lambda.chi.mean[i] <- cdfMean.inv 
    } else {
      lambda.chi.mean[i:n.group] <- lambda.chi.mean[i-1]
      break
    }
  }

  return(lambda.chi.mean)
}
