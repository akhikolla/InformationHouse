fix.MH = function(lo, hi, old1, new1, jump) {
    jold = pnorm(hi, mean = old1, sd = jump) - pnorm(lo, mean = old1, sd = jump)
    jnew = pnorm(hi, mean = new1, sd = jump) - pnorm(lo, mean = new1, sd = jump)
    log(jold) - log(jnew)  #add this to pnew
}
#----------------------------------------------------------------------------------------------
tnorm <- function(n, lo, hi, mu, sig) {
    # generates truncated normal variates based on cumulative normal distribution normal truncated lo and hi

    if (length(lo) == 1 & length(mu) > 1)
        lo <- rep(lo, length(mu))
    if (length(hi) == 1 & length(mu) > 1)
        hi <- rep(hi, length(mu))

    q1 <- pnorm(lo, mu, sig)  #cumulative distribution
    q2 <- pnorm(hi, mu, sig)  #cumulative distribution

    z <- runif(n, q1, q2)
    z <- qnorm(z, mu, sig)
    z[z == -Inf] <- lo[z == -Inf]
    z[z == Inf] <- hi[z == Inf]
    z
}
#----------------------------------------------------------------------------------------------
acceptMH <- function(p0, p1, x0, x1, BLOCK) {
    # accept for M, M-H if BLOCK, then accept as a block, otherwise, accept individually

    nz <- length(x0)  #no. to accept
    if (BLOCK)
        nz <- 1

    a <- exp(p1 - p0)  #acceptance PR
    z <- runif(nz, 0, 1)
    keep <- which(z < a)

    if (BLOCK & length(keep) > 0)
        x0 <- x1
    if (!BLOCK)
        x0[keep] <- x1[keep]
    accept <- length(keep)

    list(x = x0, accept = accept)
}
#-------------------------------
print.adapt = function(parmAccept, parmJump, accept.output, printBool = FALSE) {
    accept1 = parmAccept
    jump1 = parmJump

    for (k in 1:length(accept1)) {
        z = accept1[[k]]/accept.output
        if (printBool) {
            print(names(accept1)[k])
            print(mean(z))
            print(mean(jump1[[k]]))
        }
    }

    for (k in 1:length(jump1)) {
        cond = (accept1[[k]]/accept.output) > 0.4 & jump1[[k]] < 100
        jump1[[k]][cond] = jump1[[k]][cond] * 2
        cond = (accept1[[k]]/accept.output) < 0.2 & jump1[[k]] > 0.001
        jump1[[k]][cond] = jump1[[k]][cond] * 0.5
        accept1[[k]][] = 0
    }

    return(list(jump1 = jump1, accept1 = accept1))
}


aggregate.data = function(dat) {
    dat1 = dat[order(dat$loc.id), ]
    nloc = max(dat1$loc.id)

    ind = which(colnames(dat1) == "loc.id")

    res = matrix(NA, nloc, ncol(dat1) - 1)
    n = rep(NA, nloc)
    for (i in 1:nloc) {
        cond = dat1$loc.id == i
        dat.tmp = dat1[cond, -ind]
        n[i] = nrow(dat.tmp)
        if (n[i] > 1)
            dat.tmp = colSums(dat.tmp)
        if (n[i] == 1)
            dat.tmp = as.numeric(dat.tmp)
        res[i, ] = dat.tmp
    }
    colnames(res) = colnames(dat1[, -ind])
    loc.id = unique(dat1$loc.id)
    list(dat = res, loc.id = loc.id, n = n)
}

fix.probs = function(probs) {
    cond = probs < 1e-05
    probs[cond] = 1e-05
    cond = probs > 0.99999
    probs[cond] = 0.99999
    probs
}
#----------------------------
get.logl = function(theta, phi, y, nmat) {
    prob = fix.probs(theta %*% phi)
    dbinom(y, size = nmat, prob = prob, log = T)
}

#----------------------------------------------------
update.omegaRemote = function(remote, param, jump, ncommun, nbands, ndig.values, a.omega, b.omega) {
    omega.orig = omega.old = param$omega
    tmp = tnorm(nbands * ncommun, lo = 0, hi = 1, mu = omega.old, sig = jump)
    novos = matrix(tmp, ncommun, nbands)

    tmp = fix.MH(lo = 0, hi = 1, omega.old, novos, jump)
    ajuste = matrix(tmp, ncommun, nbands)

    prior.old = matrix(dbeta(omega.old, a.omega, b.omega, log = T), ncommun, nbands)
    prior.new = matrix(dbeta(novos, a.omega, b.omega, log = T), ncommun, nbands)

    for (i in 1:ncommun) {
        omega.new = omega.old
        omega.new[i, ] = novos[i, ]

        prob.old = param$theta %*% omega.old
        llk.old = colSums(dbinom(remote, size = ndig.values, prob = prob.old, log = T))

        prob.new = param$theta %*% omega.new
        llk.new = colSums(dbinom(remote, size = ndig.values, prob = prob.new, log = T))

        k = acceptMH(llk.old + prior.old[i, ], llk.new + prior.new[i, ] + ajuste[i, ], omega.old[i, ], omega.new[i, ], F)
        omega.old[i, ] = k$x
    }
    list(omega = omega.old, accept = omega.old != omega.orig)
}

update.thetaRemote = function(remote, param, jump, ncommun, nloc, ndig.values) {
    v.orig = v.old = param$v
    tmp = tnorm(nloc * (ncommun - 1), lo = 0, hi = 1, mu = v.old[, -ncommun], sig = jump[, -ncommun])
    novos = cbind(matrix(tmp, nloc, ncommun - 1), 1)
    ajuste = matrix(fix.MH(lo = 0, hi = 1, v.old, novos, jump), nloc, ncommun)

    prior.old = matrix(dbeta(v.old, 1, param$gamma, log = T), nloc, ncommun)
    prior.new = matrix(dbeta(novos, 1, param$gamma, log = T), nloc, ncommun)

    for (j in 1:(ncommun - 1)) {
        # last column has to be 1
        v.new = v.old
        v.new[, j] = novos[, j]

        theta.old = convertSBtoNormal(vmat = v.old, ncol = ncommun, nrow = nloc, prod = rep(1, nloc))
        theta.new = convertSBtoNormal(vmat = v.new, ncol = ncommun, nrow = nloc, prod = rep(1, nloc))

        # contribution from reflectance data
        pold = theta.old %*% param$omega
        pnew = theta.new %*% param$omega
        p1.old = rowSums(dbinom(remote, size = ndig.values, pold, log = T))
        p1.new = rowSums(dbinom(remote, size = ndig.values, pnew, log = T))

        k = acceptMH(p1.old + prior.old[, j], p1.new + prior.new[, j] + ajuste[, j], v.old[, j], v.new[, j], F)
        v.old[, j] = k$x
    }
    theta = convertSBtoNormal(vmat = v.old, ncol = ncommun, nrow = nloc, prod = rep(1, nloc))
    list(theta = theta, v = v.old, accept = v.old != v.orig)
}

#-------------------------------------------------------------------------------
# FULL CONDITIONAL DISTRIBUTION FUNCTIONS
#----------------------------
update.phiAbundanceSB = function(param, jump, ncomm, nspp, y, nmat, a.phi, b.phi) {
    phi.orig = phi.old = param$phi
    proposed = matrix(tnorm(nspp * ncomm, lo = 0, hi = 1, mu = phi.old, sig = jump), ncomm, nspp)

    for (i in 1:ncomm) {
        phi.new = phi.old
        phi.new[i, ] = proposed[i, ]
        adj = fix.MH(lo = 0, hi = 1, old1 = phi.old[i, ], new1 = phi.new[i, ], jump = jump[i, ])

        prob.old = get.logl(theta = param$theta, phi = phi.old, y = y, nmat = nmat)
        prob.new = get.logl(theta = param$theta, phi = phi.new, y = y, nmat = nmat)
        pold = colSums(prob.old) + dbeta(phi.old[i, ], a.phi, b.phi, log = T)
        pnew = colSums(prob.new) + dbeta(phi.new[i, ], a.phi, b.phi, log = T)
        k = acceptMH(p0 = pold, p1 = pnew + adj, x0 = phi.old[i, ], x1 = phi.new[i, ], BLOCK = F)
        phi.old[i, ] = k$x
    }
    list(phi = phi.old, accept = phi.orig != phi.old)
}
#-----------------------------
update.thetaAbundanceSB = function(param, jump, nloc, ncomm, y, nmat, gamma) {
    v.orig = v.old = param$vmat
    tmp = tnorm(nloc * (ncomm - 1), lo = 0, hi = 1, mu = v.old[, -ncomm], sig = jump[, -ncomm])
    novos = cbind(matrix(tmp, nloc, ncomm - 1), 1)
    ajuste = matrix(fix.MH(lo = 0, hi = 1, old1 = v.old, new1 = novos, jump = jump), nloc, ncomm)

    prior.old = matrix(dbeta(v.old, 1, gamma, log = T), nloc, ncomm)
    prior.new = matrix(dbeta(novos, 1, gamma, log = T), nloc, ncomm)

    for (j in 1:(ncomm - 1)) {
        # last column has to be 1
        v.new = v.old
        v.new[, j] = novos[, j]

        theta.old = convertSBtoNormal(vmat = v.old, ncol = ncomm, nrow = nloc, prod = rep(1, nloc))
        theta.new = convertSBtoNormal(vmat = v.new, ncol = ncomm, nrow = nloc, prod = rep(1, nloc))

        pold = get.logl(theta = theta.old, phi = param$phi, y = y, nmat = nmat)
        pnew = get.logl(theta = theta.new, phi = param$phi, y = y, nmat = nmat)
        p1.old = rowSums(pold)
        p1.new = rowSums(pnew)

        k = acceptMH(p0 = p1.old + prior.old[, j], p1 = p1.new + prior.new[, j] + ajuste[, j], x0 = v.old[, j], x1 = v.new[, j], F)
        v.old[, j] = k$x
    }
    theta = convertSBtoNormal(vmat = v.old, ncol = ncomm, nrow = nloc, prod = rep(1, nloc))
    list(theta = theta, v = v.old, accept = v.old != v.orig)
}


get.theta=function(nlk,gamma,ncomm,nloc){
  vmat=matrix(NA,nloc,ncomm)
  for (i in 1:(ncomm-1)){
    if (i==(ncomm-1)) cumsoma=nlk[,ncomm]
    if (i< (ncomm-1)) cumsoma=rowSums(nlk[,(i+1):ncomm])
    vmat[,i]=rbeta(nloc,nlk[,i]+1,cumsoma+gamma)
  }
  vmat[,ncomm]=1
  convertVtoTheta(vmat,rep(1,nloc))
}
