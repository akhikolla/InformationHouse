em_fit <- function(logfrailtypar, dist, pvfm,
                   Y, Xmat, # id,  # this is some data stuff
                   atrisk, # a list with stuff that won't change in the EM
                   basehaz_line ,  # need for log-likelihood
                   mcox = list(),
                   Cvec, lt = FALSE, Cvec_lt, # we need to start somewhere with the Cvec (E step comes first)
                   em_control, se,
                   return_loglik = TRUE
) {


  pars <- dist_to_pars(dist, logfrailtypar, pvfm)
  if (isTRUE(em_control$verbose)) {
    print(paste0(#"dist=", pars$dist,
      "logfrailtypar= ", logfrailtypar,
      " / alpha=", pars$alpha,
      " / bbeta=", pars$bbeta))
  }

  if(logfrailtypar < -100) warning("theta virtually 0; try another starting value")

  if(length(Xmat)==0) {
    g_x <- matrix(rep(0, nrow(Y)),ncol = 1)
  } else {
    g_x <- t(mcox$coefficients %*% t(Xmat))
  }

  # if the logfrailtypar is large (i.e. frailty variance 0) then just return the Cox likelihood
  if(logfrailtypar > log(em_control$upper_tol)) {

    #message("Frailty parameter very large, frailty variance close to 0")
    loglik <- mcox$loglik[length(mcox$loglik)]


    if(isTRUE(return_loglik)) {
      if(isTRUE(em_control$verbose)) print(paste("loglik = ",loglik))
      return(-loglik)
    }

  }


  loglik_old = -Inf
  ncycles <- 0

  convergence <- FALSE

  while(!isTRUE(convergence)) {


    if(isTRUE(em_control$fast_fit)) {
      e_step_val <- fast_Estep(Cvec,
                               Cvec_lt,
                               atrisk$nev_id,
                               alpha = pars$alpha,
                               bbeta = pars$bbeta,
                               pvfm = pvfm,
                               dist = pars$dist)
    } else {
      e_step_val <- Estep(Cvec,
                          Cvec_lt,
                          atrisk$nev_id,
                          alpha = pars$alpha,
                          bbeta = pars$bbeta,
                          pvfm = pvfm,
                          dist = pars$dist)
    }

     logz <- log((e_step_val[,1] / e_step_val[,2])[atrisk$order_id])


     if(!is.null(atrisk$strats))
       loglik <-
       sum(mapply(function(bhz, death) sum(log(bhz[death])),
                  basehaz_line ,
                  atrisk$death)) +
       sum(g_x[Y[,3] == 1]) +
       sum(e_step_val[,3]) +
       sum(Y[,3]) -
       sum((do.call(c, atrisk$nevent) * log(do.call(c, atrisk$nevent)))[do.call(c, atrisk$nevent) > 0]) else
         loglik <-  sum((log(basehaz_line) + g_x)[Y[,3] == 1]) + sum(e_step_val[,3]) +
       sum(Y[,3]) - sum((atrisk$nevent * log(atrisk$nevent))[atrisk$nevent > 0])



    # if this happens, then something is going very wrong
    if(loglik < loglik_old - em_control$lik_tol)
      warning(paste0("likelihood decrease of ", loglik - loglik_old ))

    if(abs(loglik - loglik_old) < em_control$eps) break

    loglik_old <- loglik


    # browser()
    # working here
    #
    # dlist <- survival::survreg.distributions$exponential
    #
    # # this has a trans aspect
    # tranfun <- dlist$trans
    # exactsurv <- (Y[, ncol(Y)] == 1)

    # if (any(exactsurv)) {
    #     logcorrect <- sum(log(dlist$dtrans(Y[exactsurv, 1])))                                                                 1])))
    # }

   #  scale <- dlist$scale
   #
   #  dlist <- survival::survreg.distributions[[dlist$dist]]
   #  Yss <- cbind(dlist$trans(Y[,1]), Y[,2])
   # mcox_par <- survival::survreg.fit(x = Xmat, y = Y, strata = atrisk$strats, offset = logz, init = NULL,
   #                              control = survival::survreg.control(), weights = NULL,
   #                              scale = scale,
   #                              dist = dlist)

    mcox <- survival::agreg.fit(x = Xmat, y = Y, strata = atrisk$strats, offset = logz, init = NULL,
                                control = survival::coxph.control(), weights = NULL,
                                method = "breslow", rownames = NULL)

    ncycles <- ncycles + 1

    # NOTE: this is what linear.predictors actually is:
    # exp(mcox$coefficients * (Xmat - mean(Xmat)) + logz)

    if(length(Xmat) == 0) {
      lp <- mcox$linear.predictors
      g_x <- t(matrix(rep(0, length(mcox$linear.predictors)), nrow = 1))
    } else {
      g_x <- t(mcox$coefficients %*% t(Xmat))
      lp <- mcox$linear.predictors + as.numeric(t(mcox$coefficients) %*% mcox$means)

    }

    explp <- exp(lp)


    if(!is.null(atrisk$strats)) {
      explp_str <- split(explp, atrisk$strats)

      nrisk <- mapply(FUN = function(explp, y) rowsum_vec(explp, y, max(y)),
                      explp_str,
                      atrisk$ord_tstop_str,
                      SIMPLIFY = FALSE)

      # nrisk  <- mapply(FUN = function(explp, y) rev(cumsum(rev(rowsum(explp, y[,2])))),
      #                  split(explp, atrisk$strats),
      #                  split.data.frame(Y, atrisk$strats),
      #                  SIMPLIFY = FALSE)

      esum <- mapply(FUN = function(explp, y) rowsum_vec(explp, y, max(y)),
                     explp_str,
                     atrisk$ord_tstart_str,
                     SIMPLIFY = FALSE)

      # esum  <-  mapply(FUN = function(explp, y) rev(cumsum(rev(rowsum(explp, y[,1])))),
      #                  split(explp, atrisk$strats),
      #                  split.data.frame(Y, atrisk$strats),
      #                  SIMPLIFY = FALSE)

      nrisk_b  <- mapply(FUN = function(nrisk, esum, indx)  nrisk - c(esum, 0,0)[indx],
                         nrisk ,
                         esum ,
                         atrisk$indx,
                         SIMPLIFY = FALSE)

      haz  <- mapply(FUN = function(nevent, nrisk) nevent/nrisk, # * newrisk,
                     atrisk$nevent,
                     nrisk_b ,
                     SIMPLIFY = FALSE)

      basehaz_line  <- mapply(FUN = function(haz, time_to_stop) haz[time_to_stop],
                              haz ,
                              atrisk$time_to_stop,
                              SIMPLIFY = FALSE)

      cumhaz  <- lapply(haz , cumsum)

      cumhaz_0_line  <- mapply(FUN = function(cumhaz, time_to_stop) cumhaz[time_to_stop],
                               cumhaz ,
                               atrisk$time_to_stop,
                               SIMPLIFY = FALSE)

      cumhaz_tstart  <- mapply(FUN = function(cumhaz, indx2) c(0, cumhaz)[indx2 + 1],
                               cumhaz ,
                               atrisk$indx2,
                               SIMPLIFY = FALSE)

      cumhaz_line  <- mapply(FUN = function(cumhaz_0_line, cumhaz_tstart)
        (cumhaz_0_line - cumhaz_tstart),
        cumhaz_0_line ,
        cumhaz_tstart ,
        SIMPLIFY = FALSE)

      cumhaz_line <- do.call(c, cumhaz_line)[order(atrisk$positions_strata)]


    } else {


      nrisk <- rowsum_vec(explp, atrisk$ord_tstop, max(atrisk$ord_tstop))
      esum <- rowsum_vec(explp, atrisk$ord_tstart, max(atrisk$ord_tstart))
      # nrisk <- rev(cumsum(rev(rowsum(explp, Y[, ncol(Y) - 1]))))
      # esum <- rev(cumsum(rev(rowsum(explp, Y[, 1]))))
      nrisk <- nrisk - c(esum, 0,0)[atrisk$indx]

      haz <- atrisk$nevent/nrisk
      basehaz_line <- haz[atrisk$time_to_stop]
      cumhaz <- cumsum(haz)

      cumhaz_0_line <- cumhaz[atrisk$time_to_stop]

      cumhaz_tstart <- c(0, cumhaz)[atrisk$indx2 + 1]

      cumhaz_line <- (cumhaz_0_line - cumhaz_tstart)

    }


    if(isTRUE(lt)) {

      if(!is.null(atrisk$strats))
        cumhaz_tstart <- do.call(c, cumhaz_tstart)[order(atrisk$positions_strata)]

      Cvec_lt <- rowsum(x = cumhaz_tstart * exp(g_x), atrisk$order_id , reorder = FALSE)

    } else {
      Cvec_lt <- 0 * Cvec
    }


    Cvec <- rowsum(cumhaz_line * exp(g_x), atrisk$order_id, reorder = FALSE)


    if(ncycles > em_control$maxit) {
      warning(paste("did not converge in ", em_control$maxit," iterations." ))
      break
    }
  }


  if(isTRUE(return_loglik)) {
    if(isTRUE(em_control$verbose)) print(paste("loglik = ",loglik))
    return(-loglik)
  }


  # From this point on, the standard errors & return object

  # atrisk

  if(!is.null(atrisk$strats)) {
    time_str <- lapply(split(Y[,2], atrisk$strats), function(x) sort(unique(x)))

    tev_str <- mapply(function(tim, hz) tim[hz>0],
                    time_str,
                    atrisk$nevent,
                    SIMPLIFY = FALSE)

  haz_tev = lapply(haz, function(x) x[x>0])

  tev <- tev_str

  nev_tp_str <- lapply(atrisk$nevent, function(x) x[x!=0])

  } else {
    tev <- atrisk$time[haz > 0]
    haz_tev = haz[haz > 0]
    nev_tp <- atrisk$nevent[atrisk$nevent!=0]

  }


  if(!isTRUE(se)) {
      # Vcov <- mcox$var
    if(!is.null(atrisk$strats)) {
      tev <- tev_str
      len <- length(do.call(c, tev_str))
    } else
      len <- length(tev)

    Vcov <- matrix(NA, len + length(mcox$coefficients),
                   len + length(mcox$coefficients))
      res = list(loglik = loglik,
               tev = tev,
               haz = haz_tev,
               nev_id = atrisk$nev_id,
               Cvec = Cvec,
               frail = e_step_val[,1] / e_step_val[,2],
               coef = mcox$coefficients,
               Vcov = Vcov,
               fitted = g_x + logz,
               cumhaz_line = cumhaz_line)
    return(res)
  }

  # Standard error calculation


  # atrisk$nevent[atrisk$nevent!=0]

  z_elp = exp(lp)
  elp = exp(lp)  / exp(logz)

  # Building E[d2l/dx^2 | Z]


  if(length(Xmat)>0) {
    x <- lapply(apply(Xmat, 1, list), function(x) x[[1]])
    x_z_elp <- Map(function(a,b) a*b, x, z_elp)
    x_z_elp_H0 <- Map(function(a,b,c) a*b*c, x, z_elp, cumhaz_line)
    x_elp_H0 <- Map(function(a,b,c) a*b*c, x, z_elp / exp(logz), cumhaz_line)

    xx <- lapply(x, function(x) x %*% t(x) )
    xx_z_elp_H0 <- Map(function(a,b, c) a * b * c, xx, z_elp, cumhaz_line)
    m_d2l_dgdg <- Reduce("+", xx_z_elp_H0)



    # which lines contain each event time within each strata

    if(!is.null(atrisk$strats)) {

    p1 <- mapply(function(tev, tstart, tstop) lapply(tev, function(tk) which(tstart < tk & tk <= tstop)),
           tev_str,
           split(Y[,1], atrisk$strats),
           split(Y[,2], atrisk$strats),
           SIMPLIFY = FALSE
    )

    # lapply(function(x) lapply(x, function(...) Reduce("+", ...),

    m_d2l_dhdg <- do.call(rbind, lapply(
      mapply(function(xzelp, p1) lapply(p1, function(x) xzelp[x]),
                             split(x_z_elp, atrisk$strats),
                             p1,
                             SIMPLIFY = FALSE
    ),

    function(x) do.call(rbind, lapply(x, function(...) Reduce("+", ...)))))
  } else

    m_d2l_dhdg <-
      do.call(rbind,
              lapply(lapply(
                lapply(tev, function(tk) which(Y[,1] < tk & tk <= Y[,2])),
                function(x) x_z_elp[x]),
                function(...) Reduce("+", ...))
              )


  } else {
    m_d2l_dgdg <- NULL
    m_d2l_dhdg <- NULL
  }


  if(!is.null(atrisk$strats))
    m_d2l_dhdh <- diag(do.call(c, nev_tp_str) / do.call(c,haz_tev)^2) else
      m_d2l_dhdh <- diag(nev_tp/haz_tev^2)



  # Building E[dl/dx (dl/dx)' | Z]

  if(isTRUE(em_control$fast_fit)) {
      z <- e_step_val[,1] / e_step_val[,2]
      zz <- e_step_val[,4]
    } else {
      estep_plusone <- Estep(Cvec,
                             Cvec_lt,
                             atrisk$nev_id+1,
                             alpha = pars$alpha,
                             bbeta = pars$bbeta,
                             pvfm = pvfm,
                             dist = pars$dist)

      # estep_again <- Estep(Cvec,
      #                      Cvec_lt,
      #                      atrisk$nev_id,
      #                      alpha = pars$alpha,
      #                      bbeta = pars$bbeta,
      #                      pvfm = pvfm,
      #                      dist = pars$dist)
      zz <- estep_plusone[,1] /e_step_val[,2]
      z <- e_step_val[,1] / e_step_val[,2]
    }




  if(!is.null(atrisk$strats)) {
  dl1_dh <- do.call(c, nev_tp_str) / do.call(c, haz_tev)

  tl_ord <- mapply(function(tstart, tev) findInterval(tstart, tev) ,
                   split(Y[,1], atrisk$strats),
                   tev_str,
                   SIMPLIFY = FALSE)

  tr_ord <- mapply(function(tstop, tev) findInterval(tstop, tev, left.open = FALSE, rightmost.closed = FALSE) ,
                   split(Y[,2], atrisk$strats),
                   tev_str,
                   SIMPLIFY = FALSE)

  # dl2_dh <- try(inf_mat_match(
  #   tl_ord,
  #   tr_ord,
  #   z_elp,
  #   length(tev)
  # ))

  dl2_dh <- mapply(function(a,b,c,d)
    cumsum_elp(a,b,c,d),
    tl_ord,
    tr_ord,
    split(z_elp, atrisk$strats),
    lapply(tev_str, length),
    SIMPLIFY = FALSE
    )

  elp_to_tev <- mapply(
    function(elp, y1, y2, ord_id, tev)
      lapply(split.data.frame(data.frame(elp, y1 = y1, y2 = y2),
                              ord_id),
             function(dat) cumsum_elp(dat$y1, dat$y2, dat$elp, length(tev))),
    split(elp, atrisk$strats),
    tl_ord,
    tr_ord,
    split(atrisk$order_id, atrisk$strats),
    tev_str,
    SIMPLIFY = FALSE
  )

} else {

  dl1_dh <- nev_tp / haz_tev

  tl_ord <- findInterval(Y[,1], tev)
  tr_ord <- findInterval(Y[,2], tev, left.open = FALSE, rightmost.closed = FALSE)

  dl2_dh <- try(cumsum_elp(
    tl_ord,
    tr_ord,
    z_elp,
    length(tev)
  ))

  elp_to_tev <-  lapply(
    split.data.frame(data.frame(elp,
                                y1 = findInterval(Y[,1], tev),
                                y2 = findInterval(Y[,2], tev, left.open = FALSE, rightmost.closed = FALSE)),
                     atrisk$order_id),
    function(dat) cumsum_elp(dat$y1, dat$y2, dat$elp, length(tev))
  )
}




  if(length(Xmat) > 0) {


    tmp1 <- rowsum(do.call(rbind, x_elp_H0), atrisk$order_id, reorder = FALSE) * sqrt(zz - z^2)
    cor_dg <- Reduce("+",lapply(split(tmp1, 1:nrow(tmp1)), function(x) x %*% t(x)))

    # I_gg_loss <- cor_dg

    if(!is.null(atrisk$strats)) {
    zz_z2_str <- lapply(split(atrisk$order_id, atrisk$strats), function(x) (zz - z^2)[x])

    # same as tmp 1 but now on strata and no sqrt in zz_z2
    tmp2 <- mapply(function(xelph0, ordid, zz_z2) rowsum(do.call(rbind, xelph0) * zz_z2, ordid, reorder = FALSE),
           split(x_elp_H0, atrisk$strats),
           split(atrisk$order_id, atrisk$strats),
           zz_z2_str,
           SIMPLIFY = FALSE
    )

    cor_dg_dh <- t(do.call(rbind,
            mapply(function(elptev,tm2) Reduce("+",
                                       Map(function(a,b) a %*% t(b) ,
                                    elptev, tm2)),
           elp_to_tev,
           lapply(tmp2, function(x) split(x, 1:nrow(x))),
           SIMPLIFY = FALSE
           )
    ))
    } else {
      cor_dg_dh <- t(Reduce("+",
                            Map(function(a,b) a %*% t(b),
                                elp_to_tev,
                                split(
                                  tmp1 *  sqrt(zz - z^2),
                                  1:nrow(tmp1))
                            )
      )
      )
    }

    I_gh_loss <-  cor_dg_dh

  } else {
    cor_dg_dh <- NULL
    cor_dg <- NULL
  }

  if(!is.null(atrisk$strats)) {
  zz_z2_clus <- lapply(
    lapply(
      split(atrisk$order_id, atrisk$strats),
      unique),
    function(x) sqrt(zz - z^2)[x])

  elptev_zzz2 <- mapply(function(elptev, zzz2) Map(function(a,b) a * b, elptev, zzz2),
       elp_to_tev, zz_z2_clus, SIMPLIFY = FALSE)

  # browser()
  cor_dh <- bdiag(
    lapply(elptev_zzz2, function(a)
      {
        m <- matrix(0, length(a[[1]]),  length(a[[1]]))
        m[upper.tri(m, diag = TRUE)] <- sumxxt(a, length(a[[1]]))
        m + t(m) -diag(diag(m), nrow = dim(m)[1],
                            ncol = dim(m)[2])
        }
      )
    )
  } else {
    a <- Map(function(a,b) a * b,
             elp_to_tev,
             sqrt(zz - z^2))
    m <- matrix(0, length(a[[1]]),  length(a[[1]]))
    m[upper.tri(m, diag = TRUE)] <- sumxxt(a, length(a[[1]]))
    cor_dh <- m + t(m) - diag(diag(m))
  }

  I_hh <- m_d2l_dhdh - cor_dh

  if(length(Xmat)>0) {
    I_gg <- m_d2l_dgdg - cor_dg
    I_hg <- m_d2l_dhdg - t(cor_dg_dh)

    if(!is.null(atrisk$strats)) tev <- do.call(c, tev_str)

    Imat <- matrix(0, ncol(Xmat) + length(tev), ncol(Xmat) + length(tev))

    Imat[1:length(mcox$coefficients), 1:length(mcox$coefficients)] <- I_gg

    Imat[(length(mcox$coefficients)+1):nrow(Imat), (length(mcox$coefficients)+1):nrow(Imat)] <- as.matrix(I_hh)

    Imat[1:length(mcox$coefficients), (length(mcox$coefficients)+1):nrow(Imat) ] <- t(I_hg)
    Imat[(length(mcox$coefficients)+1):nrow(Imat), 1:length(mcox$coefficients) ] <- I_hg

  } else Imat <- I_hh


  Vcov = try(chol2inv(chol(Imat)), silent = TRUE)


  if(!is.null(atrisk$strats)) {
    tev <- tev_str
  }


  if(!isTRUE(return_loglik)) {
    res = list(loglik = loglik,
               tev = tev,
               haz = haz_tev,
               nev_id = atrisk$nev_id,
               Cvec = Cvec,
               frail = e_step_val[,1] / e_step_val[,2],
               coef = mcox$coefficients,
               Vcov = Vcov,
               fitted = g_x + logz,
               logz = logz,
               cumhaz_line = cumhaz_line)

    res
  }


}


