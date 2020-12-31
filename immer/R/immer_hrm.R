## File Name: immer_hrm.R
## File Version: 0.962

#####################################################
# Hierarchical rater model Patz et al. (2002)
immer_hrm <- function( dat, pid, rater,
                    iter,  burnin, N.save=3000, prior=NULL,
                    est.a=FALSE, est.sigma=TRUE, est.mu=FALSE, est.phi="a", est.psi="a",
                    MHprop=NULL, theta_like=seq( -10, 10, len=30), sigma_init=1, print_iter=20 )
{

    time <- list( "start"=Sys.time() )
    useRcpp <- TRUE
    CALL <- match.call()

    sd_init <- sigma_init

    #************************************
    # rater and person identifiers
    res0 <- immer_identifiers_relabel( dat=dat, pid=pid, rater=rater )
    rater <- res0$rater
    pid <- res0$pid
    dat <- res0$dat
    pid_unique <- res0$pid_unique
    rater_unique <- res0$rater_unique
    item0 <- res0$item0
    rater_pars0 <- res0$rater_pars0

    #************************************
    # inits
    # inits theta parameters
    theta <- inits_theta_1dim( dat=dat, pid=pid, eps=.05, sd_init=sd_init )
    N <- length(theta)

    # parameters settings
    # Options are 'n' (no estimation), 'e' (set all parameters equal to each other),
    # 'i' (item wise estmation), 'r' (rater wise estimation) and
    # 'a' (all parameters are estimated independently from each other).
    est_settings <- list( est.a=est.a, est.sigma=est.sigma, est.mu=est.mu,
                        est.phi=est.phi, est.psi=est.psi )

    # inits item parameters
    res0 <- inits_itempars( dat=dat, prior=prior, sd_init=sd_init )
    b <- res0$b
    a <- res0$a
    K <- res0$K
    I <- res0$I
    maxK <- res0$maxK

    # inits rater parameters HRM
    res0 <- inits_raterpars_hrm( rater=rater, I=I, est_settings=est_settings )
    phi <- res0$phi
    psi <- res0$psi

    R <- res0$R
    dat <- as.matrix(dat)
    dat_ind <- 1 - is.na(dat)
    xi_ind <- 1 * ( rowsum( 1 - is.na(dat), pid ) > 0 )
    ND <- nrow(dat)

    #*********************************************
    # prior parameters
    prior <- prior_hrm( prior=prior, b=b, a=a, phi=phi, est_settings=est_settings, sd_init=sd_init )
    sigma <- sqrt( prior$sigma2$sig02 )
    mu <- prior$mu$M

    #**********************************************
    # objects for saving traces

    BB1 <- iter - burnin
    save_list <- seq( burnin+1, iter )

    if ( N.save > BB1 ){ BB <- BB1 }
    if ( N.save <=BB1 ){
        h1 <- ceiling( BB1 / N.save )
        save_list <- seq( burnin+1, iter, h1)
        BB <- length(save_list)
    }
    psiM <- phiM <- array( NA, dim=c(I,R,BB) )
    bM <- array(NA, dim=c(I,K,BB) )
    aM <- matrix( NA, nrow=I, ncol=BB )
    sigmaM <- rep( NA, BB )
    muM <- rep( NA, BB )
    devM <- rep(NA, BB)
    person <- data.frame( "pid"=pid_unique, "NSamp"=0, "EAP"=0, "SD.EAP"=0 )

    #**********************************************
    # Metropolis-Hastings tuning
    MHprop <- MHprop_hrm( MHprop=MHprop, b=b, a=a, phi=phi, theta=theta, iter=iter, burnin=burnin )

    #********************************************
    #  ITERATIONS
    it <- 0
    bb <- 1

    eps <- 1E-20
    eps11 <- 1E-7
    dat <- as.matrix(dat)
    dat_ind <- as.matrix(dat_ind)
    maxcat <- max(maxK)
    b <- as.matrix(b)
    phi <- as.matrix(phi)
    psi <- as.matrix(psi)
    xi_ind <- as.matrix(xi_ind)

    mcmc_start_time <- Sys.time()

    for ( it in seq(1,iter) ){

        #**** sample xsi
        xi <- sampling_hrm_xi( dat=dat, theta=theta, b=b, a=a, phi=phi, psi=psi, K=K, pid=pid,
                    rater=rater, ND=ND, dat_ind=dat_ind, N=N, I=I, maxK=maxK, useRcpp=useRcpp, xi_ind=xi_ind )

        #**** sample b
        res0 <- sampling_hrm_b( xi=xi, xi_ind=xi_ind, b=b, a=a, maxK=maxK, prior=prior, MHprop=MHprop, I=I,
                        theta=theta, useRcpp=useRcpp )
        b <- res0$b
        MHprop <- res0$MHprop

        #**** sample a
        res0 <- sampling_hrm_a( xi=xi, xi_ind=xi_ind, b=b, a=a, maxK=maxK, prior=prior, MHprop=MHprop,
                        I=I, theta=theta, useRcpp=useRcpp )
        a <- res0$a
        MHprop <- res0$MHprop

        #**** sample phi
        res0 <- sampling_hrm_phi( dat=dat, dat_ind=dat_ind, maxK=maxK, R=R, rater=rater, pid=pid,
                        phi=phi, psi=psi, prior=prior, MHprop=MHprop, I=I, xi=xi, useRcpp=useRcpp,
                        est_settings=est_settings )
        phi <- res0$phi
        MHprop <- res0$MHprop

        #**** sample psi
        res0 <- sampling_hrm_psi( dat=dat, dat_ind=dat_ind, maxK=maxK, R=R, rater=rater, pid=pid,
                        phi=phi, psi=psi, prior=prior, MHprop=MHprop, I=I, xi=xi, useRcpp=useRcpp,
                        est_settings=est_settings )
        psi <- res0$psi
        MHprop <- res0$MHprop

        #**** sample theta
        mu_theta <- rep( mu,N)
        SD_theta <- rep( sigma,N)
        res0 <- sampling_hrm_theta_1dim( theta=theta, N=N, I=I, maxK=maxK, a=a, b=b, xi=xi, xi_ind=xi_ind,
                        dat=dat, dat_ind=dat_ind, pid=pid, MHprop=MHprop, mu_theta=mu_theta,
                        SD_theta=SD_theta, useRcpp=useRcpp, eps=1E-20 )
        theta <- res0$theta
        MHprop <- res0$MHprop

        # sampling of mu
        mu <- sampling_hrm_mu_1dim( theta=theta, prior=prior, N=N )

        # sampling of sigma
        sigma <-  sqrt( immer_mcmc_draw_variance( N=1, w0=prior$sigma2$w0,
                        sig02=prior$sigma2$sig02, n=N, sig2=stats::var(theta) ) )

        #-------------- OUTPUT
        if ( it %in% save_list ){
            bM[,, bb ] <- b
            aM[,bb] <- a
            phiM[,,bb] <- phi
            psiM[,,bb] <- psi
            muM[bb] <- mu
            sigmaM[bb] <- sigma
            person$NSamp <- person$NSamp + 1
            person$EAP <- person$EAP + theta
            person$SD.EAP <- person$SD.EAP + theta^2
            ND <- nrow(dat)
            #***
            # calculate deviance
            if (FALSE){
                dev <- 1
                theta0 <- theta[ pid ]
                for (ii in 1:I){
                    K_ii <- maxK[ii]
                    pr_tot <- 0
                    for (hh in seq(0,K_ii) ){
                        x0 <- rep(hh,ND)
                        pr1_hh <- probs_gpcm( x=x0, theta=theta0, b=as.numeric(b[ii,]), a=a[ii],
                                    K=K_ii, useRcpp=FALSE)
                        pr2_x <- probs_hrm( x=dat[,ii], xi=x0, phi=phi[ ii,rater ],
                                    psi=psi[ ii, rater ], K=K_ii, useRcpp=FALSE )
                        pr_tot <- pr_tot + pr1_hh * pr2_x
                    }
                    dev <- ifelse( dat_ind[,ii]==1, dev*pr_tot, dev )
                }
                dev <- rowsum( log( dev + 1E-200 ), pid )
                dev <- sum( dev[,1] )
                devM[bb] <- dev
            }
            bb <- bb + 1
        }


        #------- update MH parameters
        if ( sum( MHprop$ITER_refreshing %in% it ) > 0  ){
            MHprop <- MHprop_refresh( MHprop=MHprop )
        }

        #------ print iteration progress
        res <- immer_hrm_verbose_iterations( iter=iter, it=it, print_iter=print_iter,
                    mcmc_start_time=mcmc_start_time )

    }
    #------------- end iterations

    #**********************************************
    # arrange output

    # item parameters
    item <- data.frame( item0, "a"=rowMeans( aM ) )
    for (kk in 1:K){
        item[, paste0("b", kk) ] <- rowMeans( bM[,kk,,drop=FALSE] )
    }

    # rater parameters
    rater_pars <- NULL
    for (ii in 1:I){
        dfr <- data.frame(
                "phi"=apply( phiM[ii,,], 1, mean, na.rm=TRUE ),
                "psi"=apply( psiM[ii,,], 1, mean, na.rm=TRUE )
                    )
        rater_pars <- rbind( rater_pars, dfr )
    }
    rater_pars <- data.frame( rater_pars0, rater_pars )

    # person parameters
    person$EAP <- person$EAP / person$NSamp
    person$SD.EAP <- sqrt( ( person$SD.EAP - person$NSamp * person$EAP^2  ) /
                                person$NSamp  )

    # EAP reliability
    EAP.rel <- stats::var( person$EAP ) /
                        ( stats::var( person$EAP ) + mean( person$SD.EAP^2 ) )
    # information criteria
    ic <- list( N=N, I=I, R=R, ND=nrow(dat), maxK=maxK, K=K )

    # list with all traces
    traces <- list( b=bM, a=aM, phi=phiM, psi=psiM, mu=muM,
                            sigma=sigmaM, deviance=devM )

    attr( traces, "NSamples" ) <- BB
    attr(traces, "burnin") <- burnin
    attr( traces, "iter" ) <- iter

    # collect traces and produce MCMC summary
    res11 <- immer_collect_traces( traces=traces, est_settings=est_settings )
    summary.mcmcobj <- sirt::mcmc.list.descriptives( res11$mcmcobj )

    #******
    # extract estimated parameters
    est_pars <- list()
    est_pars$a <- item$a
    est_pars$b <- item[, paste0("b", seq(1,K) ) ]
    items <- paste0( item$item )
    phi <- matrix( NA, nrow=I, ncol=R)
    rownames(phi) <- items
    psi <- phi
    item_index <- match( paste(rater_pars$item), items )
    phi[ cbind( item_index, rater_pars$rid ) ] <- rater_pars$phi
    psi[ cbind( item_index, rater_pars$rid ) ] <- rater_pars$psi
    est_pars$phi <- phi
    est_pars$psi <- psi
    est_pars$mu <- mean(muM)
    est_pars$sigma <- mean(sigmaM)

    #****
    # compute likelihood and posterior
    res_ll <- loglik_hrm( dat=dat, dat_ind=dat_ind, est_pars=est_pars, theta_like=theta_like,
                            rater=rater, pid=pid, maxK=maxK )

    # add information criteria
    ic$dev <- -2*res_ll$ll
    ic <- immer_ic_hrm( ic=ic, summary.mcmcobj=summary.mcmcobj )

    #*****
    # output list
    time$end <- Sys.time()
    res <- list( person=person, item=item, rater_pars=rater_pars,
                est_pars=est_pars, sigma=est_pars$sigma, mu=est_pars$mu,
                mcmcobj=res11$mcmcobj, summary.mcmcobj=summary.mcmcobj, dat=dat, pid=pid,
                rater=rater, EAP.rel=EAP.rel, ic=ic, f.yi.qk=res_ll$f.yi.qk,
                f.qk.yi=res_ll$f.qk.yi, theta_like=theta_like, pi.k=res_ll$pi.k,
                like=res_ll$ll, traces=traces, MHprop=MHprop, prior=prior,
                est_settings=est_settings, N.save=BB, iter=iter, burnin=burnin, time=time,
                CALL=CALL )
    res$description <- "Function immer_hrm() | Hierarchical Rater Model (Patz et al., 2002)"
    class(res) <- "immer_hrm"
    return(res)
}
#############################################################################
