## File Name: immer_latent_regression.R
## File Version: 0.37

immer_latent_regression <- function( like, theta=NULL, Y=NULL, group=NULL,
        weights=NULL, conv=1E-5, maxit=200, verbose=TRUE)
{
    time <- list( start=Sys.time() )
    CALL <- match.call()
    X <- Y  # relabelling
    N <- nrow(like)

    #-- theta
    if ( is.null(theta) ){
        theta <- attr(like, "theta")
    }
    TP <- length(theta)

    #-- design matrix regression
    if ( is.null(X) ){
        X <- matrix(1, nrow=N, ncol=1)
        if (! is.null(group)){
            group_u <- unique(group)
            G <- length(group_u)
            X <- matrix(0,nrow=N, ncol=G)
            for (gg in 1:G){
                X[,gg] <- 1 * ( group==group_u[gg] )
            }
            colnames(X) <- paste0("Group", group_u)
        }
    }
    X <- as.matrix(X)
    if ( is.null(colnames(X))){
        colnames(X) <- paste0("X", 1:ncol(X))
    }
    if (is.null(weights)){
        weights <- rep(1,N)
    }

    beta <- rep(0, ncol(X))
    if (is.null(group)){
        group <- rep(1,N)
    }
    group <- match( group, unique(group) )
    G <- length( unique(group) )
    gamma <- rep(1, G)
    ind_group <- list()
    for (gg in 1:G){
        ind_group[[gg]] <- which( group==gg )
    }
    Xw <- X * weights

    XtX <- MASS::ginv( t(X) %*% Xw )
    thetaM <- matrix( theta, nrow=N, ncol=TP, byrow=TRUE)
    mu <- as.vector( X %*% beta )
    sigma <- gamma[ group ]
    iter <- 1
    iterate <- TRUE

    #*** begin iterations
    while (iterate){
        beta0 <- beta
        gamma0 <- gamma

        #- calculate prior
        prior <- immer_latent_regression_prior_normal( mu=mu, sigma=sigma,
                    theta=theta )
        #-- calculate posterior
        res <- immer_latent_regression_posterior( like=like, prior=prior,
                    weights=weights )
        post <- res$post
        loglike <- res$loglike
        deviance <- -2*loglike

        #-- M-step regression
        theta_bar <- rowSums( post * thetaM )
        beta <- XtX %*% crossprod(Xw, theta_bar )
        mu <- as.vector( X %*% beta )
        resid <- thetaM - matrix( mu, nrow=N, ncol=TP)
        rp <- resid^2 * post
        for (gg in 1:G){
            ind_gg <- ind_group[[gg]]
            N_gg <- sum(weights[ind_gg])
            gamma[gg] <- sqrt( sum( rp[ind_gg,]*weights[ind_gg] ) / N_gg )
        }
        sigma <- gamma[ group ]
        parm_change <- max( abs( c( beta-beta0, gamma - gamma0) ))
        utils::flush.console()
        iter <- iter + 1
        if (iter > maxit){ iterate <- FALSE }
        if (parm_change < conv){
            iterate <- FALSE
        }
        if (verbose){
            cat("* Iteration", iter-1,    "\n" )
            v1 <- paste0("  Deviance=", round( deviance, 6 ),
                "  | Parameter change=", round(parm_change, 8) )
            cat(v1, "\n")
            utils::flush.console()
        }
    }
    #-- coefficients
    beta <- as.vector(beta)
    NB <- length(beta)
    names(beta) <- paste0("beta_",colnames(X))
    gamma <- as.vector(gamma)
    names(gamma) <- paste0("gamma_Group", 1:G)

    #-- standard errors
    h <- 1e-4
    pars <- c(beta, gamma)
    infomat <- immer_latent_regression_vcov_xpd( X=X, group=group, G=G, pars=pars,
                    theta=theta, weights=weights, like=like, h=h )
    V <- MASS::ginv(infomat)
    rownames(V) <- colnames(V) <- names(pars)
    se <- sqrt( diag(V) )
    beta_stat <- data.frame( parm=names(beta), est=pars[1:NB],
                    se=se[1:NB] )
    gamma_stat <- data.frame( parm=names(gamma), est=pars[NB+1:G],
                    se=se[NB+1:G] )

    #-- information criteria
    ic <- list(dev=deviance, n=N)
    ic$np <- length(beta) + length(gamma)
    ic$G <- G
    ic <- immer_IC_calc(ic=ic)

    #--- output
    time$end <- Sys.time()
    time$diff <- time$end - time$start
    description <- "Unidimensional latent regression"
    res <- list( coef=pars, vcov=V, beta=beta, gamma=gamma,
                    beta_stat=beta_stat, gamma_stat=gamma_stat,
                    sigma=sigma, mu=mu, sigma=sigma, ic=ic,
                    loglike=loglike, deviance=deviance, N=N, G=G, group=group,
                    iter=iter-1, time=time, CALL=CALL,
                    description=description )
    class(res) <- "immer_latent_regression"
    return(res)
}
