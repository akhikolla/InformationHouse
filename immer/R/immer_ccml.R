## File Name: immer_ccml.R
## File Version: 0.17

immer_ccml <- function( dat, weights=NULL, irtmodel="PCM", A=NULL, b_fixed=NULL, control=NULL )
{
    time <- list( start=Sys.time() )
    CALL <- match.call()

    #-- data processing
    dat0 <- dat
    I <- ncol(dat)
    dat <- as.matrix(dat)
    dat_resp <- 1 - is.na(dat)
    dat[ is.na(dat) ] <- 0
    maxK <- apply( dat, 2, max )
    K <- max(maxK)
    N <- nrow(dat)
    if (is.null(weights)){
        weights <- rep(1,N)
    }
    W <- sum(weights)

    #-- count frequencies
    dfr <- immer_ccml_proc_freq( dat=dat, dat_resp=dat_resp, K=K, weights=weights )
    dfr <- as.data.frame(dfr)
    colnames(dfr) <- c("ll_index", "item1", "item2", "score", "cat1", "cat2", "n")
    dfr <- dfr[ dfr$ll_index > 0, ]
    sz <- rowsum( dfr$n, dfr$ll_index )
    dfr$ntot <- sz[dfr$ll_index]

    #-- create design matrix if not provided
    #... INCLUDE IT LATER

    #-- initial values xsi
    xsi <- rep( 0, dim(A)[3] )
    names(xsi) <- dimnames(A)[[3]]
    if ( is.null(b_fixed) ){
        b_fixed <- matrix(0, nrow=I, ncol=K+1)
    }

    #-- preparation optimization
    ll_index1 <- dfr$ll_index - 1
    par0 <- as.vector(xsi)
    A_ <- as.vector(A)
    item10 <- dfr$item1 - 1
    item20 <- dfr$item2 - 1
    cat1 <- dfr$cat1
    cat2 <- dfr$cat2
    max_ll_index <- max(dfr$ll_index)
    n <- dfr$n
    ntot <- dfr$ntot

    #- define optimization function
    opt_fct <- function(par){
        immer_ccml_opt_function_par( b_fixed=b_fixed, A_=A_, par=par, ll_index1=ll_index1,
                item10=item10, item20=item20, cat1=cat1, cat2=cat2, n=n, ntot=ntot,
                max_ll_index=max_ll_index )
    }
    #- define gradient
    grad_fct <- function(par){
        immer_ccml_gradient_par( b_fixed=b_fixed, A_=A_, par=par, ll_index1=ll_index1, item10=item10, item20=item20,
                    cat1=cat1, cat2=cat2, n=n, ntot=ntot, max_ll_index=max_ll_index )
    }

    #--- optimization
    res <- nlminb_result <- stats::nlminb( start=par0, objective=opt_fct, gradient=grad_fct, control=control )
    par <- coef <- res$par
    names(par) <- names(coef) <- names(xsi)
    objective <- res$objective

    #-- calculate intercept matrix
    b <- immer_ccml_calc_item_intercepts( b_fixed=b_fixed, A_=A_, par=par )

    #-- calculate standard errors
    res <- immer_ccml_se( b_fixed=b_fixed, A_=A_, par=par, ll_index1=ll_index1, item10=item10, item20=item20,
                cat1=cat1, cat2=cat2, n=n, ntot=ntot, max_ll_index=max_ll_index, h=1e-4 )
    J <- res$xpd_mat
    H <- res$obs_mat
    colnames(J) <- rownames(J) <- names(xsi)
    colnames(H) <- rownames(H) <- names(xsi)
    J1 <- MASS::ginv(J)
    V <- H %*% J1 %*% H
    se <- sqrt(diag(V))

    #-- information criteria
    ic <- list( objective=objective, np=length(xsi), N=N, I=I )
    H1 <- MASS::ginv(H)
    tr_HJ <- immer_trace( J %*% H1 )
    ic$dev <- 2*ic$objective
    ic$CLAIC <- ic$dev + 2*tr_HJ
    ic$CLBIC <- ic$dev + log(ic$N)*tr_HJ
    ic$R <- NA
    ic$ND <- NA

    #-- output xsi parameters
    xsi_out <- data.frame( est=coef, se=se )
    rownames(xsi_out) <- names(xsi)

    #-- output item parameters
    item <- b[,-1,drop=FALSE]
    colnames(item) <- paste0("Cat", 1:K)
    M <- colSums( dat * dat_resp * weights ) / colSums( dat_resp * weights )
    item <- data.frame( item=colnames(dat), M=M, item )

    #--- output
    time$end <- Sys.time()
    time$diff <- time$end - time$start
    description <- "Composite Conditional Maximum Likelihood Estimation"
    res <- list( coef=coef, vcov=V, se=se, b=b, objective=objective, nlminb_result=nlminb_result, A=A,
                    weights=weights, b_fixed=b_fixed, K=K, maxK=maxK, N=N, W=W, ic=ic, suff_stat=dfr,
                    xsi_out=xsi_out, item=item, time=time,
                    CALL=CALL, description=description)
    class(res) <- "immer_ccml"
    return(res)
}
