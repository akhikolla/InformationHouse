## File Name: immer_jml.R
## File Version: 0.9636


immer_jml <- function(dat, A=NULL, maxK=NULL, center_theta=TRUE, b_fixed=NULL, irtmodel="PCM",
                pid=NULL, rater=NULL, eps=.3, est_method="eps_adj",
                maxiter=1000, conv=1E-5, max_incr=3, maxiter_update=10, maxiter_line_search=6,
                conv_update=1E-5,
                verbose=TRUE, use_Rcpp=TRUE, shortcut=TRUE )
{
    time <- list( start=Sys.time() )
    CALL <- match.call()

    #-- process rating data
    res <- immer_proc_data( dat=dat, pid=pid, rater=rater, weights=NULL, maxK=maxK)
    dat <- res$dat2.NA
    pid <- res$pid
    maxK <- res$maxK
    K <- res$K
    pseudoitems_design <- res$pseudoitems_design

    #-- shortcut for analyzing the dataset
    res <- immer_jml_proc_shortcut( dat=dat, pid=pid, shortcut=shortcut, weights=NULL)
    dat <- res$dat
    pid <- res$pid
    shortcut <- res$shortcut
    N <- res$N
    shortcut_index <- res$shortcut_index
    weights <- res$weights

    #-- sufficient statistics
    dat0 <- dat
    dat_resp <- 1 - is.na(dat)
    dat[ is.na(dat) ] <- 0
    I <- ncol(dat)
    I_adj <- mean( rowSums(dat_resp) )
    bc_adj_fac <- ( I_adj - 1 ) / I_adj

    #-- create design matrix if not provided
    res <- immer_create_design_matrix_A( maxK=maxK, A=A, b_fixed=b_fixed, irtmodel=irtmodel )
    A <- res$A
    b_fixed <- res$b_fixed
    irtmodel <- res$irtmodel

    #-- scoring
    maxK_M <- immer_matrix2(maxK, nrow=N)
    max_pers <- rowSums( maxK_M * dat_resp )
    sumscore_pers <- rowSums(dat)

    eps_adjust_pers <- FALSE
    eps_adjust_item <- FALSE
    if (est_method=="eps_adj"){
        eps_adjust_pers <- TRUE
        eps_adjust_item <- TRUE
    }

    person <- data.frame(index=1:N, sum_score=sumscore_pers, max_pers=max_pers )
    if (eps_adjust_pers){
        person$eps <- eps
        max_pers1 <- maxK_M * dat_resp
        person$N_adj <- rowSums( max_pers1 * (max_pers1 + 1 ) / 2 )
        person$eps_i <- person$eps / person$N_adj
        person$score_pers <- person$eps + ( person$max_pers - 2 * person$eps )/person$max_pers * person$sum_score
        person$score_pers_adj <- immer_score_person_adjusted( sum_score=person$sum_score,
                                max_pers=person$max_pers, eps=eps)

    } else {
        person$eps <- eps
        person$eps_i <- 0
        person$score_pers <- immer_score_person_adjusted( sum_score=person$sum_score,
                                max_pers=person$max_pers, eps=eps)
    }


    dat_score <- array( dat_resp * person$eps_i, dim=c(N,I,K+1) )
    dat_score2 <- dat_score
    for (ii in 1:I){
        if (maxK[ii] < K){
            dat_score[,ii, seq(maxK[ii]+2,K+1) ] <- 0
        }
        for (kk in 0:K){
            dat_add <- ( dat[,ii]==kk ) * ( 1 - person$eps_i * (maxK[ii]+1) ) * dat_resp[,ii]
            dat_score[,ii,kk+1] <- dat_score[,ii,kk+1] + dat_add
            dat_score2[,ii,kk+1] <- kk*dat_score[,ii,kk+1]
        }
    }

    #-- sufficient statistics
    score_pers <- rep(0,N)
    for (kk in 1:K){
        score_pers <- score_pers + kk * rowSums( dat_score[,,kk+1] )
        person$score_pers0 <- score_pers
    }
    score_pers <- person$score_pers

    score_items0 <- matrix( NA, nrow=I, ncol=K)
    score_items <- matrix( NA, nrow=I, ncol=K)
    for (kk in 1:K){
        score_items0[,kk] <- colSums( (dat==kk)*dat_resp*weights, na.rm=TRUE )
        score_items[,kk] <- colSums( dat_score[,,kk+1]*weights, na.rm=TRUE )
    }
    if ( ! eps_adjust_item){
        score_items <- score_items0
    }
    dimA <- dim(A)
    NX <- dim(A)[3]
    ItemScore <- rep(0,NX)
    for (xx in 1:NX){
        ItemScore[xx] <- sum( A[,,xx] * score_items )
    }

    #-- initial parameters
    theta <- stats::qlogis( ( score_pers + .5) / ( max_pers + 1) )
    theta <- immer_jml_center_theta( theta=theta, center_theta=center_theta )

    xsi <- rep(0,NX)
    b <- matrix(0, nrow=I, ncol=K)
    if (is.null(b_fixed)){
        b_fixed <- b
    }

    #-- process data
    N <- length(theta)
    I <- ncol(dat)

    iter <- 0
    iterate <- TRUE

    if (use_Rcpp){
        fct_item <- immer_jml_update_item_Rcpp
        fct_theta <- immer_jml_update_theta_Rcpp
    } else {
        fct_item <- immer_jml_update_item_R
        fct_theta <- immer_jml_update_theta_R
    }

    deviance <- Inf
    parm_minimal <- list( deviance=deviance )
    deviance.history <- rep(NA, maxiter)

    #--- begin algorithm
    while (iterate){

        xsi0 <- xsi
        theta0 <- theta
        deviance0 <- deviance

        #** update item parameters
        args_item <- list( score_items=score_items, I=I, K=K, b=b, A=A, xsi=xsi, theta=theta, N=N, dat_resp=dat_resp,
                    max_incr=max_incr, maxiter_update=maxiter_update, conv_update=conv_update,
                    b_fixed=b_fixed, ItemScore=ItemScore, shortcut_index=shortcut_index, weights=weights )
        res <- do.call( what=fct_item, args=args_item)
        b <- res$b
        xsi <- res$xsi
        xsi_der2 <- res$xsi_der2
        # probs <- res$probs

        #** update person parameters
        args_theta <- list( score_pers=score_pers, I=I, K=K, N=N, theta=theta, b=b, dat_resp=dat_resp,
                        maxiter_update=maxiter_update, conv_update=conv_update, center_theta=center_theta,
                        max_incr=max_incr, shortcut_index=shortcut_index )
        res <- do.call( what=fct_theta, args=args_theta)
        theta <- res$theta
        theta_der2 <- res$theta_der2
        probs <- res$probs

        #** calculate log-likelihood
        loglike <- immer_jml_calc_loglike( dat_score=dat_score, probs=probs, K=K, weights=weights)
        deviance <- -2*loglike

        if (verbose){
            cat("* Iteration ", iter+1,    "\n" )
        }

        #--- line search in case of non-decreasing deviance
        if (deviance > deviance0 ){
            immer_jml_print_progress_line_search( verbose=verbose, deviance=deviance, digits_deviance=6)
            lambda <- 1
            lambda_fac <- 2
            # lambda_fac <- 1.4
            xsi_old <- xsi0
            xsi_new <- xsi
            theta_old <- theta0
            theta_new <- theta
            iter_in <- 0
            maxiter_ls <- maxiter_line_search
            while( ( deviance > deviance0 ) & ( iter_in < maxiter_ls ) ){
                lambda <- lambda / lambda_fac
                xsi <- xsi_old + lambda * ( xsi_new - xsi_old)
                b <- immer_jml_calc_item_intercepts(A=A, xsi=xsi)
                args_theta$b <- b
                theta <- theta_old + lambda * ( theta_new - theta_old)
                args_theta$maxiter_update <- maxiter_update
                # args_theta$maxiter_update <- 1
                args_theta$theta <- theta
                # args_theta$theta <- theta_old
                args_theta$center_theta <- FALSE
                args_theta$max_incr <- 1
                res <- do.call( what=fct_theta, args=args_theta)
                theta <- res$theta
                theta_der2 <- res$theta_der2
                probs <- res$probs
                loglike <- immer_jml_calc_loglike( dat_score=dat_score, probs=probs, K=K, weights=weights)
                deviance <- -2*loglike
                iter_in <- iter_in + 1
                immer_jml_print_progress_line_search( verbose=verbose, deviance=deviance, digits_deviance=6)
            }
        }

        deviance.history[iter+1] <- deviance

        #--- collect result list
        if ( deviance < parm_minimal$deviance ){
            parm_minimal <- list( theta=theta, xsi=xsi, xsi_der2=xsi_der2,
                    probs=probs, iter_opt=iter + 1, b=b, deviance=deviance,
                    loglike=loglike)
        }

        iter <- iter + 1
        item_parm_change <- max( abs(xsi0-xsi))
        pers_parm_change <- max( abs(theta0-theta))

        if (iter >=maxiter ){ iterate <- FALSE }
        if ( (item_parm_change < conv) & (pers_parm_change < conv )){
            iterate <- FALSE
        }

        if (verbose){
            v2 <- paste0( " | Deviance change=", round( deviance0 - deviance, 6) )
            if ( iter==1){ v2 <- "" }
            v1 <- paste0("  Deviance=", round( deviance, 6 ), v2, "\n",
                        "    Item parm. change=", round(item_parm_change, 8) )
            v1 <- paste0( v1, " | Person parm. change=", round(pers_parm_change, 8) )
            cat(v1, "\n")
            utils::flush.console()
        }

    }
    #--------------

    #-- assign elements of parm_minimal
    iter_opt <- NULL
    envir <- environment()
    res <- immer_attach_list_elements(x=parm_minimal, envir=envir)

    if ( est_method !="jml_bc" ){
        bc_adj_fac <- 1
    }
    if ( est_method %in% c("jml_bc","eps_adj") ){
        xsi <- bc_adj_fac * xsi
        b <- immer_jml_calc_item_intercepts(A=A, xsi=xsi, b_fixed=b_fixed)
        args_theta$theta <- theta
        args_theta$b <- b
        if ( est_method=="eps_adj"){
            args_theta$score_pers <- person$score_pers_adj
        }
        res <- do.call( what=fct_theta, args=args_theta)
        theta <- res$theta
        theta_der2 <- res$theta_der2
        probs <- res$probs
    }

    #-- standard errors
    eps0 <- 1E-20
    # se_adj <- bc_adj_fac
    se_adj <- 1
    xsi_se <- se_adj * sqrt( abs(1 / xsi_der2 ) + eps0)
    theta_se <- sqrt( abs(1 / theta_der2 ) + eps0)

    #-- information criteria
    ic <- immer_jml_ic( loglike=loglike, N=N, center_theta=center_theta, xsi=xsi, I=I )

    #--- reorder data if shortcut is used
    if (shortcut){
        ind <- order( shortcut_index$orig )
        dat <- dat[ ind, ]
        dat_resp <- dat_resp[ ind, ]
        pid <- pid[ ind ]
        theta <- theta[ind]
        theta_se <- theta_se[ ind ]
        dat_score <- dat_score[ ind,, ]
        score_pers <- score_pers[ ind ]
    }

    #-- person parameters
    person$theta <- theta
    person$theta_se <- theta_se
    sd_obs <- stats::sd(theta)
    # average reliability
    rel <- max( 0, 1 - mean(theta_se^2) / sd_obs^2     )
    sd_lat <- sqrt( sd_obs^2 * rel )
    person_desc <- list( mean=mean(theta), sd_obs=sd_obs, mle_rel=rel,
                        sd_lat=sd_lat, min=min(theta), max=max(theta) )

    #-- item parameters
    rownames(b) <- colnames(dat)
    colnames(b) <- paste0("Cat", 1:K)
    item <- data.frame("item"=colnames(dat), "N"=colSums(dat_resp),
                    "M"=colSums( dat * dat_resp ) / colSums(dat_resp) )
    item$SD <- colSums( dat^2 * dat_resp ) / colSums(dat_resp)
    item$SD <- sqrt( item$SD - item$M^2 )
    item <- data.frame( item, b )
    names(xsi) <- dimnames(A)[[3]]
    xsi_dfr <- data.frame("par"=names(xsi), "est"=xsi, "se"=xsi_se )

    #--- output
    time$end <- Sys.time()
    time$diff <- time$end - time$start
    res <- list( b=b, item=item, theta=theta, theta_se=theta_se, xsi=xsi,
                xsi_se=xsi_se, xsi_dfr=xsi_dfr, probs=probs, person=person,
                person_desc=person_desc, shortcut_index=shortcut_index, pid=pid, dat=dat,
                dat_score=dat_score, dat_resp=dat_resp, score_pers=score_pers,
                score_items=score_items, eps=eps, eps_adjust_pers=eps_adjust_pers, eps=eps,
                eps_adjust_item=eps_adjust_item, est_method=est_method, A=A, b_fixed=b_fixed,
                bc_adj_fac=bc_adj_fac, irtmodel=irtmodel, iter=iter, iter_opt=iter_opt,
                loglike=loglike, deviance=deviance, ic=ic, deviance.history=deviance.history,
                pseudoitems_design=pseudoitems_design, CALL=CALL, time=time )
    res$description <- "Function immer_jml() | Joint maximum likelihood estimation"
    class(res) <- "immer_jml"
    return(res)
}
