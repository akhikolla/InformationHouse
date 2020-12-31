## File Name: immer_cml.R
## File Version: 0.842

########################################################
# CML function in immer package
immer_cml <- function( dat, weights=NULL, W=NULL, b_const=NULL,
        par_init=NULL, a=NULL, irtmodel=NULL, normalization="first",
        nullcats="zeroprob", diff=FALSE, use_rcpp=FALSE, ... )
{
    s1 <- Sys.time()
    CALL <- match.call()
    #---------------------
    # data preparation
    res <- lpcm_data_prep( dat=dat, weights=weights, a=a )
    N <- res$N
    NP <- res$NP
    I <- res$I
    dat <- res$dat
    weights <- res$weights
    patt <- res$patt
    suffstat <- res$suffstat
    splitvec <- res$splitvec
    item_index <- res$item_index
    parm_index <- res$parm_index
    pars_info <- res$pars_info
    maxscore <- res$maxscore
    maxK <- res$maxK
    K <- max(maxK)
    score <- res$score
    used_persons <- res$used_persons
    score_freq <- res$score_freq
    a <- res$a

    #--------------------------
    # generate W matrix and b_const if not provided
    res0 <- lpcm_generate_design( pars_info=pars_info, irtmodel=irtmodel, W=W, b_const=b_const,
                    normalization=normalization, I=I, maxK=maxK, nullcats=nullcats )
    W <- res0$W
    b_const <- res0$b_const
    irtmodel <- res0$irtmodel

    #--------------------------
    # initial parameters
    if ( is.null(par_init) ){
        par_init <- lpcm_inits( dat=dat, weights=weights, maxK=maxK, b_const=b_const, W=W, irtmodel=irtmodel,
                            normalization=normalization )
    }

    #*** compute sufficient statistics for missing data patterns
    gr1_list <- list()
    for (pp in 1:NP){
        W1 <- W[ parm_index[[pp]],, drop=FALSE ]
        gr1_list[[pp]] <- suffstat[[pp]] %*% W1
    }

    #---------------------------
    # functions passed to optimization function
    ## objective function: conditional log-likelihood

    if (use_rcpp){
        parm_index0 <- parm_index
        splitvec_len <- list()
        for (pp in 1:NP){
            parm_index0[[pp]] <- parm_index0[[pp]] - 1
            splitvec_len[[pp]] <- as.vector( table(splitvec[[pp]]) )
            gr1_list[[pp]] <- as.vector( gr1_list[[pp]])
        }
        W_logical <- W !=0
        cloglik <- function(par){
            esf_par0 <- W %*% par + b_const
            res <- immer_cml_cloglik_helper( esf_par0=esf_par0,
                        parm_index=parm_index0, splitvec_len=splitvec_len, suffstat=suffstat,
                        score_freq=score_freq, diff=diff, NP=NP )
            return(-res)
        }
    } else {
        cloglik <- function (par){
                esf_par0 <- W %*% par + b_const
                cll <- sum( unlist( sapply( 1:NP, FUN=function(pp){
                    esf_par <- esf_par0[ parm_index[[pp]], 1 ]
                    b <- esf_par
                    esf_par <- split(esf_par, splitvec[[pp]] )
                    esf <- psychotools::elementary_symmetric_functions(
                                        par=esf_par, order=0, diff=diff)[[1]]
                    res <- - sum(suffstat[[pp]] * b ) - sum(score_freq[[pp]] * log(esf))
                    return(res)
                } )   ) )
            return(-cll)
        }
    }

    if (use_rcpp){
        ## analytical gradient in Rcpp
        agrad <- function(par){
            esf_par0 <- W %*% par + b_const
            res <- immer_cml_agrad_helper( esf_par0=esf_par0,
                            parm_index=parm_index0, splitvec_len=splitvec_len, suffstat=suffstat,
                            score_freq=score_freq, diff=diff, NP=NP, W=W, W_logical=W_logical,
                            gr1_list=gr1_list )
            return(res)
        }
    } else {
        ## analytical gradient in R
        agrad <- function(par){
            esf_par0 <- W %*% par + b_const
            h1 <- sapply( 1:NP, FUN=function(pp){
                esf_par <- esf_par0[ parm_index[[pp]], 1 ]
                esf_par <- split(esf_par, splitvec[[pp]] )
                esf <- psychotools::elementary_symmetric_functions(par=esf_par,
                            order=1, diff=diff)
                gamma0 <- esf[[1]]
                gamma1 <- esf[[2]]
                W1 <- W[ parm_index[[pp]],, drop=FALSE ]
                gr1 <- gr1_list[[pp]]
                gr2 <- - colSums( ( score_freq[[pp]] * (gamma1 / gamma0))  %*% W1 )
                gr <- gr1 + gr2
                gr <- gr[1,]
                return(gr)
                } )
            if (NP>1){
                h1 <- rowSums(h1)
            }
            return(h1)
        }
    }

    #-----------------
    # optim
    opt <- stats::optim(par=par_init, fn=cloglik, gr=agrad, method="BFGS",
                        hessian=TRUE, ... )

    #-----------------
    # summaries
    par <- opt$par
    b <- W %*% par + b_const
    b <- b[,1]

    item <- matrix( NA, nrow=I, ncol=K )
    item[ cbind( pars_info$itemid, pars_info$cat ) ] <- b
    colnames(item) <- paste0("Cat", 1:K)
    item[ item==99 ] <- NA
    # compute item difficulty
    itemdiff <- item[ cbind(1:I,maxK ) ] / maxK
    M <- colSums( weights * dat, na.rm=TRUE ) /
                    colSums( weights * ( 1 - is.na(dat) ), na.rm=TRUE )

    item <- data.frame( item=colnames(dat),
                    sumwgt=colSums( weights * ( 1 - is.na(dat) ), na.rm=TRUE ),
                    M=M, a=a,  itemdiff=itemdiff, item )
    rownames(item) <- NULL
    vcov1 <- solve(opt$hessian)
    par_summary <- data.frame( par=colnames(W), est=par, se=sqrt(diag(vcov1)) )

    #------ output
    s2 <- Sys.time()
    time <- list( start=s1, end=s2 )
    res <- list( item=item, b=b, coefficients=par, vcov=vcov1,
                par_summary=par_summary, loglike=-opt$value, deviance=2*opt$value,
                result_optim=opt, W=W, b_const=b_const, par_init=par_init, suffstat=suffstat,
                score_freq=score_freq, dat=dat, used_persons=used_persons, NP=NP, N=N, I=I,
                maxK=maxK, K=K, npars=length(par), pars_info=pars_info, parm_index=parm_index,
                item_index=item_index, score=score, time=time, CALL=CALL, use_rcpp=use_rcpp,
                description='Conditional Maximum Likelihood Estimation' )
    class(res) <- "immer_cml"
    return(res)
}
########################################################
