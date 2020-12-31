## File Name: immer_cmml_prepare.R
## File Version: 0.32

immer_cmml_prepare <- function( dat, Q=NULL, weights=NULL)
{

    #-- data processing
    dat0 <- dat
    I <- ncol(dat)
    dat <- as.matrix(dat)
    items <- colnames(dat)
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
    dfr <- immer_cmml_proc_freq(dat=dat, dat_resp=dat_resp, K=K, weights=weights)
    colnames(dfr) <- c("item1","item2","cat1","cat2","freq","ll_index")
    dfr <- as.data.frame(dfr)
    sz <- rowsum( dfr$freq, dfr$ll_index )
    dfr$ntot <- sz[dfr$ll_index]
    item_table <- dfr[, 1:4]
    item_pairs <- item_table[ ( item_table$cat1==0 ) & ( item_table$cat2==0 ), ]
    item_pairs <- item_pairs[, c("item1", "item2") ]

    #----------------------------
    #-- design matrices
    cpml_design <- list()

    #*** rho
    rho <- matrix( 0.5, nrow=I, ncol=I)
    diag(rho) <- 1
    rho <- immer_matrix_add_names(x=rho, row=items, col=items)
    rho_index <- immer_cmml_proc_generate_rho_index(I=I)
    rho_index <- immer_matrix_add_names(x=rho_index, row=items, col=items)
    N_rho <- max( rho_index )
    cpml_design$rho <- rho
    cpml_design$rho_index <- rho_index
    cpml_design$N_rho <- N_rho

    #*** tau
    res <- immer_cmml_proc_generate_tau(maxK=maxK, K=K)
    thresh <- paste0("thresh",0:(K+1))
    tau <- immer_matrix_add_names(x=res$tau, row=items, col=thresh)
    tau_index <- immer_matrix_add_names(x=res$tau_index, row=items, col=thresh)
    #***
    #tau_index <- tau_index - 1
    #***
    cpml_design$tau <- tau
    cpml_design$tau_index <- tau_index
    cpml_design$N_tau <- max( tau_index )

    #*** LAM
    if ( is.null(Q) ){
        Q <- matrix( 1, nrow=I, ncol=1)
        colnames(Q) <- "F"
    }
    Q <- immer_matrix_add_names(x=Q, row=items)
    Q_col <- colnames(Q)
    D <- ncol(Q)
    cpml_design$Q <- Q
    res <- immer_cmml_proc_generate_LAM(Q=Q)
    LAM <- immer_matrix_add_names(x=res$LAM, row=items, col=Q_col)
    LAM_index <- immer_matrix_add_names(x=res$LAM_index, row=items, col=Q_col )
    LAM_init <- 0*LAM
    W_LAM <- res$W_LAM
    LAM_basispar <- res$LAM_basispar
    LAM_basispar <- immer_vector_add_names(LAM_basispar, "LAM_b")
    cpml_design$LAM <- LAM
    cpml_design$LAM_index <- LAM_index
    cpml_design$LAM_init <- LAM_init
    cpml_design$W_LAM <- W_LAM
    cpml_design$LAM_basispar <- LAM_basispar


    #*** GAM
    GAM <- tau
    GAM_index <- tau_index - 1
    GAM_init <- tau * ( abs(tau)> 100 )
    N_GAM <- max(GAM_index)
    W_GAM <- diag(N_GAM)
    GAM_basispar <- GAM[ GAM_index > 0 ]
    GAM_basispar <- immer_vector_add_names(GAM_basispar, "GAM_b")

    cpml_design$GAM <- GAM
    cpml_design$GAM_index <- GAM_index
    cpml_design$GAM_init <- GAM_init
    cpml_design$W_GAM <- W_GAM
    cpml_design$GAM_basispar <- GAM_basispar

    #*** PHI
    res <- immer_cmml_proc_generate_PHI(D=D, use_diag=TRUE)
    PHI <- immer_matrix_add_names(x=res$PHI, row=Q_col, col=Q_col)
    PHI_index <- immer_matrix_add_names(x=res$PHI_index, row=Q_col, col=Q_col)
    PHI_init <- immer_matrix_add_names(x=res$PHI_init, row=Q_col, col=Q_col)
    W_PHI <- res$W_PHI
    PHI_basispar <- immer_vector_add_names(res$PHI_basispar, "PHI_b")

    cpml_design$PHI <- PHI
    cpml_design$PHI_index <- PHI_index
    cpml_design$PHI_init <- PHI_init
    cpml_design$W_PHI <- W_PHI
    cpml_design$PHI_basispar <- PHI_basispar

    #*** PSI
    res <- immer_cmml_proc_generate_PHI(D=I, use_diag=FALSE)
    PSI <- immer_matrix_add_names(x=res$PHI, row=items, col=items)
    PSI_index <- immer_matrix_add_names(x=res$PHI_index, row=items, col=items)
    PSI_init <- immer_matrix_add_names(x=res$PHI_init, row=items, col=items)
    W_PSI <- res$W_PHI
    PSI_basispar <- immer_vector_add_names(res$PHI_basispar, "PSI_b")

    cpml_design$PSI <- PSI
    cpml_design$PSI_index <- PSI_index
    cpml_design$PSI_init <- PSI_init
    cpml_design$W_PSI <- W_PSI
    cpml_design$PSI_basispar <- PSI_basispar

    attr(dfr, "N") <- N

    #--- output
    res <- list( dat=dat, dat_resp=dat_resp, maxK=maxK, K=K, weights=weights,
                    suff_stat=dfr, item_table=as.matrix(item_table),
                    item_pairs=as.matrix(item_pairs), cpml_design=cpml_design )
    return(res)
}
