## File Name: immer_cmml.R
## File Version: 0.29


immer_cmml <- function( suff_stat, cpml_design, control=NULL )
{
    time <- list( start=Sys.time() )
    CALL <- match.call()

    #*** logarithmized parameters
    do_log_LAM <- do_log_GAM <- do_log_PHI <- do_log_PSI <- 0

    #*** epsilon value for preventing division by zero for probabilities
    eps <- 1E-6

    #**** attach arguments
    freq <- as.vector(suff_stat[,"freq"])
    item_table <- as.matrix(suff_stat[, c("item1","item2","cat1","cat2")])
    item_pairs <- as.matrix(suff_stat[ (suff_stat$cat1==0) & (suff_stat$cat2==0), c("item1", "item2") ])

    LAM_init <- cpml_design$LAM_init
    GAM_init <- cpml_design$GAM_init
    PHI_init <- cpml_design$PHI_init
    PSI_init <- cpml_design$PSI_init
    LAM_basispar <- cpml_design$LAM_basispar
    GAM_basispar <- cpml_design$GAM_basispar
    PHI_basispar <- cpml_design$PHI_basispar
    PSI_basispar <- cpml_design$PSI_basispar
    W_LAM <- cpml_design$W_LAM
    W_GAM <- cpml_design$W_GAM
    W_PHI <- cpml_design$W_PHI
    W_PSI <- cpml_design$W_PSI
    LAM_index <- cpml_design$LAM_index
    GAM_index <- cpml_design$GAM_index
    PHI_index <- cpml_design$PHI_index
    PSI_index <- cpml_design$PSI_index
    rho_index <- cpml_design$rho_index
    tau_index <- cpml_design$tau_index
    N_rho <- max(rho_index)
    N_tau <- max(tau_index)

    #-- process design matrices in case of no estimated parameters
    W_LAM <- immer_cmml_proc_parameter_index( W_par=W_LAM, par_index=LAM_index )
    W_GAM <- immer_cmml_proc_parameter_index( W_par=W_GAM, par_index=GAM_index )
    W_PHI <- immer_cmml_proc_parameter_index( W_par=W_PHI, par_index=PHI_index )
    W_PSI <- immer_cmml_proc_parameter_index( W_par=W_PSI, par_index=PSI_index )

    #**** create basis parameters
    res <- immer_cmml_proc_basispar( LAM_basispar=LAM_basispar, GAM_basispar=GAM_basispar,
                PHI_basispar=PHI_basispar, PSI_basispar=PSI_basispar, LAM_index=LAM_index,
                GAM_index=GAM_index, PHI_index=PHI_index, PSI_index=PSI_index )
    basispar <- res$basispar
    basispar_design <- res$basispar_design
    basispar_length <- res$basispar_length
    N_LAM <- res$N_LAM
    N_GAM <- res$N_GAM
    N_PHI <- res$N_PHI
    N_PSI <- res$N_PSI

    basispar_design1 <- basispar_design - 1
    LAM_index1 <- LAM_index - 1
    GAM_index1 <- GAM_index - 1
    PHI_index1 <- PHI_index - 1
    PSI_index1 <- PSI_index - 1
    rho_index1 <- rho_index - 1
    tau_index1 <- tau_index - 1

    #- define optimization function
    opt_fct <- function(par){
        res <- immer_cmml_basispar_to_probs( basispar=par, basispar_design=basispar_design1,
                basispar_length=basispar_length, W_LAM=W_LAM, LAM_init=LAM_init, LAM_index=LAM_index1,
                W_GAM=W_GAM, GAM_init=GAM_init, GAM_index=GAM_index1, W_PHI=W_PHI, PHI_init=PHI_init,
                PHI_index=PHI_index1, W_PSI=W_PSI, PSI_init=PSI_init, PSI_index=PSI_index1,
                N_LAM=N_LAM, N_GAM=N_GAM, N_PHI=N_PHI, N_PSI=N_PSI, item_pairs=item_pairs,
                item_table=item_table, tau_index=tau_index1, rho_index=rho_index1, N_rho=N_rho,
                N_tau=N_tau, do_log_LAM=do_log_LAM, do_log_GAM=do_log_GAM, do_log_PHI=do_log_PHI,
                do_log_PSI=do_log_PSI )
        probs <- res$probs + eps
        val <- - sum( freq * log( probs ) )
        return(val)
    }

    res <- immer_cmml_basispar_to_derivatives( basispar=basispar, basispar_design=basispar_design1,
        basispar_length=basispar_length, W_LAM=W_LAM, LAM_init=LAM_init, LAM_index=LAM_index1,
        W_GAM=W_GAM, GAM_init=GAM_init, GAM_index=GAM_index1, W_PHI=W_PHI, PHI_init=PHI_init,
        PHI_index=PHI_index1, W_PSI=W_PSI, PSI_init=PSI_init, PSI_index=PSI_index1,
        N_LAM=N_LAM, N_GAM=N_GAM, N_PHI=N_PHI, N_PSI=N_PSI, item_pairs=item_pairs,
        item_table=item_table, tau_index=tau_index1, rho_index=rho_index1, N_rho=N_rho,
        N_tau=N_tau, do_log_LAM=do_log_LAM, do_log_GAM=do_log_GAM, do_log_PHI=do_log_PHI,
        do_log_PSI=do_log_PSI )


    #**** rearrange objects for output
    cpml_design$W_LAM <- W_LAM
    cpml_design$W_GAM <- W_GAM
    cpml_design$W_PHI <- W_PHI
    cpml_design$W_PSI <- W_PSI

    #--- output
    time$end <- Sys.time()
    time$diff <- time$end - time$start
    description <- "Composite Marginal Maximum Likelihood Estimation"
    res <- list( suff_stat=suff_stat, cpml_design=cpml_design,
                    time=time, CALL=CALL, description=description)
    class(res) <- "immer_cmml"
    return(res)
}


