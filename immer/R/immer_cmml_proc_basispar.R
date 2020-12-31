## File Name: immer_cmml_proc_basispar.R
## File Version: 0.10


immer_cmml_proc_basispar <- function(LAM_basispar, GAM_basispar, PHI_basispar, PSI_basispar,
    LAM_index, GAM_index, PHI_index, PSI_index )
{
    basispar <- c(LAM_basispar, GAM_basispar, PHI_basispar, PSI_basispar )
    NB <- length(basispar)
    basispar_design <- matrix( 0, nrow=NB, ncol=4)
    colnames(basispar_design) <- c("LAM", "GAM", "PHI", "PSI")
    Ntemp <- 0
    vv <- 1

    #- LAM
    N_LAM <- length(LAM_basispar)
    if (N_LAM>0){
        ind <- Ntemp + 1:N_LAM
        basispar_design[ ind,vv] <- ind
        Ntemp <- Ntemp + N_LAM
    }
    vv <- vv + 1

    #- GAM
    N_GAM <-  length(GAM_basispar)
    if (N_GAM>0){
        ind <- Ntemp + 1:N_GAM
        basispar_design[ ind,vv] <- ind
        Ntemp <- Ntemp + N_GAM
    }
    vv <- vv + 1

    #- PHI
    N_PHI <- length(PHI_basispar)
    if (N_PHI>0){
        ind <- Ntemp + 1:N_PHI
        basispar_design[ ind,vv] <- ind
        Ntemp <- Ntemp + N_PHI
    }
    vv <- vv + 1

    #- PSI
    N_PSI <- length(PSI_basispar)
    if (N_PSI>0){
        ind <- Ntemp + 1:N_PSI
        basispar_design[ ind,vv] <- ind
        Ntemp <- Ntemp + N_PSI
    }
    vv <- vv + 1

    basispar_length <- colSums( basispar_design > 0 )

    #--- output
    res <- list( basispar=basispar, basispar_design=basispar_design, N_LAM=N_LAM, N_GAM=N_GAM,
                    N_PHI=N_PHI, N_PSI=N_PSI, basispar_length=basispar_length)
    return(res)
}

