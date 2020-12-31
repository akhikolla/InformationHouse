## File Name: mnlfa_penalty_values_item.R
## File Version: 0.15


mnlfa_penalty_values_item <- function(parms, parms_iterations,
    parms_regular_types, parms_regular_lam, N_item, center_group_parms)
{
    NH <- length(parms_iterations)
    NP <- length(parms)
    parm_names <- names(parms)
    parms_values <- mnlfa_create_vector_with_names(vec=parm_names)
    parms_regularized <- parms_values
    parms_estimated <- parms_values
    for (hh in seq_len(NH) ){
        parms_indices <- parms_iterations[[hh]]
        regular_type <- parms_regular_types[parms_indices][1]
        regular_lam <- parms_regular_lam[parms_indices][1]
        res <- mnlfa_penalty_values( parms=parms, parms_indices=parms_indices,
                    regular_type=regular_type, regular_lam=regular_lam, N_item=N_item,
                    center_group_parms=center_group_parms)
        parms_values[ parms_indices ] <- res$val
        parms_regularized[ parms_indices ] <- res$regularized
        parms_estimated[ parms_indices ] <- res$estimated
    }

    #--- output
    res <- list(parms_values=parms_values, parms_regularized=parms_regularized,
                parms_estimated=parms_estimated)
    return(res)
}
