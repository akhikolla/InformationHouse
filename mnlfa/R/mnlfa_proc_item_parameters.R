## File Name: mnlfa_proc_item_parameters.R
## File Version: 0.418


mnlfa_proc_item_parameters <- function(dat, formula_int, formula_slo, item_type,
    parms_regular_types, parms_regular_lam, regular_type, regular_lam, parms_iterations)
{
    I <- length(item_type)
    items <- names(item_type)

    #* lists of formulas
    formula_int <- mnlfa_proc_item_parameters_formula_list(
                            formula_parm=formula_int, items=items)
    formula_slo <- mnlfa_proc_item_parameters_formula_list(
                            formula_parm=formula_slo, items=items)

    #** parameter lists
    parm_list <- list()
    parm_Xdes <- list()
    parm_index <- list()
    create_type <- FALSE
    if (is.null(parms_regular_types)){
        parms_regular_types    <- list()
        create_type <- TRUE
    }
    create_lam <- FALSE
    if (is.null(parms_regular_lam)){
        parms_regular_lam <- list()
        create_lam <- TRUE
    }
    create_iter <- FALSE
    if (is.null(parms_iterations)){
        parms_iterations <- list()
        create_iter <- TRUE
    }

    for (ii in 1:I){
        item_type_ii <- item_type[ii]
        item_ii <- names(item_type)[ii]
        #* 2PL
        if (item_type_ii %in% c("1PL","2PL") ){
            # design matrices
            Xdes_int <- stats::model.matrix( object=formula_int[[ii]], data=dat)
            Xdes_slo <- stats::model.matrix( object=formula_slo[[ii]], data=dat)
            des1 <- list(Xdes_int=Xdes_int, Xdes_slo=Xdes_slo)
            #- inits item parameters
            nc_int <- ncol(Xdes_int)
            offset_int <- FALSE
            if (nc_int > 0){
                g1 <- rep(0,ncol(Xdes_int) )
                names(g1) <- paste0(item_ii, "_b_", colnames(Xdes_int) )
            } else {
                g1 <- c(0)
                names(g1) <- paste0(item_ii, "_b_offset" )
                des1$Xdes_int <- matrix(0, nrow=nrow(Xdes_int), ncol=1)
                colnames(des1$Xdes_int) <- paste0(item_ii, "_b_offset" )
                offset_int <- TRUE
            }
            v1 <- list( b=g1 )

            #** design matrix item slopes
            nc_slo <- ncol(Xdes_slo)
            offset_slo <- FALSE
            if (nc_slo > 0){
                g2 <- rep(0,ncol(Xdes_slo))
                g2[1] <- 1
                names(g2) <- paste0(item_ii, "_a_", colnames(Xdes_slo) )
            } else {
                g2 <- c(1)
                names(g2) <- paste0(item_ii, "_a_offset" )
                des1$Xdes_slo <- matrix(0, nrow=nrow(Xdes_slo), ncol=1)
                colnames(des1$Xdes_slo) <- paste0(item_ii, "_a_offset" )
                offset_slo <- TRUE
            }

            v1$a <- g2
            #- regularization parameters
            np_b <- length(v1$b)
            np_a <- length(v1$a)
            reg1 <- rep("none", np_b + np_a)
            reg2 <- rep(0, np_b + np_a)
            names(reg1) <- c( names(v1$b), names(v1$a) )
            names(reg2) <- names(reg1)
            reg1[ 1 + seq_len( np_b - 1 ) ] <- regular_type[1]
            reg1[ np_b + 1 + seq_len( np_a - 1 ) ] <- regular_type[2]
            reg2[ 1 + seq_len( np_b - 1 ) ] <- regular_lam[1]
            reg2[ np_b + 1 + seq_len( np_a - 1 ) ] <- regular_lam[2]
            #- parameter indices
            b_index <- seq_len(np_b)
            # if (offset_int){ b_index <- NULL; }
            a_index <- np_b + seq_len(np_a)
            h1 <- list(b_index=b_index, a_index=a_index)
        }
        parm_list[[ii]] <- v1
        parm_Xdes[[ii]] <- des1
        parm_index[[ii]] <- h1
        if (create_type){
            parms_regular_types[[ii]] <- reg1
        }
        if (create_lam){
            parms_regular_lam[[ii]] <- reg2
        }
        if (create_iter){
            if (item_type[[ii]]=="2PL"){
                parms_iterations[[ii]] <- as.list(1:(np_b+np_a))
                if (offset_int & ( ! offset_slo) ){
                    parms_iterations[[ii]] <- as.list(2:(np_b+np_a))
                }
                if ( (!offset_int) & ( offset_slo)){
                    parms_iterations[[ii]] <- as.list(1:(np_b))
                }
                if (offset_int & offset_slo){
                    parms_iterations[[ii]] <- as.list(NULL)
                }
            }
            if (item_type[[ii]]=="1PL"){
                parms_iterations[[ii]] <- as.list(1:(np_b))
                if (offset_int){
                    parms_iterations[[ii]] <- as.list(NULL)
                }
            }
        }

    }  ## end item ii

    names(parm_list) <- items
    names(parm_Xdes) <- items
    names(parms_regular_types) <- items
    names(parms_regular_lam) <- items
    names(parms_iterations) <- items
    names(parm_index) <- items

    #--- output
    res <- list( parm_list=parm_list, parm_Xdes=parm_Xdes, parms_regular_types=parms_regular_types,
                    parms_regular_lam=parms_regular_lam, parms_iterations=parms_iterations,
                    parm_index=parm_index)
    return(res)
}


