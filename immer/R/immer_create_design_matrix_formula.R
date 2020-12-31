## File Name: immer_create_design_matrix_formula.R
## File Version: 0.16


immer_create_design_matrix_formula <- function( itemtable, formulaA )
{
    #-- process table
    has_rater <- c("rater") %in% colnames(itemtable)
    itemtable <- as.data.frame(itemtable)
    itemtable$index <- seq(1, nrow(itemtable))
    itemtable <- immer_convert_factor_data_frame(x=itemtable, variable="item" )
    itemtable <- immer_convert_factor_data_frame(x=itemtable, variable="rater" )
    if ( ! ( "item_name" %in% colnames(itemtable) ) ){
        itemtable[,"item_name"] <- paste0( "item", itemtable$item )
    }
    if ( ! ( "rater_name" %in% colnames(itemtable) ) ){
        itemtable[,"rater_name"] <- paste0( "rater", itemtable$rater )
    }
    if (has_rater){
        itemtable$pseudoitem_name <- paste0( itemtable$item_name, "-", itemtable$rater_name )
    } else {
        itemtable$pseudoitem_name <- itemtable$item_name
    }

    #-- include steps
    itemtable2 <- NULL
    K <- max(itemtable$maxK)
    for (kk in 1:K){
        itemtable1 <- itemtable
        itemtable1$step <- kk
        itemtable2 <- rbind(itemtable2, itemtable1)
    }
    itemtable2 <- itemtable2[ order(itemtable2$index), ]
    itemtable2$step_orig <- itemtable2$step
    itemtable2$step_num <- itemtable2$step
    itemtable2$step <- itemtable2$step + K - itemtable2$maxK
    itemtable2 <- itemtable2[ itemtable2$step_orig <=itemtable2$maxK, ]
    itemtable2 <- immer_convert_factor_data_frame(x=itemtable2, variable="step" )
    itemtable2$item_num <- as.numeric( paste(itemtable2$item) )
    if (has_rater){
        rownames(itemtable2) <- paste0( itemtable2$item_name, "-", itemtable2$rater_name, "-",
                                    "step", itemtable2$step_orig)
        itemtable2$rater_num <- as.numeric( paste(itemtable2$rater) )
    }
    has_step <- ( K >=2 )

    #-- modify formula
    tA <- stats::terms.formula(x=formulaA)
    terms_A <- unique( c( "-1", attr(tA, "term.labels") ) )
    formulaA <- as.formula( paste0( "~ ", paste0( terms_A, collapse=" + " ) ) )

    #-- apply formula
    contrast_list <- list(item="contr.sum", rater="contr.sum", step="contr.sum")
    taf <- rownames( attr(tA,"factors") )
    if ( ( ! has_rater ) |  ( ! ( "rater" %in% taf ) ) ){
        contrast_list$rater <- NULL
    }
    if ( ( ! has_step ) |  ( ! ( "step" %in% taf ) ) ){
        contrast_list$step <- NULL
    }
    des <- stats::model.matrix( object=formulaA, data=itemtable2, contrasts=contrast_list )
    colnames(des) <- gsub( "item", "__item", colnames(des) )

    #--- multiply by step parameters
    NP <- ncol(des)
    ind <- setdiff( 1:NP, grep( "step", colnames(des) ) )
    des[,ind] <- itemtable2$step_orig * des[,ind ]

    #--- remove columns in design matrix
    i1 <- itemtable[ itemtable$maxK < K, ]
    i1 <- i1[ ! duplicated(i1$item), ]
    negsum <- ( colMeans( des <=0 )==1 )
    ind_neg <- which(negsum)
    cn_des <- colnames(des)
    NI <- nrow(i1)
    if (NI > 0){
        for (ii in 1:NI){
            item_ii <- paste0( "__item", i1$item[ii] )
            ind <- intersect( grep( ":step", cn_des ), grep( paste0( item_ii, ":"), cn_des ))
            ind <- intersect( ind, ind_neg)
            if ( length(ind) > 0 ){
                des <- des[, - ind ]
            }
        }
    }

    #--- rename parameters in design matrix
    i1 <- itemtable2[ itemtable2$maxK < K, ]
    if (has_rater){
        i1 <- i1[ i1$rater==1, ]
    }
    items <- unique( paste(i1$item_name) )
    NI <- length(items)
    if (NI>0){
        for (ii in 1:NI){
            item_ii <- items[ii]
            i1_ii <- i1[ i1$item_name %in% item_ii,, drop=FALSE]
            NX <- nrow(i1_ii)
            if (NX>0){
                item_ii0 <- paste0("__item", i1_ii$item[1] )
                for (xx in 1:NX){
                    v1 <- paste0( item_ii0, ":step", i1_ii$step[xx] )
                    ind0 <- grep( v1, colnames(des) )
                    colnames(des) <- gsub( paste0( item_ii0, ":step", i1_ii$step[xx] ),
                                paste0( item_ii0, ":step", i1_ii$step_orig[xx] ), colnames(des) )
                }
            }
        }
    }

    #--- reorder column names
    NP <- ncol(des)
    dfr_cndes <- data.frame("index"=1:NP, "cn"=colnames(des) )
    ind <- intersect( grep( "item", colnames(des) ), grep( "step", colnames(des) ) )
    if ( length(ind) > 0 ){
        dfr_cndes$item_step <- 0
        dfr_cndes$item_step[ ind ] <- 1
        dfr_cndes$item_extract[ ind ] <- gsub("__item", "", paste0( dfr_cndes[ ind, "cn" ] ) )
        dfr_cndes$item_extract[ ind ] <-  as.numeric(unlist( lapply( strsplit( paste0( dfr_cndes[ ind, "item_extract" ] ),
                                        split=":"), FUN=function(ll){ ll[1] } ) ) )
        dfr_ind <- dfr_cndes[ind,]
        dfr_cndes[ind, ] <- dfr_ind[ order( dfr_ind$item_extract), ]
        des <- des[, dfr_cndes$index ]
        colnames(des) <- dfr_cndes$cn
    }

    #--- rename parameter names
    itemtable_ren <- itemtable[ ! duplicated( itemtable$item ),, drop=FALSE]
    NI <- nrow( itemtable_ren )
    for (ii in 1:NI){
        item_ii <- paste(itemtable[ ii, "item"])
        ind <- which( colnames(des)==paste0( "__item", item_ii ) )
        if (length(ind) > 0 ){
            colnames(des)[ ind ] <- itemtable[ ii, "item_name"]
        }
        ind <- grep( paste0( "__item", item_ii, ":"), colnames(des))
        if (length(ind) > 0 ){
            colnames(des)[ ind ] <- gsub( paste0( "__item", item_ii, ":"),
                                        paste0( itemtable[ii,"item_name"], ":"), colnames(des)[ind] )
        }
    }

    #--- create matrix A
    I <- nrow(itemtable)
    A <- array( 0, dim=c(I,K+1,NP) )
    dimnames(A) <- list( itemtable$pseudoitem_name, paste0("cat", 0:K), colnames(des) )
    for (kk in 1:K){
        ind_kk <- which(itemtable2$step_orig==kk)
        itemtable2_kk <- itemtable2[ ind_kk, ]
        des_kk <- des[ ind_kk, ]
        for (pp in 1:NP){
            A[ itemtable2_kk$index, kk+1, pp ] <- des_kk[,pp]
        }
    }

    #--- output
    res <- list( itemtable=itemtable, itemtable2=itemtable2,
                    formulaA=formulaA, has_rater=has_rater, K=K, has_step=has_step, des=des,
                    tA=tA, A=A )
    return(res)
}
