## File Name: immer_jml_facets_proc_data.R
## File Version: 0.27


immer_jml_facets_proc_data <- function(y, design, weights, max_equal, eps, bc,
    center )
{
    #** weights
    if ( is.null(weights) ){
        weights <- rep(1, length(y))
    }

    #** convert center to logical
    des_names <- colnames(design)
    if (is.null(center) ){
        center <- setdiff( des_names, "item")
    }

    DN <- length(des_names)
    center0 <- center
    center <- rep(FALSE, DN)
    names(center) <- des_names
    center[ center0 ] <- TRUE

    #** compute maximum per item and per observation
    des_item <- paste(design[,"item"])
    max_item <- stats::aggregate( y, list(des_item), max )
    colnames(max_item) <- c("item", "max")
    if (max_equal){
        max_item$max <- max(max_item$max)
    }
    maxcat <- max(max_item$max)
    y_max <- max_item[ match( des_item, max_item$item ), "max" ]
    ind <- ! is.na(y)
    y <- y[ind]
    y_max <- y_max[ind]
    design <- design[ind, ]
    N <- nrow(design)

    #!! include step if it is not contained in design

    #** rearrange design matrix
    first_names <- c("person", "item", "step")
    v1 <- intersect( des_names, first_names)
    v2 <- setdiff(des_names, v1)
    des_names <- c(v1, v2)

    design <- design[, des_names]
    #*** include labels for facets
    design_labels <- list()
    for (dd in des_names){
        des_dd <- paste(design[,dd])
        des_lab_dd <- paste(sort( unique( paste( des_dd) ) ))
        ND <- length(des_lab_dd)
        dfr_dd <- data.frame( lab=paste(des_lab_dd), id=1:ND)
        attr(dfr_dd, "N") <- ND
        design_labels[[ dd ]] <- dfr_dd
        design[,dd] <- match( des_dd, des_lab_dd)
    }
    data <- data.frame(design, y=y, max=y_max, weights=weights)

    #** compute sufficient statistics
    yw <- y*weights
    y_maxw <- y_max*weights
    suff_stat <- list()
    parm <- list()
    parm_sign <- list()

    des_names2 <- setdiff( des_names, "step")
    for (dd in des_names2){
        dfr_dd <- rowsum(yw, design[,dd])
        dfr_dd <- data.frame(id=as.numeric(rownames(dfr_dd)), score_raw=dfr_dd[,1] )
        dfr_dd$score_max <- rowsum(y_maxw, design[,dd])[,1]
        dfr_dd$score_extreme <- ( dfr_dd$score_raw==0 ) + ( dfr_dd$score_raw==dfr_dd$score_max )
        dfr_dd$score <- dfr_dd$score_raw
        parm_sign_dd <- -1
        if (dd=="person"){
            parm_sign_dd <- 1
            if (!bc){
                ind <- dfr_dd$score==0
                dfr_dd[ ind, "score" ] <- eps
                ind <- dfr_dd$score==dfr_dd$score_max
                dfr_dd[ ind, "score" ] <- dfr_dd$score_max[ind] - eps
            }
            if (bc){
                dfr_dd$score <- eps + (dfr_dd$score_max - 2*eps)/ dfr_dd$score_max * dfr_dd$score
            }

        }
        dfr_dd$init <- parm_sign_dd*log( dfr_dd$score / (dfr_dd$score_max - dfr_dd$score ) )
        if ( center[dd] ){
            dfr_dd$init <- immer_center( x=dfr_dd$init )
        }
        suff_stat[[dd]] <- dfr_dd
        parm[[dd]] <- dfr_dd$init
        parm_sign[[dd]] <- parm_sign_dd
    }

    #*** sufficient statistics step parameters
    dd <- "step"
    des_dd <- design[,dd]
    N_dd <- max(des_dd)
    dfr_dd <- matrix(NA, nrow=N_dd, ncol=maxcat)
    for (mm in 1:maxcat){
        dfr_dd[,mm] <- rowsum( ( y >=mm ) * weights, des_dd)[,1]
    }
    suff_stat[[dd]] <- dfr_dd
    parm[[dd]] <- 0*dfr_dd
    parm_sign[[dd]] <- -1

    des_names1 <- setdiff(des_names, "step")
    is_step <- TRUE
    max_step <- stats::aggregate( y, list(design[,"step"]), max)
    colnames(max_step) <- c("item_step", "max")
    max_step <- as.matrix(max_step)
    # y design
    y_design <- as.matrix( data.frame(1:N, y+1) )

    #-- output
    res <- list( design=design, design_labels=design_labels, weights=weights,
                    N=N, data=data, des_names=des_names, suff_stat=suff_stat,
                    parm=parm, parm_sign=parm_sign, center=center, center0=center0,
                    maxcat=maxcat, K=maxcat+1, des_names1=des_names1, is_step=is_step,
                    max_item=max_item, max_step=max_step, y_design=y_design )
    return(res)
}
