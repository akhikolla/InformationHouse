# $Id: StoNED.R 223 2020-06-20 10:38:10Z lao $
################################################################
# Convex nonparametric least squares
# here for convex (Cost) function
# with multiplicative error term:
# Y=b'X*exp(e)
# and additive error term:
# Y=b'X + e
################################################################

###############################################################
# USAGE: X are RHS, Y are LHS,
# COST specifies whether a cost function needs is estimated (COST=1)
# or a production function (COST=0),
# MULT determines if multiplicative (MULT=1) or additive (MULT=0)
# is estimated
# RTS determines returns to scale assumption: RTS="vrs", "drs", "irs"
# and "crs" are possible for constant or variable returns to scale
# METHOD specifies the way efficiency is estimated: MM for Method
# of Momente and PSL for pseudo likelihood estimation are possible
###############################################################

require(ucminf)
require(quadprog)


stoned <- function(X, Y, RTS="vrs", COST=0, MULT=0, METHOD="MM") {  

    ###############################################
    ###### row, columns and intercept for non-crs
    ###############################################
    
    if (!is.matrix(X)) X <- as.matrix(X)
    if (!is.matrix(Y)) Y <- as.matrix(Y)
    
    n <- nrow(X)   # Number of units
    m <- ncol(X)   # Number of parameters in each unit except for intercept

    rts <- c("vrs","drs","crs","irs")
    if ( missing(RTS) ) RTS <- "vrs" 
    if ( is.numeric(RTS) )  {
        RTStemp <- rts[RTS] # there is no first fdh
        RTS <- RTStemp
    }
    RTS <- tolower(RTS)
    if ( !(RTS %in% rts) )  stop("Unknown scale of returns:", RTS)

    if (!is.null(METHOD) & !(METHOD %in% c("MM", "PSL"))) {
        stop("Unknown choice for METHOD: ", METHOD)
    }

    inter <- ifelse (RTS=="crs", 0, 1)
    if (is.null(METHOD))  METHOD <- "NONE"
    if (RTS=="crs") MULT <- 1
    if (RTS=="drs") MULT <- 1
    if (RTS=="irs") MULT <- 1
    if (RTS=="crs" | RTS=="irs" | RTS=="drs") {
        message("Multiplicative model induced for this scale assumption")
    }


    ###############################################
    ### Get data together, create A - matrix with
    ### intercept if necessary and B-Matrix
    ###############################################
    if (MULT==1) {
        if (inter==1) {
            Z <- diag(1 / c(Y))
            for (i in 1:m) {
                # Z <- cbind( Z, (diag(n) * as.matrix(X)[,i])/Y )
                Z <- cbind(Z, diag(c(X[,i] / Y)))
            }
        }
        else if (inter==0) {
            Z <- NULL
            for (i in 1:m) {
                Z <- cbind(Z, diag(c(X[,i] / Y)))
            }
        }
        B <- rep(1, n)
    }

    if (MULT==0) {
        if (inter==1) {
            Z <- diag(n)
            for (i in 1:m) {
                Z <- cbind(Z, diag(X[,i]))
            }
        }
        else if (inter==0) {
            Z <- NULL
            for (i in 1:m) {
                Z <- cbind(Z, diag(X[,i]))
            }
        }
        B <- Y
    }
    
    colnames(Z) <- outer(1:n, ifelse(inter==1, 0, 1):m, FUN = function(x, y) {
        paste(x, y, sep="_")
    })

    ###############################################
    ###### Create non-neg constraints
    ###############################################
    non_neg_beta <- cbind(matrix(0, nrow = m * n, ncol=n), diag(m * n))
    
    if (inter==1) {
        if (RTS=="vrs") {
            sign_cons <- non_neg_beta
        }
        else if (RTS=="drs") {
            RTS_cons <- cbind(ifelse(COST==0,1,-1) * diag(m * n), matrix(0, nrow = m * n, ncol = n))
            sign_cons <- rbind(non_neg_beta, RTS_cons)
        }
        else if (RTS=="irs") {
            RTS_cons <- cbind(ifelse(COST==0,-1,1) * diag(m * n), matrix(0, nrow = m * n, ncol = n))
            sign_cons <- rbind(non_neg_beta, RTS_cons)
        }
    }
    else if (inter==0) {
        ## 'crs' er uden konstantled
        sign_cons <- diag(m * n)
    }


    ###############################################
    ###### Create concavity constraints
    ###############################################
    alpha_matlist <- list()
    beta_matlist <- list()

    for (i in 1:n) {
        a_matr <- diag(n)
        a_matr[,i] <- a_matr[,i] - 1
        alpha_matlist[[i]] <- a_matr

        # Create beta
        b_matr <- matrix(0, ncol = m * n, nrow = n)

        for (j in 1:m) {
            b_matr[, i + n * (j - 1)] <- -X[i,j]
            for (k in 1:n) {
                b_matr[k, k + (j - 1) * n] <- b_matr[k, k + (j - 1) * n] + X[i,j]
            }
        }
        beta_matlist[[i]] <- b_matr
    }

    if (inter==1) {
        sspot_cons <- cbind(do.call(rbind, alpha_matlist), do.call(rbind, beta_matlist))
    }
    else if (inter==0) {
        sspot_cons <- cbind(do.call(rbind, beta_matlist))
    }
    
    
    #####################################
    # !!!!!!!!! Cost function !!!!!!!!!
    # Convexity - not concavity constr!
    #####################################
    if (COST==1) {
        sspot_cons <- -sspot_cons
    }


    #####################################
    ### Getting the constraints together
    #####################################
    G <- rbind(sign_cons, sspot_cons)
    h <- rep(0, nrow(sign_cons) + nrow(sspot_cons))

    # Set column and row names for the constrains; should reuse any existing names
    colnames(G) <- outer(1:n, ifelse(inter==1, 0, 1):m,
                        FUN = function(x, y) {paste(x, y, sep = "_")})
    rownames(G) <- names(h) <-
        c(outer(1:n, 1:m, FUN = function(x, y) { paste(x, y, sep = "_") }),
          outer(1:n, 1:n, FUN = function(x, y) { paste(y, x, sep = "")  }))

    #####################################
    ### Run
    #####################################


    # Drop restriktioner hvor hele raekken er nul
    helenul <- apply(G, 1, function(x) {all(x==0)})
    G <- G[!helenul, ]
    h <- h[!helenul]

    dvec  <- base::crossprod(Z, B) # was: t(A)%*%B
    Dmat  <- base::crossprod(Z, Z) # t(A) %*% A
    # Make sure Dmat is positive definite
    Dmat <- Dmat + diag(1e-10 * max(Dmat), dim(Dmat)[1], dim(Dmat)[2])

    z_cnls   <- solve.QP(Dmat, dvec, Amat = t(G), bvec = h)
    z_cnls$type <- "solve.QP"
    z_cnls$X <- z_cnls$solution
    value <- 0.5 * z_cnls$solution %*% Dmat %*% z_cnls$solution - c(dvec) %*% z_cnls$solution

    z_cnls$IsError <- FALSE

    slack  <-  G %*% z_cnls$X - h
    slackResidual <-  -sum(slack[slack < 0])
    z_cnls$slackNorm <- slackResidual
    eps2 <- sum((Z %*% z_cnls$X - B)^2)
    z_cnls$solutionNorm     <- eps2


    #####################################
    ### Collect results
    #####################################
    
    # First column is the intercept when 'inter==1'
    beta_matrix <- matrix(z_cnls$X, ncol = length(z_cnls$X) / n)
    if (inter==0) {
        y_hat <- diag(beta_matrix %*% t(X))
    }
    else if (inter==1) {
        y_hat <- diag(beta_matrix %*% t(cbind(1, X)))
    }

    if (MULT==1) {
        resid <- log(Y) - log(y_hat)
        # Hvorfor baade 'yhat' og 'fitted'?
        fitted <- Y / exp(resid)
    }
    else if (MULT==0) {
        resid <- Y - y_hat
        fitted <- Y - resid
    }

    # Restled paa en anden maade
    restled <- Y - Z %*% z_cnls$X
    SSR <- crossprod(restled)  # t(restled) %*% restled


    #####################################
    ### Inefficiency Estimation
    #####################################
    if (METHOD=="MM") {
        M2 <- sum((resid - mean(resid)) * (resid - mean(resid))) / length(resid)
        M3 <- sum((resid - mean(resid)) * (resid - mean(resid)) * 
                  (resid - mean(resid))) * (1 / length(resid))

        if (COST==1) M3 <- -M3
        if (M3 > 0) print("wrong skewness? Third moment of the residuals is greater than 0")
        # Sigma_u and Sigma_v as functions of M3 and M2
        sigma_u <- (M3 / (sqrt(2 / pi) * (1 - 4 / pi)))^(1 / 3)
        sigma_v <- sqrt(M2 - ((pi - 2) / pi) * sigma_u * sigma_u)
        lamda <- sigma_u / sigma_v
        myy <- (sigma_u^2 / pi)**(1 / 2)
    }  ## if (METHOD=="MM")

    if (METHOD=="PSL") {
        # minus log-likelihood funktion
        ll <- function(lmd = 1, e) {
            sc <- sqrt(2 * lmd^2 / pi / (1 + lmd^2))
            si <- sqrt(mean(e^2) / (1 - sc^2))
            mu <- si * sc
            ep <- e - mu
            likelihood <- -length(ep) * log(si) + 
                            sum(pnorm(-ep * lmd / si, log.p = TRUE)) -
                            0.5 * sum(ep^2) / si^2
            return(-likelihood)
        }

        if (COST==0) {
            sol <-  ucminf(par = 1, fn = ll, e = resid)
        }
        else if (COST==1) {
            sol <-  ucminf(par = 1, fn = ll, e = -resid)
        }
        if (sol$convergence < 0) {
            warning("Not necessarily converged, $convergence = ",
                sol$convergence, "\n", sol$message)
        }
        lamda <- sol$par

        ## Use estimate of lambda to calculate sigma^2
        sc <- sqrt(2 * lamda^2 / pi / (1 + lamda^2))
        sig.sq <- mean(resid^2) / (1 - sc^2)
        ## Now calculate sigma.u and sigma.v
        sigma_v <- sqrt(sig.sq / (1 + lamda^2))
        sigma_u <- sigma_v * lamda
    }  ## if (METHOD=="PSL")

    if (METHOD=="NONE") {
        # Hvis nedenstaaende variabler ikke saettes til NA giver det 
        # fejl ved manglende variabler senere i programmet.
        sigma_v <- sigma_u <- NA
    }

    # Shift factor
    if (MULT==0) {
        Expect_Ineff <- sigma_u * sqrt(2 / pi)
    }
    else if (MULT==1) {
        Expect_Ineff <- 1 - exp(-sigma_u * sqrt(2 / pi))
    }
    sumsigma2 <- sigma_u^2 + sigma_v^2

    #### JMLS Point estimator for inefficiency ####
    # Composite error
    if (COST==0) {
        Comp_err <- resid - sigma_u * sqrt(2 / pi)
    }
    else if (COST==1) {
        Comp_err <- resid + sigma_u * sqrt(2 / pi)
    }
    # mu star and sigma star
    mu_star <- -Comp_err * sigma_u * sigma_u / sumsigma2
    sigma_star <- sigma_u * sigma_u * sigma_v * sigma_v
    cond_mean <- mu_star + (sigma_u * sigma_v) *
        (dnorm(Comp_err / sigma_v * sigma_v) / (1 - pnorm(Comp_err / sigma_v * sigma_v)))
    # Be aware: for the calculation of the cond mean I follow Keshvari / Kuosmanen 2013
    # stating that the formula in Kuosmanen/Kortelainen is not correct!

    if (MULT==1) {
        eff_score <- exp(-cond_mean)
        if (COST==0) {
            Frontier_reference_points <- fitted * exp(sigma_u * sqrt(2 / pi))
        }
        else if (COST==1) {
            Frontier_reference_points <- fitted * exp(-sigma_u * sqrt(2 / pi))
        }
    }
    else if (MULT==0) {
        eff_score <- 1 - cond_mean / Y
        if (COST==1) {
            Frontier_reference_points <- fitted - sigma_u * sqrt(2 / pi)
        }
        else if (COST==0) {
            Frontier_reference_points <- fitted + sigma_u * sqrt(2 / pi)
        }
    }
    

    return(list(residualNorm = z_cnls$slackNorm,
        solutionNorm = z_cnls$solutionNorm,
        error = z_cnls$IsError,
        # type = z_cnls$type,
        coef = beta_matrix,
        sol = z_cnls$X,
        residuals = resid,
        restled = restled,
        fit = fitted,
        yhat = y_hat,
        eff = eff_score,
        front = Frontier_reference_points,
        sigma_u = sigma_u,
        SSR = SSR
        # G=G,
        # H=h
        # Fitted_Values = yhat
    ))
}  ## stoned
