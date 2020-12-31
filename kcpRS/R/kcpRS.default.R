
#' @export
kcpRS.default <-function(data,
                         RS_fun,
                         RS_name,
                         wsize = 25,
                         nperm = 1000,
                         Kmax = 10,
                         alpha = .05,
                         varTest = FALSE,
                         ncpu = 1) {
  myfun <-
    function(data, wsize, FUN) {
      FUN(data, wsize)
    } #functions that allows the user to input RS_fun

  data <- as.data.frame(data)
  data <- scale(data)

  if (wsize >=nrow(data)) {
    warning("Invalid window size")
  }

  if (nperm>0 & Kmax <= 2) {
    warning("K should be greater than 2 to carry out the permutation test.")
  }

  if (nperm !=0 & nperm <50 ) {
    warning("Very low number of permuted data sets selected for the permutation test.")
  }

  if (ncpu > detectCores() ) {
    warning(paste0("Maximum value for ncpu is ",detectCores(),"."))
  }

  if (ncpu<=detectCores()){

    #compute running statistic
    RS <- myfun(data, wsize, FUN = RS_fun)

    if (nrow(RS) <= Kmax) {
      warning("Not enough number of windows for the requested Kmax.")
    }

    #find change point locations given K and the H matrix
    kcp_RS_result <- kcpa(RunStat = RS, Kmax, wsize)

    alp <- ifelse(nperm == 0, NA, alpha)
    p_var_test = NA         #added functionality: if nperm is 0, permutation test is not done
    p_varDrop_test = NA
    subTest_alpha = NA
    cps <- NULL
    best_k <- 0

    #permutation test
    if (nperm > 0) {

      cl <- makeCluster(ncpu)
      registerDoParallel(cl)

      perm_test_result <- permTest(data, RS_fun, wsize, nperm, Kmax, alpha, varTest)
      p_var_test = ifelse(isTRUE(varTest), perm_test_result$p_var_test, "NA")
      p_varDrop_test = perm_test_result$p_varDrop_test
      subTest_alpha = ifelse(isTRUE(varTest), alpha / 2, alpha)

      #grid search if there is at least 1 CP
      if (perm_test_result$sig == 1) {
        RS = as.matrix(RS)
        N <- nrow(RS)
        Rmin <- kcp_RS_result$Rmin

        #compute vmax
        lower <- ceiling(.05 * nrow(RS))
        upper <- floor(.95 * nrow(RS))
        a <- cov(as.matrix(RS[1:lower,]))
        b <- cov(as.matrix(RS[upper:nrow(RS),]))
        a <- sum(diag(a))
        b <- sum(diag(b))

        if (lower == 1)
          #in this case only use the trivial bound of vmax <<- 1 bc of kernel
          a <- 1
        vmax <- max(a, b)

        #compute the increment to be used for the search by finding C_zero
        n_searches <- 10000
        c_zero <- NULL
        pt_0 <- (vmax / N) * (1 + log(N))
        Rmin_0 <- Rmin[1]

        for (k in 1:Kmax) {
          d <- k + 1
          pt_k <- ((vmax * d) / N) * (1 + log(N / d))
          Rmin_k <- Rmin[k + 1]

          c_zero[k] <- (Rmin_k - Rmin_0) / (pt_0 - pt_k)
        }

        c_zero_max <- ceiling(max(c_zero))
        #we need a c larger than c_zero to ensure that the penalty is large and
        #k<-0 will be chosen in the model selection
        inc <- c_zero_max / n_searches


        #actual grid search
        c_searched <- NULL
        c <- 1       #initial c
        L <- Kmax + 1    # max dimension
        d <- Kmax + 1    # best dimension -> k + 1, if k is the optimal no. of change points

        diff <- L - d + 1    #initial value is 1
        #search is terminated if diff<-L, implying d<-1 or k<-0 (zero change points is chosen, penalty is very large)
        while (diff < L) {
          penalized_Rmin <- NULL
          for (i in 1:L) {
            pen <- (c * vmax * i / N) * (1 + log(N / i))
            penalized_Rmin[i] <- Rmin[i] + pen
          }
          d <- which.min(penalized_Rmin) #best d for the current c
          row_temp <- c(c, d)
          c_searched <- rbind(c_searched, row_temp)
          diff <- L - d + 1   #check if search is to be terminated
          c <- c + inc      #c for the next iteration
        }


        #most stable k
        c_searched <- as.data.frame(c_searched, row.names = FALSE)
        c_searched$jump <- c(0, abs(diff(c_searched[, 2])))#compute differences to find c's that determines the shift in d
        colnames(c_searched) <- c("c", "dimension", "jump")
        c_candidates <- subset(c_searched, c_searched$jump > 0)
        n_candidate_cs <- dim(c_candidates)[1]

        if (n_candidate_cs > 1) {
          c_candidates$range <- c(abs(diff(c_candidates$c)), 0)
          best_d <-
            c_candidates$dimension[which.max(c_candidates$range)]
        } else {
          best_d <- 1
        }
        best_k <- best_d - 1


        #change points given the most stable k
        if (best_k > 0) {
          cps <- kcp_RS_result[(best_k + 1), (3:(best_k + 2))]
          if (identical(RS_fun,runAR)){cps=cps+1} #add a time point forward for cps when monitoring running AR
        } else
          cps = NULL
      }
      stopCluster(cl)
    }


    #output
    output <- list(
      "RS_name" = RS_name,
      "RS" = RS,
      "wsize" = wsize,
      "varTest" = varTest,
      "nperm" = nperm,
      "alpha" = alp,
      "subTest_alpha" = subTest_alpha,
      "BestK" = best_k,
      "changePoints" = as.numeric(cps),
      "p_var_test" = p_var_test,
      "p_varDrop_test" = p_varDrop_test,
      "CPs_given_K" = kcp_RS_result
    )

    class(output)="kcpRS"
    return(output)
  }
}
