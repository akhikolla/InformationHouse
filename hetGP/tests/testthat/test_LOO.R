library(hetGP)
context("LOO")


test_that("LOO",{
  
  set.seed(32)
  ## motorcycle data
  library(MASS)
  X <- matrix(mcycle$times, ncol = 1)
  Z <- mcycle$accel
  nvar <- 1
  
  ## Start with GP models
  for(modelfun in c("mleHomGP", "mleHetGP")){
    for(trend in c(NA, 0)){
      
      if(is.na(trend)) trend <- NULL
      ## Model fitting
      model <- match.fun(modelfun)(X = X, Z = Z, lower = rep(0.1, nvar), upper = rep(10, nvar),
                                   covtype = "Matern5_2", known = list(beta0 = trend))
      LOO_p <- LOO_preds(model)
      
      # model minus observation(s) at x_i
      d_mot <- find_reps(X, Z)
      
      LOO_ref <- matrix(NA, nrow(d_mot$X0), 2)
      for(i in 1:nrow(d_mot$X0)){
        model_i <- match.fun(modelfun)(X = list(X0 = d_mot$X0[-i,, drop = FALSE], Z0 = d_mot$Z0[-i],
                                                mult = d_mot$mult[-i]), Z = unlist(d_mot$Zlist[-i]),
                                       lower = rep(0.1, nvar), upper = rep(50, nvar), covtype = "Matern5_2",
                                       known = list(theta = model$theta, k_theta_g = model$k_theta_g, g = model$g,
                                                    Delta = model$Delta[-i], beta0 = trend))
        model_i$nu_hat <- model$nu_hat
        
        # For hetGP, need to use the same Lambdas to get the same results 
        if(modelfun == "mleHetGP"){
          model_i$Lambda <- model$Lambda[-i]
          model_i <- strip(model_i)
          model_i <- rebuild(model_i)
        }   
        
        p_i <- predict(model_i, d_mot$X0[i,,drop = FALSE])
        LOO_ref[i,] <- c(p_i$mean, p_i$sd2)
      }
      # Compare results
      expect_equal(LOO_ref[,1], as.numeric(LOO_p$mean), tol = 1e-6)
      expect_equal(LOO_ref[,2], as.numeric(LOO_p$sd2), tol = 1e-5)
    }
  }
  
  ## Then TP versions
  for(modelfun in c("mleHomTP", "mleHetTP")){

      ## Model fitting
      model <- match.fun(modelfun)(X = X, Z = Z, lower = rep(0.1, nvar), upper = rep(10, nvar),
                                   covtype = "Matern5_2", known = list(beta0 = 0))
      LOO_p <- LOO_preds(model)
      
      # model minus observation(s) at x_i
      d_mot <- find_reps(X, Z)
      
      LOO_ref <- matrix(NA, nrow(d_mot$X0), 2)
      for(i in 1:nrow(d_mot$X0)){
        model_i <- match.fun(modelfun)(X = list(X0 = d_mot$X0[-i,, drop = FALSE], Z0 = d_mot$Z0[-i],
                                                mult = d_mot$mult[-i]), Z = unlist(d_mot$Zlist[-i]),
                                       lower = rep(0.1, nvar), upper = rep(50, nvar), covtype = "Matern5_2",
                                       known = list(theta = model$theta, k_theta_g = model$k_theta_g, g = model$g,
                                                    sigma2 = model$sigma2, nu = model$nu,
                                                    Delta = model$Delta[-i], beta0 = trend))

        model_i$psi <- model$psi # psi is taken as fixed
        # For hetTP, need to use the same Lambdas and psi to get the same results 
        if(modelfun == "mleHetTP"){
          model_i$Lambda <- model$Lambda[-i]
          model_i <- strip(model_i)
          model_i <- rebuild(model_i)
        }   
        
        p_i <- predict(model_i, d_mot$X0[i,,drop = FALSE])
        LOO_ref[i,] <- c(p_i$mean, p_i$sd2)
      }
      # Compare results
      expect_equal(LOO_ref[,1], as.numeric(LOO_p$mean), tol = 1e-6)
      expect_equal(LOO_ref[,2], as.numeric(LOO_p$sd2), tol = 1e-6)
    }
  
})