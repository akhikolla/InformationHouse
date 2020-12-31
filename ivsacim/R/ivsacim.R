#' @title Fitting a Cumulative Intensity Model for Exposure Effect with Instrumental Variables
#' @description ivsacim is used to fit cumulative intensity models for exposure effects with instrumental variables.
#' @param time the censored event time
#' @param event event indicator
#' @param instrument the instrumental variable
#' @param IV_valid whether assuming IV satisfies the exclusion restriction
#' @param treatment_init the initial treatment assignment
#' @param treatment_shift_time the shift time of each subject, if no shift for a subject, set as 0
#' @param covar the baseline covariates
#' @param max_time the max time that we threshold
#' @param max_time_bet the max time that we threshold
#' @param n_sim the number of resampling, set as 0 if no resampling is needed
#' @return ivsacim returns an object of class "tivsacim".
#' An object of class "ivsacim" is a list containing the following components:
#' \item{stime}{an estimate of the baseline hazards function}
#' \item{dB_D}{an estimate of the baseline hazards function}
#' \item{B_D}{an estimate of the coefficients}
#' \item{beta}{an estimate of the baseline hazards function}
#' \item{B_D_se}{an estimate of the variance covariance matrix of coef}
#' \item{beta_se}{an estimate of the baseline hazards function}
#' \item{by_prod}{a byproduct, that will used by other functions}
#' @importFrom stats fitted glm lm predict qnorm residuals vcov rnorm
#' @importFrom lava iid
#' @importFrom survival survfit
#' @export ivsacim
#' @examples
#' n = 200
#' event = rbinom(n, 1, 0.8)
#' IV = rbinom(n, 1, 0.5)
#' trt_init = IV
#' trt_shift = rep(0, n)
#' time = rexp(n)/(0.5 + trt_init * 0.2)
#' max_t = 3
#' max_t_bet = 3
#' n_sim = 0
#' fit <- ivsacim(time, event, IV, IV_valid = TRUE, trt_init, 
#' trt_shift, covar = NULL, max_t, max_t_bet, n_sim)
ivsacim <- function (time, 
                     event, 
                     instrument,
                     IV_valid = TRUE,
                     treatment_init, 
                     treatment_shift_time = NULL, 
                     covar = NULL, 
                     max_time = NULL, 
                     max_time_bet = NULL,
                     n_sim = 0) {
  
  if (is.null(max_time)) {
    max_time = max(time)
  }
  
  if (is.null(treatment_shift_time)) {
    treatment_shift_time = rep(0, length(time))
  }
  
  stime1 = sort(time)
  event_new = event
  event_new[time > max_time] = 0
  stime = sort(time[event_new == 1])
  
  iv_centered <- IV_center(instrument, covar)
  Zc <- iv_centered$Zc
  epstheta <- iv_centered$epstheta
  Edot <- iv_centered$Edot
  pdim <- iv_centered$pdim
  n = length(time)
  k = length(stime)
  D_status <- treatment_status(n, k, stime, treatment_init, treatment_shift_time, max_time) 
  #pdim <- nrow()
  
  if (IV_valid) {
    res <- ivsacim_est(time, event, stime, Zc, D_status, epstheta, Edot)
  } 
  else {
    print("invalid_IV")
  }
  
  eps = t(res$by_prod$eps)
  k = length(stime < max_time_bet)
  eps = apply(eps, 2, cumsum)
  eps.beta = res$by_prod$eps_beta
  
  if (n_sim > 0){
    ant.resamp = n_sim
    GOF.resam0 = matrix(0, nrow = ant.resamp, ncol = k)
    max.proc0 = numeric(ant.resamp)
    ## Const. effect
    GOF.resam = matrix(0, nrow = ant.resamp, ncol = k)
    max.proc = CvM.proc = numeric(ant.resamp)
    eps.const.eff = eps[, ]
    for(j in 1:k){
      eps.const.eff[j, ] = eps[j, ] - stime[j] * eps.beta
    }
    
    for(j1 in 1:ant.resamp){
      G = rnorm(n, 0, 1)
      tmp.mat0 = eps %*% matrix(G, n, 1)
      GOF.resam0[j1, ] = c(tmp.mat0)
      max.proc0[j1] = max(abs(GOF.resam0[j1, ]))
      tmp.mat = eps.const.eff %*% matrix(G, n, 1)
      GOF.resam[j1, ] = c(tmp.mat)
      max.proc[j1] = max(abs(GOF.resam[j1, ]))
      CvM.proc[j1] = sum(GOF.resam[j1, ]^2 * c(diff(stime), max_time_bet - stime[k]))
    }
    
    GOF.proc0 = res$B_D[1:k]
    max.obs0 = max(abs(GOF.proc0))
    GOF.proc = res$B_D[1:k] - res$beta * stime
    max.obs = max(abs(GOF.proc))
    CvM.obs = sum(GOF.proc^2 * c(diff(stime), max_time_bet - stime[k]))
    pval_0 = sum(max.proc0 > max.obs0)/ant.resamp
    pval.sup = sum(max.proc > max.obs)/ant.resamp
    pval.CvM = sum(CvM.proc > CvM.obs)/ant.resamp
    res$pval_0 = pval_0
    res$pval_GOF_sup = pval.sup
    res$pval_GOF_CvM = pval.CvM
    if (n_sim > 50){
      res$GOF.resamp = rbind(stime[1:k], GOF.proc, GOF.resam[1:50, ])
    }
  }
  res$by_prod$noresampling = (n_sim <= 0)
  res$by_prod$max_time = max_time
  res$by_prod$max_time_bet = max_time_bet
  res$by_prod$D_status = D_status
  res$by_prod$k = k
  class(res) <- "ivsacim"
  return(res)
  
}




#' @title Summarizing Cumulative Intensity Function of Treatment with Instrumental Variables 
#' Estimation Using Structural Additive Cumulative Intensity Models 
#' @param object an object of class "ivsacim", usually, a result of a call to ivsacim.
#' @param x an object of class "summary.ivsacim", usually, a result of a call to summary.ivsacim.
#' @param ... further arguments passed to or from other methods.
#' @description summary method for class "ivsacim".
#' @export ivsacim
#' @export summary.ivsacim
#' @importFrom stats pnorm pchisq
#' @method summary ivsacim
#' @S3method summary ivsacim
#' @examples
#' n = 200
#' event = rbinom(n, 1, 0.8)
#' IV = rbinom(n, 1, 0.5)
#' trt_init = IV
#' trt_shift = rep(0, n)
#' time = rexp(n)/(0.5 + trt_init * 0.2)
#' max_t = 3
#' max_t_bet = 3
#' n_sim = 0
#' fit <- ivsacim(time, event, IV, IV_valid = TRUE, trt_init, 
#' trt_shift, covar = NULL, max_t, max_t_bet, n_sim)
#' summary(fit)
#' @details print.summary.ivsacim tries to be smart about formatting coefficients, an estimated variance covariance matrix of
#' the coeffieients, Z-values and the corresponding P-values.
#' @return The function summary.ivsacim computes and returns a list of summary statistics of the fitted model given in object.
summary.ivsacim <- function(object, digits = 3, ...){
  obj <- object
  chisq_beta = (obj$beta/obj$beta_se)^2
  pval_beta = 1 - pchisq(chisq_beta, df = 1)
  
  res <- list()
  res$pval_0 = obj$pval_0
  res$beta = obj$beta
  res$beta_se = obj$beta_se
  res$pval_beta = pval_beta
  res$pval_GOF_sup = obj$pval_GOF_sup
  if (obj$by_prod$noresampling) {
    res$noresult = TRUE
  }
  else {
    res$noresult = FALSE
  }
  class(res) = "summary.ivsacim"
  res
}


#' @rdname summary.ivsacim
#' @param digits number of digits we want to show
#' @S3method print summary.ivsacim
print.summary.ivsacim <- function(x, digits = 3, ...){
  obj <- x
  # We print information about object:  
  cat("   \n")
  cat("Instrumental Variables Structural Additive Cumulative Intensity Model \n\n")
  if (!obj$noresult) { 
    
    #if (sum(obj$conf.band)==FALSE)  mtest<-FALSE else mtest<-TRUE; 
    #if (mtest==FALSE) cat("Test not computed, sim=0 \n\n")
    #if (mtest==TRUE) { 
    test0 <- cbind(obj$pval_0)
    #testC<-cbind(obj$obs.testBeqC,obj$pval.testBeqC) 
    colnames(test0) <- c("Supremum-test pval")
    #rownames(test0)<-colnames(prop.excess.object$cum)[-1]
    rownames(test0) <- c("Exposure")
    #rownames(test0)<- c("Exposure")
    #colnames(testC)<- c("sup| B(t) - (t/tau)B(tau)|","p-value H_0: B(t)=b t")  
    cat("Test for non-significant exposure effect. H_0: B_D(t) = 0 \n")
    cat("   \n")
    #cat("Test for non-parametric term, H_0: B(t)=0  \n")
    prmatrix(signif(test0, digits))
    cat("   \n")
    cat("\n")
    test_gof = cbind(obj$pval_GOF_sup)
    colnames(test_gof) <- c("Supremum-test pval")
    rownames(test_gof) <- c("        ")
    cat("Goodness-of-fit test for constant effects model\n")
    prmatrix(signif(test_gof, digits))
    cat("   \n")
    
  }
  
  cat("Constant effect model  \n")
  res <- cbind(obj$beta, obj$beta_se)#,diag(obj$robvar.gamma)^.5,diag(obj$D2linv)^.5)
  z <- obj$beta/obj$beta_se
  pval <- obj$pval_beta
  res <- cbind(res, z, pval)
  rownames(res) <- c("Exposure     ")
  colnames(res) <- c("coef", "se(coef)", "z", "p")#,#Robust Std.  Error","D2log(L)^-(1/2)")  
  prmatrix(signif(res, digits))
  cat("   \n") 
  
  
  #cat("  Call: ")
  #dput(attr(obj, "Call"))
  #cat("\n")
  invisible()
}

#' @title Plotting Estimated Cumulative Intensity function with Pointwise Confidence Intervals
#' @param x the fitting object after fitting IVSACIM model
#' @param gof whether to draw the goodness-of-fit plot
#' @param ... the other arguments you want to put in the built-in plot function
#' @description The function will plot the estimated cumulative intensity function of the treatment after fitting.
#' Corresponding pointwise confidence intervals at level alpha are also included.
#' @export ivsacim
#' @export plot.ivsacim
#' @importFrom graphics plot lines abline
#' @method plot ivsacim
#' @S3method plot ivsacim
#' @return No return value, called for side effects
#' @examples
#' n = 200
#' event = rbinom(n, 1, 0.8)
#' IV = rbinom(n, 1, 0.5)
#' trt_init = IV
#' trt_shift = rep(0, n)
#' time = rexp(n)/(0.5 + trt_init * 0.2)
#' max_t = 3
#' max_t_bet = 3
#' n_sim = 100
#' fit <- ivsacim(time, event, IV, IV_valid = TRUE, trt_init, 
#' trt_shift, covar = NULL, max_t, max_t_bet, n_sim)
#' plot(fit, main = "", xlab = "Time", ylab = "Cumulative Intensity Function")
#' plot(fit, gof = TRUE, xlab = "Time", ylab = "")
plot.ivsacim <- function (x, gof = FALSE, ...){
  
  obj <- x
  if (!inherits(obj, 'ivsacim')) 
    stop ("Must be an ivsacim object")
  
  if(gof){
    #par(mfrow=c(1,1))
    minv = min(c(obj$GOF.resamp[2, ], obj$GOF.resamp[2:22, ]))
    maxv = max(c(obj$GOF.resamp[2, ], obj$GOF.resamp[2:22, ]))
    plot(obj$GOF.resamp[1, ], obj$GOF.resamp[2, ], type = "n", ylim = c(minv, maxv), main = "Goodness of fit process")
    for(j in 1:20){
      lines(obj$GOF.resamp[1, ], obj$GOF.resamp[2 + j, ], type = "s", col = "grey")   
    }
    lines(obj$GOF.resamp[1, ], obj$GOF.resamp[2, ], type = "s", lty = 1, lwd = 3)
    abline(0, 0)	
    
  }else{
    B_D <- obj$B_D
    B_D_se <- obj$B_D_se
    beta = obj$beta
    beta_se = obj$beta_se
    stime = obj$stime
    minv1 = min(c(B_D - 1.96 * B_D_se, beta * stime))
    maxv1 = max(c(B_D + 1.96 * B_D_se, beta * stime))
    
    #par(mfrow=c(1,1))
    plot(stime, B_D, type = "s", ylim = c(minv1, maxv1))
    lines(stime, c(B_D + 1.96 * B_D_se), type = 's', lty = 2)
    lines(stime, c(B_D - 1.96 * B_D_se), type = 's', lty = 2)
    abline(0, 0)
    lines(stime, beta * stime)
    lines(stime, (beta + 1.96 * beta_se) * stime, lty = 2)
    lines(stime, (beta - 1.96 * beta_se) * stime, lty = 2)
  }
  invisible()
}


