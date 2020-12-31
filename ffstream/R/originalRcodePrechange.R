#' CUSUM change detection for a stream in R with known prechange parameters
#'
#' Original implementation in R of CUSUM change detector, but now
#' expecting the prechange mean and variance to be specified.
#'
#'
#' @param stream The stream of observations.
#'
#' @param BL The burn-in length - this won't actually be used, but is kept
#'           for historical reasons.
#'
#' @param params A list of parameters for the CUSUM algorithm. Consists of
#'              \describe{
#'                  \item{\code{d}}{A control parameter also known as 
#'                                  \eqn{k}.}
#'
#'                  \item{\code{B}}{A control parameter also known as 
#'                                  \eqn{h}.}
#'              }
#'
#' @param mu0 The prechange mean, which is assumed known in this context
#'
#' @param sigma0 The prechange standard deviation, which is assumed known 
#'              in this context
#'
#'
#' @return A vector of estimated changepoints.
#'
#'
#' @section Author:
#' Dean Bodenham
#'
#'
#' @section References:
#' D. A. Bodenham and N. M. Adams (2016) 
#' \emph{Continuous monitoring for changepoints in data 
#' streams using adaptive estimation}. 
#' Statistics and Computing  
#' doi:10.1007/s11222-016-9684-8
#'
#'
#' @keywords internal
CUSUM_stream_jumpdetect_prechange <- function(stream, BL, params, mu0, sigma0){

	d <- params[[1]]
	B <- params[[2]]
	
	#will run until stream finished...
	streampos <- 1
	N <- length(stream)
	detected_count <- 0
	#vector for saving jumps, will trim later
	#just a quick way of saving space - can't have more than BL*det_count obs in stream,
    #only one changepoint
	M <- 1
	detect_pos_vec <- rep(0, M)
	
    #jump_found is a boolean that flags if a jump is detected
    jump_found <- FALSE
	while ((streampos < N) && (jump_found==FALSE)){
		
		#set values for CUSUM mean and variance
		mean_j <- mu0
		var_j <- sigma0^2
		sd_j <- sqrt(var_j)
		
		nu_j <- d*sd_j
		control_j <- B*sd_j
		
		#----------end of Phase 1: burn-in----------#
		#--------------------------------------------------------------------#
		#Phase 2: detect change
		#CUSUM starts, sample from distribution 1
		S_j <- 0
		T_j <- 0
		isOutOfLimitsBool <- FALSE
		
		#j starts at 0
		j <- 0
		timedetected <- N
		while ((jump_found==FALSE) && (streampos < N)) {
			
			#update j, then last run is when j=N
			#which starts when j=N-1
			j <- j+1
			
			#get the next observation from stream
			x_new <- get_nextobs_fromstream(stream, streampos)
			streampos <- update_streampos(streampos)
			
			#DO NOT update control parameters
			#update cusums
			S_j <- S_j + x_new - mean_j - nu_j
			S_j <- S_j*(S_j>0)
			
			T_j <- T_j + mean_j - x_new - nu_j
			T_j <- (T_j)*(T_j>0)
			
			#check if outside control	
			S_isOut <- (S_j > control_j)
			T_isOut <- (T_j > control_j)
			
			#check if there is a jump
			jump_found <- S_isOut | T_isOut
			
			if (jump_found==TRUE){
				detected_count <- detected_count + 1
                #should use streampos-1
				detect_pos_vec[detected_count] <- streampos
			}


		} # end of while ((jump_found==FALSE) && (streampos < N))
		
		#end of Phase 2 - now restarting burn-in
		
	} #end while (streampos < N)
	
	#trim detect_pos_vec
	detect_pos_vec <- detect_pos_vec[1:detected_count]
	return(detect_pos_vec)
	
} # end of CUSUM detect



#' EWMA change detection for a stream in R with known prechange parameters
#'
#' Original implementation in R of EWMA change detector, but now
#' expecting the prechange mean and variance to be specified.
#'
#'
#' @param stream The stream of observations.
#'
#' @param BL The burn-in length - this won't actually be used, but is
#'           kept for historical reasons.
#'
#' @param params A list of parameters for the EWMA algorithm. Consists of
#'              \describe{
#'                  \item{\code{r}}{A control parameter which controls
#'                                  the rate of downweighting.}
#'
#'                  \item{\code{L}}{A control parameter which determines 
#'                                  the width of the control limits.}
#'              }
#'
#' @param mu0 The prechange mean, which is assumed known in this context
#'
#' @param sigma0 The prechange standard deviation, which is assumed known 
#'              in this context
#'
#'
#' @return A vector of estimated changepoints.
#'
#'
#' @keywords internal
EWMA_stream_jumpdetect_prechange <- function(stream, BL, params, mu0, sigma0){

	#ewma params: c(r_1, L_1)
	#also r =alpha, L=beta
	r <- params[[1]]
	L <- params[2]
	
	#will run until stream finished...
	streampos <- 1
	N <- length(stream)
	detected_count <- 0
	#vector for saving jumps, will trim later
	#just a quick way of saving space - can't have more than BL*det_count obs in stream,
	#since every det restarts BL
	detect_pos_vec <- rep(0, 1)
	
	delta <- r / (2-r)
    rFactorSigmaZ <- 1
    sigmaZ <- 1

    jump_found <- FALSE
	while ((streampos < N) && (jump_found==FALSE)){
		
		#set values for EWMA mean and variance
		mu_1 <- mu0
		sigma1_sq <- sigma0^2
		sigma_1 <- sqrt(sigma1_sq)
		#set EWMA control limits
		#UL: upperlimit
		UL <- mu_1 + L * sigma_1 * delta
		#LL: lowerlimit
		LL <- mu_1 - L * sigma_1 * delta
		
		#----------end of Phase 1: burn-in----------#
		#--------------------------------------------------------------------#
		#Phase 2: detect change
		#jump_found is a boolean that flags if a jump is detected
		jump_found <- FALSE
		isOutOfLimitsBool <- FALSE
		
		timedetected <- N
        #W <- mean_new
        W <- 0
        delta <- r/(2-r)
        rFactorSigmaZ <- 1
        sigmaZ <- 1
		while ((jump_found==FALSE) && (streampos < N)) {
			
			#get the next observation from stream
			x_new <- get_nextobs_fromstream(stream, streampos)
			streampos <- update_streampos(streampos)
			
			#DO NOT update control parameters
			#the ewma updating step
			W <- r * x_new + (1-r)*W
            rFactorSigmaZ <- rFactorSigmaZ * (1-r)^2
            sigmaZ <- sqrt( delta * (1-rFactorSigmaZ) ) * sigma_1
            UL <- mu_1 + L * sigmaZ
            #LL: lowerlimit
            LL <- mu_1 - L * sigmaZ



			
			#check if there is a jump
			if((W > UL) | (W < LL)){
				jump_found <- TRUE
				detected_count <- detected_count + 1
				detect_pos_vec[detected_count] <- streampos
			}
			
		} # end of while ((jump_found==FALSE) && (streampos < N))
		
		#end of Phase 2 - now restarting burn-in
		
	} #end while (streampos < N)
	
	#trim detect_pos_vec
	detect_pos_vec <- detect_pos_vec[1:detected_count]
	return(detect_pos_vec)	
}#end of EWMA jump detect






#-------------------------------------------------------------------------#
#' Change detection using the Fixed Forgetting Factor method, prechange known
#'
#' Original implementation in R of FFF change detector, but now prechange is 
#' known  
#'
#' @param stream The stream of observations.
#'
#' @param BL The burn-in length - this won't actually be used, but is kept
#'           for historical reasons.
#'
#' @param ffparams An \emph{unnamed} list of parameters for the FFF algorithm. 
#'                 Consists of:
#'              \describe{
#'                  \item{\code{lambda}}{The value of the fixed forgetting
#'                                       factor (FFF). Should be in the range 
#'                                       [0,1].}
#'
#'                  \item{\code{p}}{The value of the significance threshold, 
#'                                  which was later renamed \code{alpha} 
#'                                  (in the paper, not in this function).}
#'              
#'                  \item{\code{resettozero}}{A flag; if it zero, then the
#'                                            ffmean will be reset to zero  
#'                                            after each change. Usually set 
#'                                            to 1 (i.e. do not reset).}
#'              
#'                  \item{\code{u_init}}{The initial value of \code{u}. 
#'                                       Should be set to 0.}
#'              
#'                  \item{\code{v_init}}{The initial value of \code{v}. 
#'                                       Should be set to 0.}
#'              
#'                  \item{\code{w_init}}{The initial value of \code{w}. 
#'                                       Should be set to 0.}
#'              
#'                  \item{\code{ffmean_init}}{The initial value of the 
#'                                            forgetting factor mean,
#'                                            \code{ffmean}.
#'                                            Should be set to 0.}
#'              
#'                  \item{\code{ffvar_init}}{The initial value of the 
#'                                            forgetting factor variance, 
#'                                            \code{ffvar}.
#'                                            Should be set to 0.}
#'              
#'                      }
#'
#' @param mu0 The prechange mean, which is assumed known in this context
#'
#' @param sigma0 The prechange standard deviation, which is assumed known 
#'              in this context
#'
#' @return A vector of estimated changepoints.
#'
#'
#' @keywords internal
FFF_stream_jumpdetect_prechange <- function(stream, BL, ffparams, mu0, sigma0){

	lambda <- ffparams[[1]]
	p <- ffparams[[2]]
	resettozero <- ffparams[[3]]
	u_new <- ffparams[[4]]
	v_new <- ffparams[[5]]
	v_old <- v_new
	w_new <- ffparams[[6]]
	ffmean_new <- ffparams[[7]]
	ffvar_new <- ffparams[[8]]
	

	#will run until stream finished...
	streampos <- 1
	N <- length(stream)
	detected_count <- 0
	detect_pos_vec <- rep(0, 1)
	
    jump_found <- FALSE
	while ( (streampos < N) && (jump_found==FALSE)){
		#begin the process of burning in and then detecting
		#--------------------------------------------------------------------#
		#Phase 1: burn-in; estimating parameters
		
		#need to estimate mean and variance in burn in period
		#we set the forgetting factor to be 1
		lambda_beforeBL <- 1
		w_BL_new <- 0
		u_BL_new <- 0
		v_BL_new <- 0
		v_BL_old <- 0
		#initialize mean and variance
		mean_BL_new <- 0
		var_BL_new <- 0
		
		
		#set values for BL mean and variance
		mean_burn <- mu0
		var_burn <- sigma0^2
		burnvec <- c(p, mean_burn, var_burn)
		
		
		#--------------------------------------------------------------------#
		#Phase 2: between burn-in and change
		if (resettozero==0){
			w_new <- 0
			u_new <- 0
			v_new <- 0
			v_old <- 0
			
			#initialize mean and variance
			ffmean_new <- mu0
			ffvar_new <- sigma0^2
			
			#cat("resetting ff arl1...\n")
			
		}
		
		#--------------------------------------------------------------------#
		#Phase 2: after BL; monitoring for jump
		
		#jump_found is a boolean that flags if a jump is detected
		#jump_found <- FALSE
		
		thelimits <- c(0,0)
		
		#a boolean saying if the jump has been detected
		while((jump_found==FALSE) && (streampos < N)){
			#get the next observation from stream
			x_new <- get_nextobs_fromstream(stream, streampos)
			streampos <- update_streampos(streampos)
			
			#easier to update weight than to recalculate it
			#get BL+1 w and u
			w_new <- update_w(w_new, lambda)
			u_new <- update_u(u_new, w_new)
			#no need to update v, but we do anyway...
			v_old <- v_new
			v_new <- getv_from_u_and_w(u_new, w_new)
			
			#update mean; old mean saved for variance later
			#alternatively, could update variance first
			#just the way the derivations work out...(deeper meaning, actually...)
			ffmean_old <- ffmean_new
			ffmean_new <- update_ffmean(ffmean_old, w_new, x_new)
			
			#update variance
			ffvar_old <- ffvar_new
			ffvar_new <- update_var(ffvar_old, ffmean_old, lambda, v_old, v_new, w_new, x_new)
			
			#new function used to calculate limits efficiently, in one function, instead of in several
			#using Chebyshev
			sigma_sq_u <- u_new*var_burn
#			thelimits <- getmeanchebylimits(mean_burn, sigma_sq_u, p)
			thelimits <- getNormalLimitsFromList(burnvec, u_new)
			
			#boolean of whether or not ffmean is int he limits
			inLimitsBool <- isInLimitsNormal(ffmean_new, thelimits)
			#check for jump
			if(inLimitsBool == FALSE)
			{
				jump_found <- TRUE
				detected_count <- detected_count + 1
#				detect_pos_vec[detected_count] <- streampos
                #need to -1
				detect_pos_vec[detected_count] <- streampos-1
			}
			
		} #end of while (jumpdetected==FALSE) - will restart process, if (streampos < N) 
		
		#end of Phase 2 - now either restart or if stream is finsihed return the detect_pos_vec
		
	} #end while (streampos < N)
	
	#trim detect_pos_vec
	detect_pos_vec <- detect_pos_vec[1:detected_count]
	return(detect_pos_vec)
}
#end FFF





#-------------------------------------------------------------------------#
#' Change detection using the AFF method, using prechange mean and vairance
#'
#' Original implementation in R of the AFF, with prechange parameters
#'
#'
#' @param stream The stream of observations.
#'
#' @param BL The burn-in length - this won't actually be used, but is kept
#'           for historical reasons.
#'
#' @param affparams An \emph{unnamed} list of parameters for the FFF algorithm.
#'                 Consists of:
#'              \describe{
#'                  \item{\code{lambda}}{The value of the fixed forgetting
#'                                       factor (FFF). Should be in the range 
#'                                       [0,1].}
#'
#'                  \item{\code{p}}{The value of the significance threshold, 
#'                                  which was later renamed \code{alpha} 
#'                                  (in the paper, not in this function).}
#'              
#'                  \item{\code{resettozero}}{A flag; if it zero, then the
#'                                            ffmean will be reset to zero  
#'                                            after each change. Usually set 
#'                                            to 1 (i.e. do not reset).}
#'              
#'                  \item{\code{u_init}}{The initial value of \code{u}. 
#'                                       Should be set to 0.}
#'              
#'                  \item{\code{v_init}}{The initial value of \code{v}. 
#'                                       Should be set to 0.}
#'              
#'                  \item{\code{w_init}}{The initial value of \code{w}. 
#'                                       Should be set to 0.}
#'              
#'                  \item{\code{affmean_init}}{The initial value of the 
#'                                            forgetting factor mean,
#'                                            \code{ffmean}.
#'                                            Should be set to 0.}
#'              
#'                  \item{\code{affvar_init}}{The initial value of the 
#'                                            forgetting factor variance, 
#'                                            \code{ffvar}.
#'                                            Should be set to 0.}
#'              
#'                  \item{\code{low_bound}}{The lower bound for \code{lambda}. 
#'                                       Usually set to \code{0.6}.}
#'              
#'                  \item{\code{up_bound}}{The upper bound for \code{lambda}. 
#'                                       Usually set to \code{1}.}
#'              
#'                  \item{\code{signchosen}}{The sign used in the gradient. 
#'                                           descent. Usually set to 
#'                                           \code{-1}.}
#'              
#'                  \item{\code{alpha}}{The value of the step size in
#'                                      the gradient descent step. In
#'                                      the paper it is referred to
#'                                      as \eqn{\epsilon}.
#'                                      Usually \code{0.01}, or otherwise
#'                                      \code{0.1} or \code{0.001}.}
#'              
#'                      }
#'
#' @param mu0 The prechange mean, which is assumed known in this context
#'
#' @param sigma0 The prechange standard deviation, which is assumed known 
#'              in this context
#'
#'
#' @return A vector with the values of the adaptive forgetting factor
#'         \eqn{\overrightarrow{\lambda}}.
#'
#'
#' @keywords internal
AFF_scaled_stream_jumpdetect_prechange <- function(stream, BL, affparams, 
                                                    mu0, sigma0){

	#parameters needed to run AFF algorithm
	#alpha <- 0.01
	
	lambda <- affparams[[1]]
	p <- affparams[[2]]
	resettozero <- affparams[[3]]
	#cat("okay...")
	u_new <- affparams[[4]]
	v_new <- affparams[[5]]
	v_old <- v_new
	w_new <- affparams[[6]]
	affmean_new <- affparams[[7]]
	#need m now...
	m_new <- affmean_new*w_new 
	affvar_new <- affparams[[8]]
	low_bound <- affparams[[9]]
	up_bound  <- affparams[[10]]
	#should be minus 1, I think
	signchosen <- affparams[[11]]
	alpha <- affparams[[12]]
	#init this somehow?
	xbar_new_deriv <- 0
	Delta_new <- 0
	Omega_new <- 0
	
	#will run until stream finished...
	streampos <- 0
	N <- length(stream)
	detected_count <- 0
	#vector for saving jumps, will trim later
	#just a quick way of saving space - can't have more than BL*det_count obs in stream,
	#since every det restarts BL
	detect_pos_vec <- rep(0, 1)
	
	#this is the quantity by which we multiply L_deriv_new - the inverse of the affderiv average
	inv_AFFderiv_estim_BL <- 0
    inv_var_burn <- 0
    lambda_vec <- rep(0, N)
	
	#a vector for saving AFF
	#we deleted the saving of the adlambdavec - see testaff1.R for a version with it saved
    jumpdetected <- FALSE
	while ( (streampos < N) && (jumpdetected==FALSE) ){
		#begin the process of burning in and then detecting
		#--------------------------------------------------------------------#
		#Phase 1: burn-in; estimating parameters
		
		#we are only resetting the burn-in estimators
		
		#need to estimate mean and variance in burn in period
		#we set the forgetting factor to be 1
		lambda_beforeBL <- 1
		w_BL_new <- 0
		u_BL_new <- 0
		v_BL_new <- 0
		v_BL_old <- 0
		#initialize mean and variance
		mean_BL_new <- 0
		var_BL_new <- 0
		
		#set the estimated AFF deriv to 0
		AFFderiv_estim_BL <- 0
		
		#reset BLcount
		BLcount <- 0
		#cat("streampos: ", streampos,  ", AFFderiv_estim_BL: ", AFFderiv_estim_BL, "\n")
		
		#set values for BL mean and variance
		mean_burn <- mu0
		var_burn <- sigma0^2 
		burnvec <- c(p, mean_burn, var_burn)
		#end of Phase 1
		#--------------------------------------------------------------------#
		
		
        if (var_burn > 0){
            inv_var_burn <- 1 / var_burn
        } else {
            inv_var_burn <- 1
        }
		
		#--------------------------------------------------------------------#
		#Phase 2: after BL; monitoring for jump
		
		#jump_found is a boolean that flags if a jump is detected
		#jump_found <- FALSE
		
		thelimits <- c(0,0)
		
		#a boolean saying if the jump has been detected
		jumpdetected <- FALSE
		while((jumpdetected==FALSE) && (streampos < N)){
			#get the next observation from stream
			streampos <- update_streampos(streampos)
			x_new <- get_nextobs_fromstream(stream, streampos)
			#cat("streampos: ", streampos, ", lambda: ", lambda, "\n")
			

			#derivatives
			Delta_new <- update_Delta(Delta_new, lambda, m_new)
			Omega_new <- update_Omega(Omega_new, lambda, w_new)

			#easier to update weight than to recalculate it
			#get BL+1 w and u
			m_new <- update_m(m_new, lambda, x_new)
			
			w_new <- update_w(w_new, lambda)
			u_new <- update_u(u_new, w_new)
			#no need to update v, but we do anyway...
			v_old <- v_new
			v_new <- getv_from_u_and_w(u_new, w_new)
			
			#update mean
			affmean_old <- affmean_new
			affmean_new <- get_xbar(m_new, w_new)
			

			xbar_old_deriv <- xbar_new_deriv
			xbar_new_deriv <- get_xbar_deriv(Delta_new, affmean_new, Omega_new, w_new)
			
			#update variance
			affvar_old <- affvar_new
			affvar_new <- update_var(affvar_old, affmean_old, lambda, v_old, v_new, w_new, x_new)
			
			#new function used to calculate limits efficiently, in one function, instead of in several
			thelimits <- getNormalLimitsFromList(burnvec, u_new)
			
			#boolean of whether or not ffmean is int he limits
			inLimitsBool <- isInLimitsNormal(affmean_new, thelimits)
			#check for jump
			if(inLimitsBool == FALSE)
			{
				jumpdetected <- TRUE
				detected_count <- detected_count + 1
				detect_pos_vec[detected_count] <- streampos
			}
			#updating lambda
			L_deriv_new <- get_L_deriv(affmean_old, xbar_old_deriv, x_new)
			

			
			#time to scale it:
			L_deriv_scaled <- L_deriv_new*inv_var_burn
			
			lambda <- update_lambda(lambda, signchosen, alpha, L_deriv_scaled)
			lambda <- inbounds(lambda, low_bound, up_bound)


            #save!
            lambda_vec[streampos]  <- lambda
			


		} #end of while (jumpdetected==FALSE) - will restart process, if (streampos < N) 
		
	} #end while (streampos < N)
	#trim detect_pos_vec
	detect_pos_vec <- detect_pos_vec[1:detected_count]
	return( list(tau=detect_pos_vec, lambda=lambda_vec) )
}
#end AFF
