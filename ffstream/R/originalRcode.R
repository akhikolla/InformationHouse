#' CUSUM change detection for a stream in R 
#'
#' Original implementation in R of CUSUM change detector, now with
#' documentation.
#'
#'
#' @param stream The stream of observations.
#'
#' @param BL The burn-in length, used to estimate the mean and variance.
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
CUSUM_stream_jumpdetect <- function(stream, BL, params){

	d <- params[[1]]
	B <- params[[2]]
	
	#will run until stream finished...
	streampos <- 1
	N <- length(stream)
	detected_count <- 0
	#vector for saving jumps, will trim later
	#just a quick way of saving space - can't have more than BL*det_count obs in stream,
	#since every det restarts BL
	M <- trunc(N/BL) + 1
	detect_pos_vec <- rep(0, M)
	
	while (streampos < N){
		#begin the process of burning in and then detecting
		#--------------------------------------------------------------------#
		#Phase 1: burn-in; estimating parameters
		
		#need to estimate mean and variance in burn in period
		#we set the forgetting factor to be 1
		lambda <- 1
		w_new <- 0
		u_new <- 0
		v_new <- 1
		v_old <- v_new
		#initialize mean and variance
		mean_new <- 0
		var_new <- 0
		
		for (k in 1:BL){
			#get the next observation from stream
			x_new <- get_nextobs_fromstream(stream, streampos)
			streampos <- update_streampos(streampos)
			
			#easier to update weight than to recalculate it
			#get BL+1 w and u
			w_new <- update_w(w_new, lambda)
			u_new <- update_u(u_new, w_new)
			#update v
			v_old <- v_new
			v_new <- getv_from_u_and_w(u_new, w_new)
			
			#update mean; old mean saved for variance later
			#alternatively, could update variance first
			#just the way the derivations work out...(deeper meaning, actually...)
			mean_old <- mean_new
			mean_new <- update_mean(mean_old, w_new, x_new)
			
			#update variance
			if (w_new==1){
				#this is the case k=1, to elimate repeated code, we introduce the if statement
				var_new <- 0
			} else {
				var_old <- var_new
				var_new <- update_var(var_old, mean_old, lambda, v_old, v_new, w_new, x_new)
			}

		} #end of for BL
		
		#set values for CUSUM mean and variance
		mean_j <- mean_new
		var_j <- var_new
		sd_j <- sqrt(var_j)
		
		nu_j <- d*sd_j
		control_j <- B*sd_j
		
		#----------end of Phase 1: burn-in----------#
		#--------------------------------------------------------------------#
		#Phase 2: detect change
		#CUSUM starts, sample from distribution 1
		S_j <- 0
		T_j <- 0
		#jump_found is a boolean that flags if a jump is detected
		jump_found <- FALSE
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







#' EWMA change detection for a stream in R with 
#'
#' Original implementation in R of EWMA change detector, now with
#' documentation.
#'
#'
#' @param stream The stream of observations.
#'
#' @param BL The burn-in length, used to estimate the mean and variance.
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
#'
#' @return A vector of estimated changepoints.
#'
#'
#' @keywords internal
EWMA_stream_jumpdetect <- function(stream, BL, params){

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
	M <- trunc(N/BL) + 1
	detect_pos_vec <- rep(0, M)
	
	delta <- r / (2-r)
    rFactorSigmaZ <- 1
    sigmaZ <- 1

	while (streampos < N){
		#begin the process of burning in and then detecting
		#--------------------------------------------------------------------#
		#Phase 1: burn-in; estimating parameters
		
		#need to have estimates of parameters before can run ewma
		#initialize W
		#don't know what mu_1 is
		W <- 0
		
		#need to estimate mean and variance in burn in period
		#we set the forgetting factor to be 1
		lambda <- 1
		w_new <- 0
		u_new <- 0
		v_new <- 1
		v_old <- v_new
		#initialize mean and variance
		mean_new <- 0
		var_new <- 0

		
		for (k in 1:BL){
			#get the next observation from stream
			x_new <- get_nextobs_fromstream(stream, streampos)
			streampos <- update_streampos(streampos)
			
			#the ewma updating step
			W <- r * x_new + (1-r)*W
			
			#easier to update weight than to recalculate it
			#get BL+1 w and u
			w_new <- update_w(w_new, lambda)
			u_new <- update_u(u_new, w_new)
			#update v
			v_old <- v_new
			v_new <- getv_from_u_and_w(u_new, w_new)
			
			#update mean; old mean saved for variance later
			#alternatively, could update variance first
			#just the way the derivations work out...(deeper meaning, actually...)
			mean_old <- mean_new
			mean_new <- update_mean(mean_old, w_new, x_new)
			
			#update variance
			if (w_new==1){
				#this is the case k=1, to elimate repeated code, we introduce the if statement
				var_new <- 0
			} else {
				var_old <- var_new
				var_new <- update_var(var_old, mean_old, lambda, v_old, v_new, w_new, x_new)
			}

		} #end of for BL
		
		#set values for EWMA mean and variance
		mu_1 <- mean_new
		sigma1_sq <- var_new
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
        W <- mean_new
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
#' Update the stream position
#'
#' Updates the index keeping track of the stream position
#'
#' @param streampos The index of the stream position.
#'
#' @keywords internal
update_streampos <- function(streampos){
	return(streampos+1)
}


#-------------------------------------------------------------------------#
#' Obtain the next observation 
#'
#' Obtain the next observation from the stream 
#' 
#' @param stream The stream (vector) of observations.
#'
#' @param streampos The index of the stream position.
#'
#' @keywords internal
get_nextobs_fromstream <- function(stream, streampos){
	return(stream[streampos])
}




#-------------------------------------------------------------------------#
#' Change detection using the Fixed Forgetting Factor method
#'
#' Original implementation in R of FFF change detector, but now
#'
#'
#' @param stream The stream of observations.
#'
#' @param BL The burn-in length, used to estimate the mean and variance.
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
#'
#' @return A vector of estimated changepoints.
#'
#'
#' @keywords internal
FFF_stream_jumpdetect <- function(stream, BL, ffparams){

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
	#vector for saving jumps, will trim later
	#just a quick way of saving space - can't have more than BL*det_count obs in stream,
	#since every det restarts BL
	M <- trunc(N/BL) + 1
	detect_pos_vec <- rep(0, M)
	
	while (streampos < N){
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
		
		#update two sets of means; one for params (mean_new, w_BL_new)
		#and one for ff (ffmean_new), along with u, v, w, and var
		for (k in 1:BL){
			#get the next observation from stream
			x_new <- get_nextobs_fromstream(stream, streampos)
			streampos <- update_streampos(streampos)
			
			#easier to update weight than to recalculate it
			#get BL+1 w and u
			w_BL_new <- update_w(w_BL_new, lambda_beforeBL)
			w_new <- update_w(w_new, lambda)
			
			u_BL_new <- update_u(u_BL_new, w_BL_new)
			u_new <- update_u(u_new, w_new)
			
			#update v
			v_BL_old <- v_BL_new
			v_BL_new <- getv_from_u_and_w(u_BL_new, w_BL_new)
			v_old <- v_new
			v_new <- getv_from_u_and_w(u_new, w_new)
			
			#update mean; old mean saved for variance later
			#alternatively, could update variance first
			#just the way the derivations work out...(deeper meaning, actually...)
			mean_BL_old <- mean_BL_new
			mean_BL_new <- update_mean(mean_BL_old, w_BL_new, x_new)
			ffmean_old <- ffmean_new
			ffmean_new <- update_ffmean(ffmean_old, w_new, x_new)
			
			#update variance
			if (w_BL_new==1){
				#this is the case k=1, to elimate repeated code, we introduce the if statement
				var_BL_new <- 0
			} else {
				var_BL_old <- var_BL_new
				var_BL_new <- update_var(var_BL_old, mean_BL_old, lambda_beforeBL, v_BL_old, v_BL_new, w_BL_new, x_new)
			}
			
			#update variance
			if (w_new==1){
				#this is the case k=1, to elimate repeated code, we introduce the if statement
				ffvar_new <- 0
			} else {
				ffvar_old <- ffvar_new
				ffvar_new <- update_var(ffvar_old, ffmean_old, lambda, v_old, v_new, w_new, x_new)
			}
			
			
		} #end of for BL
		
		#set values for BL mean and variance
		mean_burn <- mean_BL_new
		var_burn <- var_BL_new
		burnvec <- c(p, mean_burn, var_burn)
		
		
		#--------------------------------------------------------------------#
		#Phase 2: between burn-in and change
		if (resettozero==0){
			w_new <- 0
			u_new <- 0
			v_new <- 0
			v_old <- 0
			
			#initialize mean and variance
			ffmean_new <- mean_BL_new
			ffvar_new <- var_BL_new
			
			#cat("resetting ff arl1...\n")
			
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
				jumpdetected <- TRUE
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
#' Check if in the limits
#'
#' Checks if a value is in normal limits
#'
#' @param x The value to check.
#'
#' @param limitsvec The vector of limits (two values).
#'
#' @keywords internal
isInLimitsNormal <- function(x, limitsvec){
	tempbool <- TRUE
	if ( (x < limitsvec[1]) || (x > limitsvec[2])){
		tempbool <- FALSE
	}
	return(tempbool)
}


#--------------------------isInLimits--------------------------#
#isInLimits <- function(x, limitsvec){
#	tempbool <- TRUE
#	if ( (x < limitsvec[1]) || (x > limitsvec[2])){
#		tempbool <- FALSE
#	}
#	return(tempbool)
#}


#----------------update_ffvar-------------------------------#
#used to be update_ffvar; removed ff
#-------------------------------------------------------------------------#
#' Update the FF variance
#'
#' Function to sequentially update the forgetting factor variance
#'
#' @param ffvar The new (current) value of the forgetting factor variance.
#'
#' @param ffmean_old The old (previous) value of the forgetting factor mean.
#'                   In fact, the formula requires this value, rather than
#'                   the current value of the FF mean.
#'
#' @param lambda The value of the forgetting factor.
#'
#' @param v_old The old (previous) value of \code{v}.
#'
#' @param v_new The new (current) value of \code{v}.
#'
#' @param w_new The new (current) value of \code{w}.
#'
#' @param x_new The new (current) observation.
#'
#'
#' @keywords internal
update_var <- function(ffvar, ffmean_old, lambda, v_old, v_new, w_new, x_new){
	frac <- (w_new - 1)/w_new
	temp <- lambda*v_old*ffvar + frac * (ffmean_old-x_new)^2
	temp <- temp/v_new
	return(temp)
}

#-------------------------------------------------------------------------#
#' Update the weight
#'
#' Function to sequentially update the weight used in the computation of
#' the forgetting factor mean
#'
#'
#' @param w_old The old (previous) value of \code{w}.
#'
#' @param lambda_old The value of the forgetting factor. No real need to
#'                   add 'old' here, in the adaptive scheme, these are the
#'                   the values from the previous iteration.
#'
#' @keywords internal
update_w <- function(w_old, lambda_old){
	return(lambda_old * w_old + 1)
}



#-------------------------------------------------------------------------#
#' Update \code{u}
#'
#' Function to sequentially update the weight used to scale the forgetting
#' factor variance, \code{u}.
#'
#'
#' @param u_old The old (previous) value of \code{u}.
#'
#' @param w_new The new (current) value of \code{w}.
#'
#'
#' @keywords internal
update_u <- function(u_old, w_new){
	#(1- 1/w_N+1)^2 * u_N + (1/w_N+1)^2
	a <- (w_new-1)/w_new
	b <- 1/w_new
	return(a^2*u_old + b^2)
}



#-------------------------------------------------------------------------#
#' Calculate \code{v}
#'
#' Function to calculate \code{v} from \code{u} and \code{w}, according 
#' to the formula \eqn{v= w(1-u)}.
#'
#'
#' @param u The current value of \code{u}.
#'
#' @param w The current value of \code{w}.
#'
#'
#' @keywords internal
getv_from_u_and_w <- function(u, w){
	return( w*(1-u) )
}


#-------------------------------------------------------------------------#
#' Calculate \code{m}
#'
#' Function to calculate \code{m} 
#'
#'
#' @param m_old The old (previous) value of \code{m}.
#'
#' @param lambda_old The old (previous) value of \code{lambda}.
#'
#' @param x_new The current observation.
#'
#'
#' @keywords internal
update_m <- function(m_old, lambda_old, x_new){
	return(lambda_old * m_old + x_new)
}


#-------------------------------------------------------------------------#
#' Compute \code{xbar} 
#'
#' Function to calculate \code{xbar} from \code{m} and \code{w}
#'
#'
#' @param m The current value of \code{m}.
#'
#' @param w The current value of \code{w}.
#'
#'
#' @keywords internal
get_xbar <- function(m, w){
	return(m/w)
}

#-------------------------------------------------------------------------#
#' Update the FF mean
#'
#' Function to sequentially update the FF mean
#'
#'
#' @param ffmean_old The old (previous) value of \code{ffmean}.
#'
#' @param w_new The new (current) value of \code{w}.
#'
#' @param x_new The current observation.
#'
#'
#' @keywords internal
update_mean <- function(ffmean_old, w_new, x_new){
	temp <- ( (w_new-1)*ffmean_old + x_new )/w_new
	return(temp)
}


#-------------------------------------------------------------------------#
#' Update the FF mean (duplicate)
#'
#' Function to sequentially update the FF mean
#'
#'
#' @param ffmean_old The old (previous) value of \code{ffmean}.
#'
#' @param w_new The new (current) value of \code{w}.
#'
#' @param x_new The current observation.
#'
#' @seealso \code{\link{update_mean}}
#'
#'
#' @keywords internal
update_ffmean <- function(ffmean_old, w_new, x_new){
	temp <- ( (w_new-1)*ffmean_old + x_new )/w_new
	return(temp)
}



#-------------------------------------------------------------------------#
#' Update \code{Delta}
#'
#' Function to sequentially update the value of \code{Delta}
#'
#'
#' @param Delta_old The old (previous) value of \code{Delta}.
#'
#' @param lambda_old The value of the forgetting factor. 
#'
#' @param m_new The new (current) value of \code{m}.
#'
#'
#' @keywords internal
update_Delta <- function(Delta_old, lambda_old, m_new){
	return(lambda_old*Delta_old + m_new)
}


#-------------------------------------------------------------------------#
#' Update \code{Omega}
#'
#' Function to sequentially update the value of \code{Omega}
#'
#'
#' @param Omega_old The old (previous) value of \code{Omega}.
#'
#' @param lambda_old The value of the forgetting factor.
#'
#' @param m_new The new (current) value of \code{m}.
#'
#'
#' @keywords internal
update_Omega <- function(Omega_old, lambda_old, w_new){
	return(lambda_old*Omega_old + w_new)
}


#-------------------------------------------------------------------------#
#' Compute the derivative of \code{xbar}
#'
#' Function to compute the derivative of \code{xbar}.
#'
#'
#' @param Delta The value of \code{Delta}, the derivative of \code{m}.
#'
#' @param xbar The value of the forgetting factor mean.
#'
#' @param Omega The value of \code{Omega}, the derivative of \code{w}.
#'
#' @param w The value of the weight used to compute \code{xbar}.
#'
#'
#' @keywords internal
get_xbar_deriv <- function(Delta, xbar, Omega, w){
	return( (Delta - xbar*Omega)/w )
}






#-------------------------------------------------------------------------#
#' Check if a value is inside bounds
#'
#' Function to check that a value \code{x} is between two values
#'
#'
#' @param x The value to check.
#'
#' @param low_bound The lower bound.
#'
#' @param up_bound The upper bound.
#'
#'
#' @keywords internal
inbounds <- function(x, low_bound, up_bound){
	temp <- x
	if(x < low_bound){
		temp <- low_bound
	} 
	else if (x > up_bound){
		temp <- up_bound
	}
	return(temp)
	
}


#-------------------------------------------------------------------------#
#' Update lambda
#'
#' Update the forgetting factor \code{lambda} using the derivative of
#' a particular cost function
#'
#'
#' @param lambda_old The old (previous) value of \code{lambda}.
#'
#' @param signchosen Either \code{+1} or \code{-1}. Usually \code{-1}.
#'
#' @param alpha The value of the step-size in the gradient descent.
#'              In the paper this is named \eqn{\epsilon}.
#'
#' @param deriv_new The value of the derivative at the current time.
#'
#'
#' @keywords internal
update_lambda <- function(lambda_old, signchosen, alpha, deriv_new){
	return(lambda_old + signchosen*alpha*deriv_new)
}




#-------------------------------------------------------------------------#
#' Compute the derivative of L
#'
#' Compute the derivative of the cost function L, defined as:
#' \deqn{L_{N} = [\bar{x}_{N-1, \overrightarrow{\lambda}} - x_{N}}
#'
#'
#' @param xbar_old The old (previous) value of \code{xbar}.
#'
#' @param xbar_old_deriv The old (previous) value of the derivative
#'                       of \code{xbar}.
#'
#' @param x_new The value of the current observation. 
#'
#'
#' @keywords internal
get_L_deriv <- function(xbar_old, xbar_old_deriv, x_new){
	return ( -2*xbar_old_deriv*(x_new-xbar_old) )
}



#-------------------------------------------------------------------------#
#' Compute confidence limits assuming the normal distribution
#'
#' Gets the 100p% confidence interval limits for N(mu, stddev).
#' NB: stddev used, not variance!
#'
#' @param burnvec A vector of \emph{unnamed} values, consisting of:
#'                \describe{
#'                          \item{p}{First component - value of threshold.}
#'
#'                          \item{mu}{Second component - value of the mean.}
#'
#'                          \item{sd}{Third component - value of the standardi
#'                                     deviation.}
#'
#'                         }
#'
#' @param u_new The value of \code{u}, the quantity used to scale the 
#'              forgetting factor variance.
#'
#'
#' @keywords internal
getNormalLimitsFromList <- function(burnvec, u_new){
	p <- burnvec[1]
	mu <- burnvec[2]
	stddev <- sqrt(u_new*burnvec[3])
	
	p1 <- (1-p)/2
	p2 <- p+p1 
	q1 <- qnorm(p1, mean = mu, sd=stddev)
	q2 <- qnorm(p2, mean = mu, sd=stddev)
	return(c(q1, q2))	
}


# #----------------------getmeanchebylimits----------------------#
# getmeanchebylimits <- function(mu, sigma_sq_u, p){
#     #cheby eqn:
#     # Pr(|X-mu| >= k*sigma) <= 1/(k^2)
#     #i.e. Pr(|xbar-mu| >= k*sigma) <= 1/(k^2)
#     sigma_sqrt_u <- sqrt(sigma_sq_u)
#     k <- sqrt(2/(1-p))
#     cwidth <- k*sigma_sqrt_u
#     a <- mu - cwidth
#     b <- mu + cwidth
#     thelimits <- c(a,b)
#     return(thelimits)
# }






#-------------------------------------------------------------------------#
#' Compute the adaptive forgetting factor - no change detection
#'
#' Original implementation in R of the AFF
#'
#'
#' @param stream The stream of observations.
#'
#' @param BL The burn-in length, used to estimate the mean and variance.
#'
#' @param params An \emph{unnamed} list of parameters for the FFF algorithm. 
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
#'
#' @return A vector with the values of the adaptive forgetting factor
#'         \eqn{\overrightarrow{\lambda}}.
#'
#'
#' @keywords internal
#no change detection, and returns vector of lambdao values
AFF_scaled_stream_no_change_detection <- function(stream, BL, params){

	#parameters needed to run AFF algorithm
	lambda_init <- 1
	#alpha <- 0.01
	low_bound <- 0.6
	up_bound <- 1
	signchosen <- -1
	
	lambda <- params[[1]]
	p <- params[[2]]
	resettozero <- params[[3]]
	#cat("okay...")
	u_new <- params[[4]]
	v_new <- params[[5]]
	v_old <- v_new
	w_new <- params[[6]]
	affmean_new <- params[[7]]
	#need m now...
	m_new <- affmean_new*w_new 
	affvar_new <- params[[8]]
	low_bound <- params[[9]]
	up_bound  <- params[[10]]
	#should be minus 1, I think
	signchosen <- params[[11]]
	alpha <- params[[12]]
	#init this somehow?
	xbar_new_deriv <- 0
	Delta_new <- 0
	Omega_new <- 0
	Delta_old <- 0
	Omega_old <- 0
	
	#will run until stream finished...
	streampos <- 0
	N <- length(stream)

    lambda_vec <- rep(0, N)

    xbar_old_deriv <- xbar_new_deriv

	detected_count <- 0
	#vector for saving jumps, will trim later
	#just a quick way of saving space - can't have more than BL*det_count obs in stream,
	#since every det restarts BL
	
	mean_burn <- 0
	var_burn <- 1
	
	
	#this is the quantity by which we multiply L_deriv_new - the inverse of the affderiv average
	inv_AFFderiv_estim_BL <- 0
	inv_var_burn <- 1
	
	#a vector for saving AFF
	#we deleted the saving of the adlambdavec - see testaff1.R for a version with it saved
	while (streampos < N){
		#begin the process of burning in and then detecting
		#--------------------------------------------------------------------#
		#Phase 1: burn-in; estimating parameters
		
		#we are only resetting the burn-in estimators
		
		#need to estimate mean and variance in burn in period
		#we set the forgetting factor to be 1
		lambda_beforeBL <- 1

        #lambda_beforeBL <- 0.95
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
		#update two sets of means; one for params (mean_new, w_BL_new)
		#and one for ff (ffmean_new), along with u, v, w, and var
		#this while loop is very important!
		while ((BLcount < BL) && (streampos < N )){
			#increase BLcount, before we forget...
			BLcount <- BLcount + 1
			#cat("streampos: ", streampos, ", lambda: ", lambda, "\n")
			
			#get the next observation from stream
			streampos <- update_streampos(streampos)
			x_new <- get_nextobs_fromstream(stream, streampos)
            
            Delta_old <- Delta_new
			Delta_new <- update_Delta(Delta_new, lambda, m_new)
            Omega_old <- Omega_new
			Omega_new <- update_Omega(Omega_new, lambda, w_new)
			
			
			#easier to update weight than to recalculate it
			#get BL+1 w and u
			#cat("m_new: ", m_new, ", x_new: ", x_new, "\n")
			m_new <- update_m(m_new, lambda, x_new)
			
			w_BL_new <- update_w(w_BL_new, lambda_beforeBL)
			w_new <- update_w(w_new, lambda)
			
			u_BL_new <- update_u(u_BL_new, w_BL_new)
			u_new <- update_u(u_new, w_new)
			
			#update v
			v_BL_old <- v_BL_new
			v_BL_new <- getv_from_u_and_w(u_BL_new, w_BL_new)
			v_old <- v_new
			v_new <- getv_from_u_and_w(u_new, w_new)
			
			#update mean; old mean saved for variance later
			#alternatively, could update variance first
			#just the way the derivations work out...(deeper meaning, actually...)
			mean_BL_old <- mean_BL_new
			mean_BL_new <- update_mean(mean_BL_old, w_BL_new, x_new)
			
			#update mean
			affmean_old <- affmean_new
			affmean_new <- get_xbar(m_new, w_new)
			
			xbar_old_deriv <- xbar_new_deriv
			xbar_new_deriv <- get_xbar_deriv(Delta_new, affmean_new, Omega_new, w_new)
			
			
			#update variance
			if (w_BL_new==1){
				#this is the case k=1, to elimate repeated code, we introduce the if statement
				var_BL_new <- 0
			} else {
				var_BL_old <- var_BL_new
				var_BL_new <- update_var(var_BL_old, mean_BL_old, lambda_beforeBL, v_BL_old, v_BL_new, w_BL_new, x_new)
			}
			
			#update variance
			if (w_new==1){
				#this is the case k=1, to elimate repeated code, we introduce the if statement
				affvar_new <- 0
			} else {
				affvar_old <- affvar_new
				affvar_new <- update_var(affvar_old, affmean_old, lambda, v_old, v_new, w_new, x_new)
			}
			
			#we are updating lambda here (or not)
			#---------------------#
			#updating lambda
			L_deriv_new <- get_L_deriv(affmean_old, xbar_old_deriv, x_new)
			
			#time to scale it:
			L_deriv_scaled <- L_deriv_new*inv_var_burn
            #old: 
#			L_deriv_scaled <- L_deriv_new*inv_AFFderiv_estim_BL
			
			#lambda <- update_lambda(lambda, signchosen, alpha, L_deriv_scaled)
			#lambda <- inbounds(lambda, low_bound, up_bound)
			#---------------------#
			
			
			#updates with abs value of L_deriv. Later used to help scale derivative.
			AFFderiv_estim_BL <- update_mean(AFFderiv_estim_BL, w_BL_new, abs(L_deriv_new ))
			
            lambda_vec[streampos]  <- lambda
			
		} #end of for BL
		#cat("streampos: ", streampos,  ", AFFderiv_estim_BL: ", AFFderiv_estim_BL, "\n")
		
		#set values for BL mean and variance
		mean_burn <- mean_BL_new
		var_burn <- var_BL_new
		burnvec <- c(p, mean_burn, var_burn)
		#end of Phase 1
		#--------------------------------------------------------------------#
		
		
		#reinitialise this quantity after burn-in
		inv_AFFderiv_estim_BL <- 1/AFFderiv_estim_BL 
		
		if (var_burn>0){	
			inv_var_burn <- 1/var_burn
		}
		
		#--------------------------------------------------------------------#
		#Phase 2: after BL; monitoring for jump
		
		#jump_found is a boolean that flags if a jump is detected
		#jump_found <- FALSE
		
		thelimits <- c(0,0)
		
		#a boolean saying if the jump has been detected
		jumpdetected <- FALSE
		while((jumpdetected==FALSE) && (streampos < N)  ){
			#get the next observation from stream
			streampos <- update_streampos(streampos)
			x_new <- get_nextobs_fromstream(stream, streampos)


			#derivatives
            Delta_old <- Delta_new
			Delta_new <- update_Delta(Delta_new, lambda, m_new)
            Omega_old <- Omega_new
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
			#no checking for jump!
			
			#updating lambda
			L_deriv_new <- get_L_deriv(affmean_old, xbar_old_deriv, x_new)
			
			#time to scale it:
			#L_deriv_scaled <- L_deriv_new*inv_AFFderiv_estim_BL
			L_deriv_scaled <- L_deriv_new*inv_var_burn
			
			lambda <- update_lambda(lambda, signchosen, alpha, L_deriv_scaled)
			lambda <- inbounds(lambda, low_bound, up_bound)
			
            #save!
            lambda_vec[streampos]  <- lambda

		} #end of while (jumpdetected==FALSE) - will restart process, if (streampos < N) 
		
		#end of Phase 2 - now either restart or if stream is finsihed return the detect_pos_vec
		
	} #end while (streampos < N)
	
	return(lambda_vec)		
}






#-------------------------------------------------------------------------#
#' Change detection using the AFF method
#'
#' Original implementation in R of the AFF
#'
#'
#' @param stream The stream of observations.
#'
#' @param BL The burn-in length, used to estimate the mean and variance.
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
#'
#' @return A vector with the values of the adaptive forgetting factor
#'         \eqn{\overrightarrow{\lambda}}.
#'
#'
#' @keywords internal
AFF_scaled_stream_jumpdetect <- function(stream, BL, affparams){
	#parameters needed to run AFF algorithm
	lambda_init <- 1
	#alpha <- 0.01
	low_bound <- 0.6
	up_bound <- 1
	signchosen <- -1
	
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
	M <- trunc(N/BL) + 1
	detect_pos_vec <- rep(0, M)
	
	#this is the quantity by which we multiply L_deriv_new - the inverse of the affderiv average
	inv_AFFderiv_estim_BL <- 0
    inv_var_burn <- 0
    lambda_vec <- rep(0, N)
	
	#a vector for saving AFF
	#we deleted the saving of the adlambdavec - see testaff1.R for a version with it saved
	while (streampos < N){
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
		#update two sets of means; one for params (mean_new, w_BL_new)
		#and one for ff (ffmean_new), along with u, v, w, and var
		#this while loop is very important!
		while ((BLcount < BL) && (streampos < N)){
			#increase BLcount, before we forget...
			BLcount <- BLcount + 1
			#cat("streampos: ", streampos, ", lambda: ", lambda, "\n")
			
			#get the next observation from stream
			streampos <- update_streampos(streampos)
			x_new <- get_nextobs_fromstream(stream, streampos)
            

			#derivatives
			Delta_new <- update_Delta(Delta_new, lambda, m_new)
			Omega_new <- update_Omega(Omega_new, lambda, w_new)
			
			
			#easier to update weight than to recalculate it
			#get BL+1 w and u
			#cat("m_new: ", m_new, ", x_new: ", x_new, "\n")
			m_new <- update_m(m_new, lambda, x_new)
			
			w_BL_new <- update_w(w_BL_new, lambda_beforeBL)
			w_new <- update_w(w_new, lambda)
			
			u_BL_new <- update_u(u_BL_new, w_BL_new)
			u_new <- update_u(u_new, w_new)
			
			#update v
			v_BL_old <- v_BL_new
			v_BL_new <- getv_from_u_and_w(u_BL_new, w_BL_new)
			v_old <- v_new
			v_new <- getv_from_u_and_w(u_new, w_new)
			
			#update mean; old mean saved for variance later
			#alternatively, could update variance first
			#just the way the derivations work out...(deeper meaning, actually...)
			mean_BL_old <- mean_BL_new
			mean_BL_new <- update_mean(mean_BL_old, w_BL_new, x_new)
			
			#update mean
			affmean_old <- affmean_new
			affmean_new <- get_xbar(m_new, w_new)
			
			xbar_old_deriv <- xbar_new_deriv
			xbar_new_deriv <- get_xbar_deriv(Delta_new, affmean_new, Omega_new, w_new)
			
			
			#update variance
			if (w_BL_new==1){
				#this is the case k=1, to elimate repeated code, we introduce the if statement
				var_BL_new <- 0
			} else {
				var_BL_old <- var_BL_new
				var_BL_new <- update_var(var_BL_old, mean_BL_old, lambda_beforeBL, v_BL_old, v_BL_new, w_BL_new, x_new)
			}
			
			#update variance
			if (w_new==1){
				#this is the case k=1, to elimate repeated code, we introduce the if statement
				affvar_new <- 0
			} else {
				affvar_old <- affvar_new
				affvar_new <- update_var(affvar_old, affmean_old, lambda, v_old, v_new, w_new, x_new)
			}
			
			#we are updating lambda here (or not)
			#---------------------#
			#updating lambda
			L_deriv_new <- get_L_deriv(affmean_old, xbar_old_deriv, x_new)
			
			#time to scale it:
			#L_deriv_scaled <- L_deriv_new*inv_AFFderiv_estim_BL		
			
			#lambda <- update_lambda(lambda, signchosen, alpha, L_deriv_scaled)
			#lambda <- inbounds(lambda, low_bound, up_bound)
			#---------------------#


            #save!
            lambda_vec[streampos]  <- lambda
			
			
			#updates with abs value of L_deriv. Later used to help scale derivative.
			AFFderiv_estim_BL <- update_mean(AFFderiv_estim_BL, w_BL_new, abs(L_deriv_new ))


			
		} #end of for BL
		
		#set values for BL mean and variance
		mean_burn <- mean_BL_new
		var_burn <- var_BL_new
		burnvec <- c(p, mean_burn, var_burn)
		#end of Phase 1
		#--------------------------------------------------------------------#
		
		
		#reinitialise this quantity after burn-in
		inv_AFFderiv_estim_BL <- 1/AFFderiv_estim_BL 

        if (var_burn > 0){
            inv_var_burn <- 1 / var_burn
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
#                cat("\nCHANGED DETECTED: affmean=", affmean_new, ", limits=(", thelimits[1], ", ", thelimits[2], ")\n")
				jumpdetected <- TRUE
				detected_count <- detected_count + 1
				detect_pos_vec[detected_count] <- streampos
			}
			
			#updating lambda
			L_deriv_new <- get_L_deriv(affmean_old, xbar_old_deriv, x_new)
			
			#time to scale it:
#			L_deriv_scaled <- L_deriv_new*inv_AFFderiv_estim_BL		
			L_deriv_scaled <- L_deriv_new*inv_var_burn
			
			lambda <- update_lambda(lambda, signchosen, alpha, L_deriv_scaled)
			lambda <- inbounds(lambda, low_bound, up_bound)


            #save!
            lambda_vec[streampos]  <- lambda
			


		} #end of while (jumpdetected==FALSE) - will restart process, if (streampos < N) 
		
		#end of Phase 2 - now either restart or if stream is finsihed return the detect_pos_vec
		
	} #end while (streampos < N)
	
	
	#trim detect_pos_vec
	detect_pos_vec <- detect_pos_vec[1:detected_count]
	return( list(tau=detect_pos_vec, lambda=lambda_vec) )
}

