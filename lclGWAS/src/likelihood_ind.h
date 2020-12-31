/**********************************/
/* Jiaxing Lin                    */
/* 11-29-2016                     */
/* likelihood for individual      */
/**********************************/


#include <math.h>
#include "global.h"

using namespace std;

/* Likelihood for patient i */
inline double likelihood_i(double log_alpha_v[],
                	         int 	dtime, 		
							 int 	delta, 
                      		 double exp_part,	
							 double prob_dtime)
{
    //adjust index for R to C++ for dtime.
    dtime -=1;
    // For the save of computational cost, the exp(\sum(log(alpha)),exp_part)
    // instead of \product(pow(alpha, exp_part)) is adapted. 
    // Since exp is cheaper than pow and the log(alpha) can be precomputed 
    // and shared as globol variables.
    double log_product = 0;

    if(dtime >=1){
        for(int j = 0; j <= dtime-1; j++)
            log_product += log_alpha_v[j];
        log_product = exp(log_product*exp_part);
        return log_product-delta*log_product*prob_dtime;
    }
    else
        return 1-delta*prob_dtime;
}


/* Integrand function for the trio case.*/ 
/* provided by Tracy*/
/* modified by Jiaxing*/
/* the integration range is set to be (0,0,0) to (1,1,1)*/

static double likelihood_ind() 
{
	double prod_l_i 		= 1.0;
	
	// precompute some qualities to save repeating calculating in
	double exp_part 		= 0.0;
	double prob_dtime 		= 0.0; 
    	
    int k = 0;
	// precompute the exp(beta*G+b)
	exp_part 	= exp( global_beta_[0]
				 * global_G_[k]);
	
		// Precompute the survival probability at dtime.last.
		// Use exp(log(alpha)*exp_part) rather than ppw(alpha,exp_part)
		// to save cpu time cost.(This is significant 30% cost reduction).
		if( global_log_alpha_v_[global_Dtime_[k]-1] == -INFINITY ||
			global_log_alpha_v_[global_Dtime_[k]-1] == INFINITY)
        	prob_dtime = 0;
    	else
    		prob_dtime = exp(global_log_alpha_v_[global_Dtime_[k]-1]
							 *exp_part);	
	
		// Compute the likelihood function for each person and mutiple up.
		// exp_part == INFINITY means extreme small likelihood.
		// Thus,it is safe to set the corresponding likelihood to zero.
		if(exp_part == INFINITY)
			prod_l_i *= 0;
		else 		
			prod_l_i = prod_l_i
                      * likelihood_i(global_log_alpha_v_,
                    	             global_Dtime_[k],	
									 global_Delta_[k],
									 exp_part,
									 prob_dtime);
                   
   	
	return prod_l_i;
}

