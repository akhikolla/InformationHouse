/**********************************/
/* Jiaxing Lin                    */
/* 08-29-2016                     */
/* Integrand for the likelihood   */
/* for trios families partial     */
/* derivate to sigma square       */
/**********************************/


#include <math.h>
#include "global.h"

using namespace std;

/* Likelihood for patient i */
inline double likelihood_i_c(double log_alpha_v[],
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


static void Integrand(const int    *ndim, 
					  const double xx[],
                      const int    *ncomp, 
					  double 	   ff[])
{
	int ndim_ = *ndim;

	double b_trans[4]; 		//Transform xx into b
	double txmo[4];			//2x-1	
	double usq[4]; 			//(2x-1)^2

	double prod_l_i 		= 1.0;
	
	// precompute some qualities to save repeating calculating in
	double exp_part 		= 0.0;
	double prob_dtime 		= 0.0; 
  
	for(int k = 0; k < ndim_; k++){
    	
		// compute the Jacobian element
		// This part should be modified when the member in a family 
		// is not three.
		txmo[k] 	= 2*xx[k] - 1;
		usq[k]      = txmo[k]*txmo[k];
        b_trans[k] 	= txmo[k]/(1 - usq[k]);
			
		// precompute the exp(beta*G+b)
		exp_part 	= exp( global_beta_[0]
					 * global_G_[k]
					 + b_trans[k]);
	
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
            	    * (usq[k]+1)
                	/ (usq[k]-1)/(usq[k]-1)
                      * likelihood_i_c(global_log_alpha_v_,
                    	             global_Dtime_[k],	
									 global_Delta_[k],
									 exp_part,
									 prob_dtime);
                   
	}
   	
	// Compute the transformation of pdf for changing of variables. 
	// using switch for different family size.
	double pdf_comp = 0.0;
    double bkb = 0.0; 
	
	// exponential part for distribution function regarding to different 
	// dimensions.
	
	switch(ndim_)
	{
	case 3:
	  bkb =	
	    (1/((*global_sigma2_)*2))
	  *(
	      2*b_trans[0]*(b_trans[1]+b_trans[2]- b_trans[0])
	  -b_trans[1]*b_trans[2]
	  -1.5*(b_trans[1]*b_trans[1]+b_trans[2]*b_trans[2])
	  );
	  break;
	case 2:
	  bkb = 
	    (1/((*global_sigma2_)*2))
	  *1.33334*(
	      -b_trans[0]*b_trans[0]
	      +b_trans[0]*b_trans[1]
	      -b_trans[1]*b_trans[1]
	  );
	  break;
	case 4:
	  bkb = 
	    -(1/((*global_sigma2_)*2))
	    *2*(
	        +b_trans[0]*b_trans[0]
	        +b_trans[1]*b_trans[1]
	        +b_trans[2]*b_trans[2]
	        +b_trans[3]*b_trans[3]
            -b_trans[0]*b_trans[2]
	        -b_trans[0]*b_trans[3]
	        -b_trans[1]*b_trans[2]
	        -b_trans[1]*b_trans[3]
	        +b_trans[2]*b_trans[3]
	    );
	  break;
	  
	}	
	//distribution function exponential part
	pdf_comp = exp(bkb);

	// The resulted values are store into ff vector.
	// ff[0] is the likelihood function
	
    ff[0] = prod_l_i*pdf_comp;
	return;
}

