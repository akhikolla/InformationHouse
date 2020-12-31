#ifndef __global_h__
#define __global_h__

/* Global parameters that passing values for Integrand functions */
/* The reason why using global parameter is that the passing of 
   parameter for the integrand function in this version of Cuba 
   is still not avaiable.*/
static double    *global_alpha_v_,
				 *global_log_alpha_v_,
				 *global_beta_,
                 *global_sigma2_;

static int       *global_Dtime_,
                 *global_Delta_,
                 *global_G_;

				 
#endif
