/*
	This is the routine calling the 
	Cuba Cuhre doing the integration
	for the likelihood function for
	mutiple trio families.
	
	Jiaxing 20160312

	Note: eventually, we will remove 
	those temp input argument dt[]
	Delta[] and G[].
*/

#include <math.h>
#include "cuba.h"
#include "integrand.h"
#include "likelihood_ind.h"
#include "global.h"
#include <iostream>

#define NCOMP 1
#define PI  M_PI

double fam_LLVar(double sigma2,
			  int fam_size[],
              int    dt[],
              int Delta[],
              int     G[],
 		      int     m )
{
  
  //need to get beta to the global here.
  global_sigma2_ = &sigma2;
  //current pid index
  int curPid = 0;
	
  int  nregions, neval, fail;
  double integral[NCOMP],            /* output argument for calling Cuhre*/
         error[NCOMP],
         prob[NCOMP];
  
  double sum_LL = 0.00;
  for(int k = 0; k < m; k++){
    for(int j = 0; j < fam_size[k]; j++){
      global_Dtime_[j] = dt[j+curPid];
       				     /* 
					grep dtime, G and Delta for each 
					family for the evaluation of the 
					integration.
 				     */
      global_G_[j]     = G[j+curPid];
      global_Delta_[j] = Delta[j+curPid];
    }

	curPid +=fam_size[k]; // update the current pid.
    
    if(fam_size[k] !=1){
    Cuhre(fam_size[k], NCOMP,
      Integrand,
      1e-3, 1e-12,
      0, 0, 50000,
      0,
      &nregions, &neval, &fail,
      integral, error, prob
      );

    sum_LL += log((double)integral[0]/(*global_sigma2_)/sqrt((*global_sigma2_)));
    }
    else
	   sum_LL += log(likelihood_ind());    
 }
 
  return -sum_LL; 		// return the negative log likelihood 
}


