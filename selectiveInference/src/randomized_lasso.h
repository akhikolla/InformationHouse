#ifdef __cplusplus
extern "C"
{
#endif /* __cplusplus */

double log_density_gaussian(double noise_scale,             // Scale of randomization
			    int ndim,                       // Number of features -- "p"
			    int ninternal,                  // Dimension of internal data representation often 1
			    int noptimization,              // Dimension of optimization variables -- "p"
			    double *internal_linear,        // A_D -- linear part for data
			    double *internal_state,         // D -- data state
			    double *optimization_linear,    // A_O -- linear part for optimization variables
			    double *optimization_state,     // O -- optimization state
			    double *offset);                // h -- offset in affine transform -- "p" dimensional 
  
double log_density_laplace(double noise_scale,             // Scale of randomization
			   int ndim,                       // Number of features -- "p"
			   int ninternal,                  // Dimension of internal data representation often 1
			   int noptimization,              // Dimension of optimization variables -- "p"
			   double *internal_linear,        // A_D -- linear part for data
			   double *internal_state,         // D -- data state
			   double *optimization_linear,    // A_O -- linear part for optimization variables
			   double *optimization_state,     // O -- optimization state
			   double *offset);                // h -- offset in affine transform -- "p" dimensional 

double log_density_gaussian_conditional(double noise_scale,             // Scale of randomization
					int ndim,                       // Number of features -- "p"
					int noptimization,              // Dimension of optimization variables -- "p"
					double *optimization_linear,    // A_O -- linear part for optimization variables
					double *optimization_state,     // O -- optimization state
					double *offset);                // h -- offset in affine transform -- "p" dimensional 

double log_density_laplace_conditional(double noise_scale,             // Scale of randomization
				       int ndim,                       // Number of features -- "p"
				       int noptimization,              // Dimension of optimization variables -- "p"
				       double *optimization_linear,    // A_O -- linear part for optimization variables
				       double *optimization_state,     // O -- optimization state
				       double *offset);                // h -- offset in affine transform -- "p" dimensional 


#ifdef __cplusplus
}  /* extern "C" */
#endif /* __cplusplus */
