#ifdef __cplusplus
extern "C"
{
#endif /* __cplusplus */

double barrier_solve(double *gradient,                   // Gradient vector
		     double *opt_variable,               // Optimization variable
		     double *opt_proposed,               // New value of optimization variable
		     double *conjugate_arg,              // Argument to conjugate of Gaussian
		     double *precision,                  // Precision matrix of Gaussian
		     double *scaling,                    // Diagonal scaling matrix for log barrier
		     int ndim,                           // Dimension of conjugate_arg, precision
		     int max_iter,                       // Maximum number of iterations
		     int min_iter,                       // Minimum number of iterations
		     double value_tol,                   // Tolerance for convergence based on value
		     double initial_step);               // Initial stepsize 

double barrier_solve_affine(double *gradient,                   // Gradient vector
			    double *opt_variable,               // Optimization variable
			    double *opt_proposed,               // New value of optimization variable
			    double *conjugate_arg,              // Argument to conjugate of Gaussian
			    double *precision,                  // Precision matrix of Gaussian
			    double *scaling,                    // Diagonal scaling matrix for log barrier
			    double *linear_term,                // Matrix A in constraint Au \leq b
			    double *offset,                     // Offset b in constraint Au \leq b
			    double *affine_term,                // Should be equal to b - A.dot(opt_variable)    
			    int ndim,                           // Dimension of conjugate_arg, precision
			    int ncon,                           // Number of constraints
			    int max_iter,                       // Maximum number of iterations
			    int min_iter,                       // Minimum number of iterations
			    double value_tol,                   // Tolerance for convergence based on value
			    double initial_step);               // Initial step size

#ifdef __cplusplus
}  /* extern "C" */
#endif /* __cplusplus */
