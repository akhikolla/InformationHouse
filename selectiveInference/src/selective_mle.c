#include <math.h> // for fabs

// Barrier for smooth MLE approximation
// after Gaussian randomization

// Described in https://arxiv.org/abs/1703.06176

//    objective = lambda u: -u.T.dot(conjugate_arg) + u.T.dot(precision).dot(u)/2. + np.log(1.+ 1./(u / scaling)).sum()

double barrier_objective(double *opt_variable,               // Optimization variable
			 double *conjugate_arg,              // Argument to conjugate of Gaussian
			 double *precision,                  // Precision matrix of Gaussian
			 double *scaling,                    // Diagonal scaling matrix for log barrier
			 int ndim)                           // Dimension of conjugate_arg, precision
{
  int idim, jdim;
  double *conjugate_arg_ptr;
  double *opt_variable_ptr;
  double *precision_ptr;
  double *scaling_ptr;
  double product_entry, value;

  value = 0.;
  for (idim=0; idim<ndim; idim++) {
    
    // Compute an entry of precision.dot(conjugate_arg)

    product_entry = 0;
    for (jdim=0; jdim<ndim; jdim++) {
    
      precision_ptr = ((double *) precision + idim * ndim + jdim); // precision is a symmetric matrix
      opt_variable_ptr = ((double *) opt_variable + jdim);
      product_entry += (*precision_ptr) * (*opt_variable_ptr);
    }
      
    opt_variable_ptr = ((double *) opt_variable + idim);
    value += 0.5 * (*opt_variable_ptr) * product_entry;

    // now add linear term

    conjugate_arg_ptr = ((double *) conjugate_arg + idim);
    value -= (*conjugate_arg_ptr) * (*opt_variable_ptr);

    // now log term

    scaling_ptr = ((double *) scaling + idim);
    value += log(((*opt_variable_ptr) + (*scaling_ptr)) / (*opt_variable_ptr));

  }

  return(value);
}

//    grad = lambda u: -conjugate_arg + precision.dot(u) + (1./(scaling + u) - 1./u)

void barrier_gradient(double *gradient,                   // Gradient vector
		      double *opt_variable,               // Optimization variable
		      double *conjugate_arg,              // Argument to conjugate of Gaussian
		      double *precision,                  // Precision matrix of Gaussian
		      double *scaling,                    // Diagonal scaling matrix for log barrier
		      int ndim)                           // Dimension of conjugate_arg, precision
{
  int idim, jdim;
  double *gradient_ptr;
  double *conjugate_arg_ptr;
  double *opt_variable_ptr;
  double *precision_ptr;
  double *scaling_ptr;
  double product_entry;

  for (idim=0; idim<ndim; idim++) {
    
    gradient_ptr = ((double *) gradient + idim);

    // Compute an entry of precision.dot(conjugate_arg)

    product_entry = 0;
    for (jdim=0; jdim<ndim; jdim++) {
    
      precision_ptr = ((double *) precision + idim * ndim + jdim); // precision is a symmetric matrix
      opt_variable_ptr = ((double *) opt_variable + jdim);
      product_entry += (*precision_ptr) * (*opt_variable_ptr);
    }
      
    opt_variable_ptr = ((double *) opt_variable + idim);
    *gradient_ptr = product_entry;

    // now add linear term

    conjugate_arg_ptr = ((double *) conjugate_arg + idim);
    *gradient_ptr -= (*conjugate_arg_ptr);

    // now log term

    scaling_ptr = ((double *) scaling + idim);
    *gradient_ptr += 1. / ((*opt_variable_ptr) + (*scaling_ptr)) - 1. / (*opt_variable_ptr);
  }

}

double barrier_gradient_step(double *gradient,                   // Gradient vector
			     double *opt_variable,               // Optimization variable
			     double *opt_proposed,               // Proposed value of optimization variable
			     double *conjugate_arg,              // Argument to conjugate of Gaussian
			     double *precision,                  // Precision matrix of Gaussian
			     double *scaling,                    // Diagonal scaling matrix for log barrier
			     double step,                        // Step size for gradient step
			     int ndim)                           // Dimension of conjugate_arg, precision
{
  int idim;
  double *gradient_ptr;
  double *opt_variable_ptr;
  double *opt_proposed_ptr;

  for (idim=0; idim<ndim; idim++) {
    opt_variable_ptr = ((double *) opt_variable + idim);
    opt_proposed_ptr = ((double *) opt_proposed + idim);
    gradient_ptr = ((double *) gradient + idim);
    *opt_proposed_ptr = (*opt_variable_ptr) - step * (*gradient_ptr);
  }

  return barrier_objective(opt_proposed,
			   conjugate_arg,         
			   precision,             
			   scaling,               
			   ndim);
}

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
		     double initial_step)                // Initial step size
{
  int iter, idim, istep;
  double *gradient_ptr;
  double *opt_variable_ptr;
  double *opt_proposed_ptr;
  double current_value = barrier_objective(opt_variable,
					   conjugate_arg,         
					   precision,             
					   scaling,               
					   ndim);
  double proposed_value;
  double step = initial_step;
  int any_negative;

  for (iter=0; iter<max_iter; iter++) {
    
    // Compute the gradient

    barrier_gradient(gradient,    
		     opt_variable,
		     conjugate_arg,
		     precision,    
		     scaling,
		     ndim);
    
    // Find a valid step size

    istep = 0;
    while (1) {
      any_negative = 0;
      for (idim=0; idim<ndim; idim++) {
	opt_variable_ptr = ((double *) opt_variable + idim);
	gradient_ptr = ((double *) gradient + idim);
	if ((*opt_variable) - step * (*gradient_ptr) < 0) {
	  any_negative += 1;
	}
      }
      if (any_negative == 0) {
	break;
      }
      step = step * 0.5;
      istep++;
      if (istep == 50) {
	break;            // Terminate eventually -- but this will be a failure.
	                  // Should mean that opt_variable is in feasible.
      }
    }

    // Find a step size that is a descent

    istep = 0;
    while (1) {
      proposed_value = barrier_gradient_step(gradient,
					     opt_variable,
					     opt_proposed,
					     conjugate_arg,         
					     precision,             
					     scaling,               
					     step,
					     ndim);
      if (proposed_value < current_value) {
	for (idim=0; idim<ndim; idim++) {
	  opt_variable_ptr = ((double *) opt_variable + idim);
	  opt_proposed_ptr = ((double *) opt_proposed + idim);
	  *opt_variable_ptr = *opt_proposed_ptr;
	  current_value = proposed_value;
	}
	break;
      }
      step = step * 0.5;
      istep++;
      if (istep == 50) {
	break;            // Terminate eventually -- this will mean no descent.
	                  // We've solved the problem
      }
    }

    if ((fabs(current_value - proposed_value) < value_tol * fmax(fabs(current_value), 1)) & 
	(iter >= min_iter)) {
      current_value = proposed_value;
      break;
    }
    current_value = proposed_value;
    
  }

  return(current_value);

}

// Below, the constraints are composed with an affine function, i.e. Au \leq b, rather than u \geq 0

void set_affine_term(double *opt_variable,               // Optimization variable
		     double *linear_term,                // Matrix A in constraint Au \leq b
		     double *offset,                     // Offset b in constraint Au \leq b
		     double *affine_term,                // Should be equal to b - A.dot(opt_variable)    
		     int ndim,                           // Dimension of conjugate_arg, precision
		     int ncon)                           // Number of constraints
{
  int idim, icon;
  double *linear_term_ptr;
  double *affine_term_ptr;
  double *offset_ptr;
  double *opt_variable_ptr;
  double affine_entry;

  for (icon=0; icon<ncon; icon++) {

    affine_entry = 0;

    for (idim=0; idim<ndim; idim++) {
      linear_term_ptr = ((double *) linear_term + icon + idim * ncon);  // matrices are in column-major order for R!
      opt_variable_ptr = ((double *) opt_variable + idim);
      affine_entry -= (*linear_term_ptr) * (*opt_variable_ptr);
    }

    offset_ptr = ((double *) offset + icon);
    affine_entry += (*offset_ptr);      // one entry of b-Au

    affine_term_ptr = ((double *) affine_term + icon);
    *affine_term_ptr = affine_entry;
  }

}

//    objective = lambda u: -u.T.dot(conjugate_arg) + u.T.dot(precision).dot(u)/2. + np.log(1.+ 1./(b - A.dot(u) / scaling)).sum()

double barrier_objective_affine(double *opt_variable,               // Optimization variable
				double *conjugate_arg,              // Argument to conjugate of Gaussian
				double *precision,                  // Precision matrix of Gaussian
				double *scaling,                    // Diagonal scaling matrix for log barrier
				double *linear_term,                // Matrix A in constraint Au \leq b
				double *offset,                     // Offset b in constraint Au \leq b
				double *affine_term,                // Should be equal to b - A.dot(opt_variable)    
				int ndim,                           // Dimension of conjugate_arg, precision
				int ncon)                           // Number of constraints
{
  int idim, jdim, icon;
  double *conjugate_arg_ptr;
  double *opt_variable_ptr;
  double *precision_ptr;
  double *scaling_ptr;
  double *affine_term_ptr;
  double product_entry, value;

  set_affine_term(opt_variable,
		  linear_term,             
		  offset,                  
		  affine_term,             
		  ndim,                        
		  ncon);

  value = 0.;
  for (idim=0; idim<ndim; idim++) {
    
    // Compute an entry of precision.dot(conjugate_arg)

    product_entry = 0;
    for (jdim=0; jdim<ndim; jdim++) {
    
      precision_ptr = ((double *) precision + idim * ndim + jdim); // precision is a symmetric matrix
      opt_variable_ptr = ((double *) opt_variable + jdim);
      product_entry += (*precision_ptr) * (*opt_variable_ptr);
    }
      
    opt_variable_ptr = ((double *) opt_variable + idim);
    value += 0.5 * (*opt_variable_ptr) * product_entry;

    // now add linear term

    conjugate_arg_ptr = ((double *) conjugate_arg + idim);
    value -= (*conjugate_arg_ptr) * (*opt_variable_ptr);

  }

  // now log term

  for (icon=0; icon<ncon; icon++) {
    scaling_ptr = ((double *) scaling + icon);
    affine_term_ptr = ((double *) affine_term + icon);
    value += log(((*affine_term_ptr) + (*scaling_ptr)) / (*affine_term_ptr));
  }
  
  return(value);
}

//    grad = lambda u: -conjugate_arg + precision.dot(u) + (1./(scaling + u) - 1./u)

void barrier_gradient_affine(double *gradient,                   // Gradient vector
			     double *opt_variable,               // Optimization variable
			     double *conjugate_arg,              // Argument to conjugate of Gaussian
			     double *precision,                  // Precision matrix of Gaussian
			     double *scaling,                    // Diagonal scaling matrix for log barrier
			     double *linear_term,                // Matrix A in constraint Au \leq b
			     double *offset,                     // Offset b in constraint Au \leq b
			     double *affine_term,                // Should be equal to b - A.dot(opt_variable)    
			     int ndim,                           // Dimension of conjugate_arg, precision
			     int ncon)                           // Number of constraints
{
  int idim, jdim, icon;
  double *gradient_ptr;
  double *conjugate_arg_ptr;
  double *opt_variable_ptr;
  double *precision_ptr;
  double *scaling_ptr;
  double *affine_term_ptr;
  double *linear_term_ptr;
  double product_entry;

  set_affine_term(opt_variable,
		  linear_term,             
		  offset,                  
		  affine_term,             
		  ndim,                        
		  ncon);

  for (idim=0; idim<ndim; idim++) {
    
    gradient_ptr = ((double *) gradient + idim);

    // Compute an entry of precision.dot(conjugate_arg)

    product_entry = 0;
    for (jdim=0; jdim<ndim; jdim++) {
    
      precision_ptr = ((double *) precision + idim * ndim + jdim); // precision is a symmetric matrix
      opt_variable_ptr = ((double *) opt_variable + jdim);
      product_entry += (*precision_ptr) * (*opt_variable_ptr);
    }
      
    *gradient_ptr = product_entry;

    // now add linear term

    conjugate_arg_ptr = ((double *) conjugate_arg + idim);
    *gradient_ptr -= (*conjugate_arg_ptr);

    // now log term: A.T.dot(gradient of barrier)

    for (icon=0; icon<ncon; icon++) {
      scaling_ptr = ((double *) scaling + icon);
      affine_term_ptr = ((double *) affine_term + icon);
      linear_term_ptr = ((double *) linear_term + ncon * idim + icon);   // matrices are in column-major order for R!
      *gradient_ptr -= (*linear_term_ptr) * (1. / ((*affine_term_ptr) + (*scaling_ptr)) - 1. / (*affine_term_ptr));
    }
  }
}

double barrier_gradient_step_affine(double *gradient,                   // Gradient vector
				    double *opt_variable,               // Optimization variable
				    double *opt_proposed,               // Proposed value of optimization variable
				    double *conjugate_arg,              // Argument to conjugate of Gaussian
				    double *precision,                  // Precision matrix of Gaussian
				    double *scaling,                    // Diagonal scaling matrix for log barrier
				    double *linear_term,                // Matrix A in constraint Au \leq b
				    double *offset,                     // Offset b in constraint Au \leq b
				    double *affine_term,                // Should be equal to b - A.dot(opt_variable)   
				    double step,                        // Step size for gradient step
				    int ndim,                           // Dimension of conjugate_arg, precision
				    int ncon)                           // Number of constraints
{
  int idim;
  double value;
  double *gradient_ptr;
  double *opt_variable_ptr;
  double *opt_proposed_ptr;

  // Gradient step stored in opt_proposed

  for (idim=0; idim<ndim; idim++) {
    opt_variable_ptr = ((double *) opt_variable + idim);
    opt_proposed_ptr = ((double *) opt_proposed + idim);
    gradient_ptr = ((double *) gradient + idim);
    *opt_proposed_ptr = (*opt_variable_ptr) - step * (*gradient_ptr);
  }

  value = barrier_objective_affine(opt_proposed,
				   conjugate_arg,         
				   precision,             
				   scaling,               
				   linear_term,
				   offset,
				   affine_term,
				   ndim,
				   ncon);

  return(value);
}

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
			    double initial_step)                // Initial step size
{
  int iter, idim, istep, icon;
  double *gradient_ptr;
  double *opt_variable_ptr;
  double *opt_proposed_ptr;
  double current_value = barrier_objective_affine(opt_variable,
						  conjugate_arg,         
						  precision,             
						  scaling,               
						  linear_term,
						  offset,
						  affine_term,
						  ndim,
						  ncon);
  double proposed_value;
  double step = initial_step;
  double affine_updated;
  double *linear_term_ptr;
  double *affine_term_ptr;
  int any_negative;

  for (iter=0; iter<max_iter; iter++) {
    
    // Compute the gradient

    barrier_gradient_affine(gradient,    
			    opt_variable,
			    conjugate_arg,
			    precision,    
			    scaling,
			    linear_term,
			    offset,
			    affine_term,
			    ndim,
			    ncon);
    
    // Find a valid step size

    istep = 0;
    while (1) {
      any_negative = 0;

      for (icon=0; icon<ncon; icon++) {
	affine_term_ptr = ((double *) affine_term + icon);  // an entry of b-A.dot(opt)
	affine_updated = (*affine_term_ptr);

	for (idim=0; idim<ndim; idim++) {
	  linear_term_ptr = ((double *) linear_term + icon + idim * ncon);
	  gradient_ptr = ((double *) gradient + idim);
	  affine_updated += (*linear_term_ptr) * (*gradient_ptr) * step;
	}
	if (affine_updated < 0) {
	  any_negative += 1;
	}
      }

      if (any_negative == 0) {
	break;
      }
      step = step * 0.5;
      istep++;

      if (istep == 50) {
	break;            // Terminate eventually -- but this will be a failure.
	                  // Should mean that opt_variable is in feasible.
      }
    }

    // Find a step size that is a descent

    istep = 0;
    while (1) {
      proposed_value = barrier_gradient_step_affine(gradient,
						    opt_variable,
						    opt_proposed,
						    conjugate_arg,         
						    precision,             
						    scaling,               
						    linear_term,
						    offset,
						    affine_term,
						    step,
						    ndim,
						    ncon);

      if (proposed_value < current_value) {
	for (idim=0; idim<ndim; idim++) {
	  opt_variable_ptr = ((double *) opt_variable + idim);
	  opt_proposed_ptr = ((double *) opt_proposed + idim);
	  *opt_variable_ptr = *opt_proposed_ptr;
	}
	current_value = proposed_value;
	break;
      }
      step = step * 0.5;
      istep++;
      if (istep == 50) {
	break;            // Terminate eventually -- this will mean no descent.
	                  // We've solved the problem
      }
    }

    if ((fabs(current_value - proposed_value) < value_tol * fmax(fabs(current_value), 1)) & 
	(iter >= min_iter)) {
      current_value = proposed_value;
      break;
    }
    
  }

  return(current_value);

}

