#ifndef __twosided_noncrossing_probability_h__
#define __twosided_noncrossing_probability_h__

#include <vector>

// Compute the probability that a homogeneous Poisson process in [0,1] of a given intensity will stay within
// lower and upper boundary functions.

std::vector<double> poisson_process_noncrossing_probability(double intensity, const std::vector<double>& lower_bound_steps, const std::vector<double>& upper_bound_steps, bool use_fft);

// Compute the probability that an empirical process eta_{n}(t) with n points in [0,1] will stay within lower and upper boundary functions:


double ecdf_noncrossing_probability(int n, const std::vector<double>& lower_bound_steps, const std::vector<double>& upper_bound_steps, bool use_fft);



#endif
