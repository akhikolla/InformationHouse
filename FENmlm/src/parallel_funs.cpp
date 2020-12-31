#include <cmath>
#include <Rcpp.h>
#include <stdio.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <Rmath.h>
using namespace Rcpp;

// [[Rcpp::plugins(openmp)]]

// This file contains femlm functions that will be parallelized with the omp library

double cpppar_abs(double x){
	//simple function to compute the absolute value
	if(x > 0){
		return(x);
	} else {
		return(-x);
	}
}

// [[Rcpp::export]]
NumericVector cpppar_exp(NumericVector x, int nthreads){
	// parallel exponentiation using omp

	int n = x.length();
	NumericVector res(n);

	#pragma omp parallel num_threads(nthreads)
	{
		#pragma omp for
		for(int i = 0 ; i < n ; i++) {
			res[i] = exp(x[i]);
		}
	}

	return(res);
}

// [[Rcpp::export]]
NumericVector cpppar_log(NumericVector x, int nthreads){
	// parallel exponentiation using omp

	int n = x.length();
	NumericVector res(n);

	#pragma omp parallel num_threads(nthreads)
	{
		#pragma omp for
		for(int i = 0 ; i < n ; i++) {
			res[i] = log(x[i]);
		}
	}

	return(res);
}

// [[Rcpp::export]]
NumericVector cpppar_log_a_exp(int nthreads, double a, NumericVector mu, NumericVector exp_mu){
	// faster this way

	int n = mu.length();
	NumericVector res(n);

	#pragma omp parallel num_threads(nthreads)
	{
		#pragma omp for
		for(int i=0 ; i<n ; i++) {
			if(mu[i] < 200){
				res[i] = log(a + exp_mu[i]);
			} else {
				res[i] = mu[i];
			}
		}
	}

	return(res);
}

// [[Rcpp::export]]
NumericVector cpppar_lgamma(NumericVector x, int nthreads){
	// parallel lgamma using omp

	int n = x.length();
	NumericVector res(n);

	#pragma omp parallel num_threads(nthreads)
	{
		#pragma omp for
		for(int i = 0 ; i < n ; i++) {
			res[i] = lgamma(x[i]);
		}
	}

	return(res);
}

// [[Rcpp::export]]
NumericVector cpppar_digamma(NumericVector x, int nthreads){
	// parallel digamma using omp

	int n = x.length();
	NumericVector res(n);

	#pragma omp parallel num_threads(nthreads)
	{
		#pragma omp for
		for(int i = 0 ; i < n ; i++) {
			res[i] = R::digamma(x[i]);
		}
	}

	return(res);
}

// [[Rcpp::export]]
NumericVector cpppar_trigamma(NumericVector x, int nthreads){
	// parallel trigamma using omp

	int n = x.length();
	NumericVector res(n);

	#pragma omp parallel num_threads(nthreads)
	{
		#pragma omp for
		for(int i = 0 ; i < n ; i++) {
			res[i] = R::trigamma(x[i]);
		}
	}

	return(res);
}







































































