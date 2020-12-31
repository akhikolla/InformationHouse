#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

//[[Rcpp::export]]
NumericMatrix logliknormLR(NumericVector yy, NumericMatrix XX_aug, double alpha, int init){

	// Inputs: 
		// y: time series vector [y1 ... yT]
		// X: matrix of covariates, not including the intercept column which will be added.
		//	   	If X=FALSE, no covariates - intercept model only.
	// alpha: PWD parameter we're calculating the marginal predictive loglikelihood for.
	// y[1] should be y_1 (the time point furthest into the past)
	// init: time point to begin computing marginal predictive likelihoods.

	// Cast NumericMatrix etc. to arma types
	vec y = as<arma::vec>(yy);
	mat X_aug = as<arma::mat>(XX_aug);

	// Initialize number of time points and predictors
	int n = X_aug.n_rows;
	int p = X_aug.n_cols;

	// Initialize what will get recursively updated
	double sum_a = 1; 
	double yay = y(0)*y(0);
	arma::mat xay = y(0)*X_aug.row(0).t();
	arma::mat xax = X_aug.row(0).t() * X_aug.row(0);
	arma::mat output(1,1);
	output.zeros();

	// Other initializations to be used within for loop below
	double yj = 0;
	arma::mat xj(p,1);
	xj.zeros();
	arma::mat xax_inv(p,p);
	xax_inv.zeros();
	arma::mat bhat(p,1);
	bhat.zeros();
	double nu = 0;
	arma::mat S_j(1,1);
	S_j.zeros();
	arma::mat S_hat(1,1);
	S_hat.zeros();
	arma::mat inner_term(1,1);
	inner_term.zeros();
	arma::mat mu_hat(1,1);
	mu_hat.zeros();
	arma::mat yaxb(1,1);
	yaxb.zeros();
	arma::mat x_star(p,1);
	bhat.zeros();

	for (int j = 1; j < init-1; j++){
		sum_a = sum_a * alpha + 1;
		yj = y(j);
		xj = X_aug.row(j).t();
		yay = yay * alpha + yj*yj;
		xay = xay * alpha + yj*xj;
		xax = xax * alpha + xj*xj.t();
		}
	for (int j = init-1; j < n-1; j++ ){
		sum_a = sum_a * alpha + 1; 
		yj = y(j); 
		xj = X_aug.row(j).t();
		yay = yay * alpha + yj*yj;
		xay = xay * alpha + yj*xj;
		xax = xax * alpha + xj*xj.t();
		xax_inv = inv(xax);
		bhat = xax_inv * xay;
		nu = sum_a - p;
		yaxb = xay.t() * bhat;
		S_hat = yay - yaxb;
		x_star = X_aug.row(j+1).t();
		mu_hat = x_star.t()*bhat;
		S_j = S_hat * (1 + x_star.t() * xax_inv * x_star);
		inner_term = (y(j+1) - mu_hat)*(y(j+1) - mu_hat)/S_j;
		output += lgamma((nu+1)/2) - lgamma(nu/2) - log(S_j)/2 - (nu+1)*log(1+inner_term)/2;
		}	
	return as<NumericMatrix>(wrap(output));
	}
