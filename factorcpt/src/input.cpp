#include <Rcpp.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

NumericVector stl_sort(NumericVector x) {

   NumericVector y = clone(x);
   std::sort(y.begin(), y.end());
   return y;
}

NumericVector rcpp_rev(NumericVector x) {

   NumericVector revX = clone<NumericVector>(x);
   std::reverse(revX.begin(), revX.end());
   ::Rf_copyMostAttrib(x, revX);
   return revX;
}

// [[Rcpp::export]]
List func_dc(NumericMatrix z){

  int n = z.nrow();
  int len = z.ncol();
  int nt = 2*n;
  int i, j;
  NumericMatrix cs(n, len-1);
  NumericMatrix acs(n, len-1);
  NumericVector res(len-1);
  NumericMatrix mat(n, len-1);
  NumericVector ref;
  NumericVector col(n);
  NumericVector iplus(n);
  NumericVector iminus(n);
  double factor;

  factor = (len-1)/(double)len;
  factor = sqrt(factor);
  iplus = z(_, 0);
  for(j=0; j<n; j++){
    iminus(j) = sum(z(j, _))-z(j, 0);
    cs(j, 0) = factor*(iplus(j)-iminus(j)/(len-1));
    acs(j, 0) = fabs(cs(j, 0));
  }
  if(len-2>0){
    for(i=1; i<len-1; i++){
      factor = (i+1)*(len-i-1)/(double)len;
      factor = sqrt(factor);
      for(j=0; j<n; j++){
        iplus(j) = iplus(j)+z(j, i);
        iminus(j) = iminus(j)-z(j, i);
        cs(j, i) = factor*(iplus(j)/(i+1)-iminus(j)/(len-i-1));
        acs(j, i) = fabs(cs(j, i));
      }
    }
  }

  for(i=0; i<len-1; i++){
    ref = stl_sort(acs(_, i));
    ref = rcpp_rev(ref);

    factor = (nt-1)/(double)nt;
    factor = sqrt(factor);
    iplus(0) = ref(0);
    iminus(0) = sum(ref)-ref(0);
    col(0) = iplus(0)-iminus(0)/(nt-1);
    col(0) = factor*col(0);
    for(j=1; j<n; j++){
      factor = (j+1)*(nt-j-1)/(double)nt;
      factor = sqrt(factor);
      iplus(j) = iplus(j-1)+ref(j);
      iminus(j) = iminus(j-1)-ref(j);
      col(j) = iplus(j)/(j+1)-iminus(j)/(nt-j-1);
      col(j) = factor*col(j);
    }
    mat(_, i) = col; res(i) = max(col);
  }

  return(Rcpp::List::create(Rcpp::Named("cs")=cs, Rcpp::Named("acs")=acs, Rcpp::Named("res")=res, Rcpp::Named("mat")=mat));
}

// [[Rcpp::export]]
List func_dc_by(NumericMatrix z, double dmby, double dtby){

	int n = z.nrow();
	int nb = ceil(n/dmby);
	int len = z.ncol();
	int lenb = ceil((len-1)/dtby);
	int nt = 2*n;
	int i, j, k;
	int mby = dmby;
	int tby = dtby;

	NumericMatrix cs(n, len-1);
	NumericMatrix acs(n, len-1);
	NumericVector res(len-1);
	NumericMatrix mat(n, len-1);
	NumericVector ref;
	NumericVector col(n);
	NumericVector iplus(n);
	NumericVector iminus(n);

	double factor, plus, minus;

	factor = (len-1)/(double)len;
	factor = sqrt(factor);
	iplus = z(_, 0);
	for(j=0; j<n; j++){
		iminus(j) = sum(z(j, _))-z(j, 0);
		cs(j, 0) = factor*(iplus(j)-iminus(j)/(len-1));
		acs(j, 0) = fabs(cs(j, 0));
	}
	if(len-2>0){
		for(i=1; i<len-1; i++){
			factor = (i+1)*(len-i-1)/(double)len;
		  factor = sqrt(factor);
			for(j=0; j<n; j++){
				iplus(j) = iplus(j)+z(j, i);
				iminus(j) = iminus(j)-z(j, i);
				cs(j, i) = factor*(iplus(j)/(i+1)-iminus(j)/(len-i-1));
				acs(j, i) = fabs(cs(j, i));
			}
		}
	}

	for(i=0; i<lenb; i++){
		ref = stl_sort(acs(_, i*tby));
		ref = rcpp_rev(ref);

		factor = (nt-1)/(double)nt;
		factor = sqrt(factor);
		plus = ref(0);
		minus = sum(ref)-ref(0);
		col(0) = plus-minus/(nt-1);
		col(0) = factor*col(0);
		for(k=1; k<mby; k++){
			col(k) = col(0);
			plus += ref(k); minus -= ref(k);
		}
		for(j=1; j<nb; j++){
			factor = (j*mby+1)*(nt-j*mby-1)/(double)nt;
		  factor = sqrt(factor);
			plus += ref(j*mby);
			minus -= ref(j*mby);
			col(j*mby) = plus/(j*mby+1)-minus/(nt-j*mby-1);
			col(j*mby) = factor*col(j*mby);
			for(k=1; k<mby; k++){
				if(j*mby+k < n){
					col(j*mby+k) = col(j*mby);
					plus += ref(j*mby+k); minus -= ref(j*mby+k);
				}
			}
		}
		mat(_, i*tby) = col; res(i*tby) = max(col);

		for(k=1; k<tby; k++){
			if(i*tby+k < len-1){
		   		mat(_, i*tby+k) = col; res(i*tby+k) = res(i*tby);
			}
		}
	}

	return(Rcpp::List::create(Rcpp::Named("cs")=cs, Rcpp::Named("acs")=acs, Rcpp::Named("res")=res, Rcpp::Named("mat")=mat));
}

// [[Rcpp::export]]
NumericMatrix func_coef(NumericMatrix z, int scale){
	int n = z.nrow();
	int len = z.ncol();
	int lenw = pow(2.0, -scale);
	int t, i, j;
	NumericMatrix coef(n, len-lenw+1);
	NumericVector wave(lenw);

	for(j=0; j<lenw/2; j++){
		wave(j) = sqrt(pow(2.0, scale));
		wave(j+lenw/2) = -wave(j);
	}

	for(i=0; i<n; i++){
		for(t=0; t<len-lenw+1; t++){
			for(j=0; j<lenw; j++) coef(i, t) += z(i, t+j)*wave(j);
		}
	}

 	return(coef);
}

// [[Rcpp::export]]
NumericMatrix func_input(NumericMatrix coef, NumericMatrix sgn){
	int n = coef.nrow();
	int len = coef.ncol();
	int t, i, j, k, sg;
	int d = n*(n+1)/2;
	double avg;
	NumericMatrix input(len, d);

	k = 0;
	for(i=0; i<n; i++){
		for(t=0; t<len; t++) input(t, k) = pow(coef(i, t), 2);
		avg = mean(input(_, k));
		for(t=0; t<len; t++){
		  input(t, k) /= avg;
		  input(t, k) = std::sqrt(input(t, k));
		}

		k += 1;
		for(j=i+1; j<n; j++){
			sg = sgn(i, j);
			for(t=0; t<len; t++){
				input(t, k) = pow(coef(i, t)-sg*coef(j, t), 2);
			}
			avg = mean(input(_, k));
			for(t=0; t<len; t++){
				input(t, k) /= avg;
			  input(t, k) = std::sqrt(input(t, k));
			}
			k += 1;
		}
	}

 	return(input);
}

// [[Rcpp::export]]
NumericMatrix func_input_on(NumericMatrix coef){
	int n = coef.nrow();
	int len = coef.ncol();
	int t, i;
	NumericVector avg(n);
	NumericMatrix input(len, n);

	for(i=0; i<n; i++){
		for(t=0; t<len; t++){
			input(t, i) = pow(coef(i, t), 2);
		}
		avg(i) = mean(input(_, i));
		for(t=0; t<len; t++){
			input(t, i) /= avg(i);
		  input(t, i) = std::sqrt(input(t, i));
		}
	}
 	return(input);

}
