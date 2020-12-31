#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

/* c_bin_alloc identifies the bin membership for a vector x
* based on equal width bins.
* Arguments:
* x = vector of points to allocate to bins
* nbin = number of bins
* min, max = lower and upper bound for the first and last bin
*/

// [[Rcpp::export]]


NumericVector cbin_alloc(NumericVector x, int nbin, double min, double max){
  /* setup parameters and object to be returned*/
  int n = x.size();
  NumericVector output(n);
  double skip = (nbin - 1.0) / (max - min); /* skip is the width of a bin */

  /* i-th output is based on the floor of the multiple of bins from the lower bound, min */
  for(int i = 0; i < n; i++){
    output[i] = floor((x[i] - min) * skip) + 1;
    if(output[i] < 0) output[i] = 0;
    if(output[i] > n) output[i] = n;
  }
  return output;
}





/* sm_bin_wts computes the (smoothed) aggregate values from vector y falling in
* each bin, based on bin membership of associated vector x.
* Arguments:
* x = vector of points to determine allocation to bins
* y = vector of values to distribute among bins
* nbin = number of bins
* min, max = lower and upper bound for the first and last bin
*/


// [[Rcpp::export]]


NumericVector sm_bin_wts(NumericVector x, NumericVector y, int nbin, double min, double max){
  /* setup parameters and object to be returned */
  int n = x.size();
  NumericVector output(nbin);
  double skip = (nbin - 1.0) / (max - min);  /* skip is the width of a bin */

  int fl = 0; /* fl used to capture the bin to the left of a point */
  double xe; /* (xe - fl) determines weights for the adjacent bins */

  /* loop over points in x, and split associated values in y among adjacent bins */
  for(int i = 0; i < n; i++){
    xe = (x[i] - min) * skip;
    fl = floor(xe);

    if(fl < (nbin - 1) && fl >= 0){
      output[fl+1] += (xe - fl) * y[i];
      output[fl] += (fl + 1 - xe) * y[i];
    }
    else if(fl >= (nbin - 1)) output[nbin-1] += y[i];
    else output[0] += y[i];
  }
  return output;
}


/* bin_wts is as above, but with hard allocation of points to bins
* Interpretation is as above, but no sharing of points
* between bins is performed
*/


// [[Rcpp::export]]


NumericVector bin_wts(NumericVector x, NumericVector y, int nbin, double min, double max){
  int n = x.size();
  NumericVector output(nbin);
  double skip = (nbin - 1.0) / (max - min);
  int fl = 0;
  double xe;
  for(int i = 0; i < n; i++){
    xe = (x[i] - min) * skip;
    fl = floor(xe);
    if(fl < (nbin - 1) && fl >= 0){
      output[fl] += y[i];
    }
    else if(fl >= (nbin - 1)) output[nbin-1] += y[i];
    else output[0] += y[i];
  }
  return output;
}



/* ksum computes kernel sums at locations x_eval based on points x with coefficients/weights y.
* pair (x, y) should be sorted according to values in x (non-decreasing).
* Arguments:
* x = kernel locations. Usually the observed sample, or bin locations for binned summing
* y = corresponding coefficients to be added in kernel weighted sums
* x_eval = points at which to evaluate the sums
* h = bandwidth value. Must be positive
* betas = vector of kernel coefficients
* Counts = optional, identifying the location of x_eval w.r.t. x. That is, Counts[i] is the number of values in x
* less than or equal to x_eval[i]. If omitted then these will be computed. If omitted, then x_eval must also be
* sorted in non-decreasing order.
*/



// [[Rcpp::export]]

NumericVector ksum(NumericVector x, NumericVector y, NumericVector x_eval, double h, NumericVector betas, NumericVector Counts = NumericVector(1)){
  /* setup parameters and object to be returned */
  int n = x.size();
  int n_eval = x_eval.size();
  int ord = betas.size() - 1;
  NumericVector output(n_eval);
  double denom;
  double exp_mult;

  /* Ly and Ry store the recursive sums used in fast kernel computations. See Fast Exact Evaluation of Univariate Kernel Sums (Hofmeyr, 2019)
  * for details
  */
  NumericMatrix Ly(ord + 1, n);
  NumericMatrix Ry(ord + 1, n);
  for(int i = 0; i <= ord; i++) Ly(i,0) = pow(-x[0], i) * y[0];
  for(int i = 1; i < n; i++){
    for(int j = 0; j <= ord; j++){
      Ly(j, i) = pow(-x[i], j) * y[i] + exp((x[i-1] - x[i]) / h) * Ly(j, i - 1);
      Ry(j, n - i - 1) = exp((x[n - i - 1] - x[n - i]) / h)*(pow(x[n - i], j) * y[n-i] + Ry(j, n - i));
    }
  }

  /* the function operates differently depending on whether or not Counts have been supplied. */
  if(Counts.size() == 1){
    /* Counts not supplied so need to be computed, now as counts */
    NumericVector counts(n_eval);
    int count = 0;
    /* counts determined by looping, and since x and x_eval are sorted
    * counts[i+1] >= counts[i], meaning counting does not need to
    * restart for each
    */
    for(int i = 0; i < n_eval; i++){
      if(x_eval[i] >= x[n - 1]){
        counts[i] = n;
      }
      else{
        while(count < n && x[count] <= x_eval[i]){
          count += 1;
        }
        counts[i] = count;
      }
    }

    /* next loop over the terms in the polynomial in the kernel expression.
    * orddo represents the exponent, increasing to the order of the polynomial
    */
    for(int orddo = 0; orddo <= ord; orddo++){
      NumericVector coefs(orddo + 1); /* coefs are the combinatorial coefficients in the binomial expansion */
      coefs[0] = coefs[orddo] = 1;
      if(orddo>1){
        double num = 1;
        for(int j = 2; j <= orddo; j++) num *= j;
        double denom1 = 1;
        double denom2 = num / orddo;
        for(int i = 2; i <= orddo; i++){
          coefs[i - 1] = num / denom1 / denom2;
          denom1 *= i;
          denom2 /= (orddo - i + 1);
        }
      }
      denom = pow(h, orddo); /* denominator for corresponding exponent is h^orddo */
      int ix; /* ix is the location of x_eval values in the x values, i.e., term in counts */

      /* next loop over evaluation points and compute the contribution to the kernel sum
       * from current polynomial term in the kernel
       */
      for(int i = 0; i < n_eval; i++){
        ix = round(counts[i]);
        if(ix == 0){
          exp_mult = exp((x_eval[i] - x[0]) / h);
          output[i] += betas[orddo] * pow(x[0] - x_eval[i], orddo) / denom * exp_mult * y[0];
          for(int j = 0; j <= orddo; j++) output[i] += betas[orddo] * coefs[j] * pow(-x_eval[i], orddo - j) * Ry(j, 0) / denom * exp_mult;
        }
        else{
          exp_mult = exp((x[ix-1] - x_eval[i]) / h);
          for(int j = 0; j <= orddo; j++) output[i] += betas[orddo] * coefs[j] * (pow(x_eval[i], orddo - j) * Ly(j, ix - 1) * exp_mult + pow(-x_eval[i], orddo - j) * Ry(j, ix - 1) / exp_mult) / denom;
        }
      }
    }
  }
  else{
    /* if Counts are supplied, then do not need to be determined.
     * Ecerything in this block is exactly as above, but with Counts
     * instead of counts (and without determination of counts)
     */
    for(int orddo = 0; orddo <= ord; orddo++){
      NumericVector coefs(orddo + 1);
      coefs[0] = coefs[orddo] = 1;
      if(orddo > 1){
        double num = 1;
        for(int j = 2; j <= orddo; j++) num *= j;
        double denom1 = 1;
        double denom2 = num / orddo;
        for(int i = 2; i <= orddo; i++){
          coefs[i - 1] = num / denom1 / denom2;
          denom1 *= i;
          denom2 /= (orddo - i + 1);
        }
      }
      denom = pow(h, orddo);
      int ix;
      for(int i = 0; i < n_eval; i++){
        ix = round(Counts[i]);
        if(ix == 0){
          exp_mult = exp((x_eval[i] - x[0]) / h);
          output[i] += betas[orddo] * pow(x[0] - x_eval[i], orddo) / denom * exp_mult * y[0];
          for(int j = 0; j <= orddo; j++) output[i] += betas[orddo] * coefs[j] * pow(-x_eval[i], orddo - j) * Ry(j, 0) / denom * exp_mult;
        }
        else{
          exp_mult = exp((x[ix-1] - x_eval[i]) / h);
          for(int j = 0; j <= orddo; j++) output[i] += betas[orddo] * coefs[j] * (pow(x_eval[i], orddo - j) * Ly(j, ix - 1) * exp_mult + pow(-x_eval[i], orddo - j) * Ry(j, ix - 1) / exp_mult) / denom;
        }
      }
    }
  }
  return output;
}



/* dksum computes kernel derivative sums at locations x_eval based on points x with coefficients/weights y.
*
* The interpretation of the entire function and its arguments is exactly the same as in ksum, except where indicated
*/

// [[Rcpp::export]]

NumericVector dksum(NumericVector x, NumericVector y, NumericVector x_eval, double h, NumericVector betas, NumericVector Counts = NumericVector(1)){
  int n = x.size();
  int n_eval = x_eval.size();
  int ord = betas.size() - 1;
  NumericVector output(n_eval);
  double denom;
  double exp_mult;

  /* vector tbetas is the corresponding set of coefficients for the derivative of the kernel */
  NumericVector tbetas(ord + 1);
  for(int k = 0; k < ord; k++) tbetas[k] = (k + 1) * betas[k + 1] - betas[k];
  tbetas[ord] = -betas[ord];
  NumericMatrix Ly(ord + 1, n);
  NumericMatrix Ry(ord + 1, n);
  for(int i = 0; i <= ord; i++) Ly(i, 0) = pow(-x[0], i) * y[0];
  for(int i = 1; i < n; i++){
    for(int j = 0; j <= ord; j++){
      Ly(j, i) = pow(-x[i], j) * y[i] + exp((x[i - 1] - x[i]) / h) * Ly(j, i - 1);
      Ry(j, n - i - 1) = exp((x[n - i - 1] - x[n - i]) / h) * (pow(x[n - i], j) * y[n - i] + Ry(j, n - i));
    }
  }
  if(Counts.size() == 1){
    int count = 0;
    NumericVector counts(n_eval);
    for(int i = 0; i < n_eval; i++){
      if(x_eval[i] >= x[n - 1]){
        counts[i] = n;
      }
      else{
        while(count < n && x[count] <= x_eval[i]){
          count += 1;
        }
        counts[i] = count;
      }
    }
    for(int orddo = 0; orddo <= ord; orddo++){
      NumericVector coefs(orddo + 1);
      coefs[0] = coefs[orddo] = 1;
      if(orddo > 1){
        double num = 1;
        for(int j = 2; j <= orddo; j++) num *= j;
        double denom1 = 1;
        double denom2 = num / orddo;
        for(int i = 2; i <= orddo; i++){
          coefs[i - 1] = num / denom1 / denom2;
          denom1 *= i;
          denom2 /= (orddo - i + 1);
        }
      }
      denom = pow(h, orddo);
      int ix;
      for(int i = 0; i < n_eval; i++){
        ix = round(counts[i]);
        if(ix==0){
          exp_mult = exp((x_eval[i] - x[0]) / h);
          output[i] += tbetas[orddo] * pow(x[0] - x_eval[i], orddo) / denom * exp_mult * y[0];
          for(int j = 0; j <= orddo; j++) output[i] += tbetas[orddo] * coefs[j] * pow(-x_eval[i], orddo - j) * Ry(j, 0) / denom * exp_mult;
        }
        else{
          exp_mult = exp((x[ix - 1] - x_eval[i]) / h);
          for(int j = 0; j <= orddo; j++) output[i] -= tbetas[orddo] * coefs[j] * (pow(x_eval[i], orddo - j) * Ly(j, ix - 1) * exp_mult - pow(-x_eval[i], orddo - j) * Ry(j, ix - 1) / exp_mult) / denom;
        }
      }
    }
  }
  else{
    for(int orddo = 0; orddo <= ord; orddo++){
      NumericVector coefs(orddo + 1);
      coefs[0] = coefs[orddo] = 1;
      if(orddo > 1){
        double num = 1;
        for(int j = 2; j <= orddo; j++) num *= j;
        double denom1 = 1;
        double denom2 = num / orddo;
        for(int i = 2; i <= orddo; i++){
          coefs[i - 1] = num / denom1 / denom2;
          denom1 *= i;
          denom2 /= (orddo - i + 1);
        }
      }
      denom = pow(h, orddo);
      int ix;
      for(int i = 0; i < n_eval; i++){
        ix = round(Counts[i]);
        if(ix == 0){
          exp_mult = exp((x_eval[i] - x[0]) / h);
          output[i] += tbetas[orddo] * pow(x[0] - x_eval[i], orddo) / denom * exp_mult * y[0];
          for(int j = 0; j <= orddo; j++) output[i] += tbetas[orddo] * coefs[j] * pow(-x_eval[i], orddo - j) * Ry(j, 0) / denom * exp_mult;
        }
        else{
          exp_mult = exp((x[ix - 1] - x_eval[i]) / h);
          for(int j = 0; j <= orddo; j++) output[i] -= tbetas[orddo] * coefs[j] * (pow(x_eval[i], orddo - j) * Ly(j, ix - 1) * exp_mult - pow(-x_eval[i], orddo - j) * Ry(j, ix - 1) / exp_mult) / denom;
        }
      }
    }
  }
  return output;
}



/* kndksum computes both kernel and kernel derivative sums at locations x_eval based on points x with coefficients/weights y.
*
* The interpretation of the entire function and its arguments is exactly the same as in ksum and dksum.
*/

// [[Rcpp::export]]

NumericMatrix kndksum(NumericVector x, NumericVector y, NumericVector x_eval, double h, NumericVector betas, NumericVector Counts = NumericVector(0)){
  int n = x.size();
  int n_eval = x_eval.size();
  int ord = betas.size() - 1;
  NumericMatrix output(n_eval, 2); /* now return a matrix containing both kernel sums and the kernel derivative sums */
  double denom;
  double exp_mult;
  NumericMatrix Ly(ord + 1, n);
  NumericMatrix Ry(ord + 1, n);
  NumericVector tbetas(ord+1);
  for(int k=0; k<ord; k++) tbetas[k] = (k + 1) * betas[k + 1] - betas[k];
  tbetas[ord] = -betas[ord];
  for(int i = 0; i <= ord; i++) Ly(i, 0) = pow(-x[0], i) * y[0];
  for(int i = 1; i < n; i++){
    for(int j = 0; j <= ord; j++){
      Ly(j, i) = pow(-x[i], j) * y[i] + exp((x[i - 1] - x[i]) / h) * Ly(j, i - 1);
      Ry(j, n - i - 1) = exp((x[n - i - 1] - x[n - i]) / h) * (pow(x[n - i], j) * y[n - i] + Ry(j, n - i));
    }
  }
  if(Counts.size() == 0){
    int count;
    NumericVector counts(n_eval);
    count = 0;
    for(int i = 0; i < n_eval; i++){
      if(x_eval[i] >= x[n - 1]){
        counts[i] = n;
      }
      else{
        while(count < n && x[count] <= x_eval[i]){
          count += 1;
        }
        counts[i] = count;
      }
    }
    for(int orddo = 0; orddo <= ord; orddo++){
      NumericVector coefs(orddo + 1);
      coefs[0] = coefs[orddo] = 1;
      if(orddo > 1){
        double num = 1;
        for(int j = 2; j <= orddo; j++) num *= j;
        double denom1 = 1;
        double denom2 = num / orddo;
        for(int i = 2; i <= orddo; i++){
          coefs[i - 1] = num / denom1 / denom2;
          denom1 *= i;
          denom2 /= (orddo - i + 1);
        }
      }
      denom = pow(h, orddo);
      int ix;
      for(int i = 0; i < n_eval; i++){
        ix = round(counts[i]);
        if(ix == 0){
          exp_mult = exp((x_eval[i] - x[0]) / h);
          output(i, 0) += betas[orddo] * pow(x[0] - x_eval[i], orddo) / denom * exp_mult * y[0];
          output(i, 1) += tbetas[orddo] * pow(x[0] - x_eval[i], orddo) / denom * exp_mult * y[0];
          for(int j = 0; j <= orddo; j++){
            output(i, 0) += betas[orddo] * coefs[j] * pow(-x_eval[i], orddo - j) * Ry(j, 0) / denom * exp_mult;
            output(i, 1) += tbetas[orddo] * coefs[j] * pow(-x_eval[i], orddo - j) * Ry(j, 0) / denom * exp_mult;
          }
        }
        else{
          exp_mult = exp((x[ix - 1] - x_eval[i]) / h);
          for(int j = 0; j <= orddo; j++){
            output(i, 0) += betas[orddo] * coefs[j] * (pow(x_eval[i], orddo - j) * Ly(j, ix - 1) * exp_mult + pow(-x_eval[i], orddo - j) * Ry(j, ix - 1) / exp_mult) / denom;
            output(i, 1) -= tbetas[orddo] * coefs[j] * (pow(x_eval[i], orddo - j) * Ly(j, ix - 1) * exp_mult - pow(-x_eval[i], orddo - j) * Ry(j, ix - 1) / exp_mult) / denom;
          }
        }
      }
    }
  }
  else{
    for(int orddo = 0; orddo <= ord; orddo++){
      NumericVector coefs(orddo + 1);
      coefs[0] = coefs[orddo] = 1;
      if(orddo > 1){
        double num = 1;
        for(int j = 2; j <= orddo; j++) num *= j;
        double denom1 = 1;
        double denom2 = num / orddo;
        for(int i = 2; i <= orddo; i++){
          coefs[i - 1] = num / denom1 / denom2;
          denom1 *= i;
          denom2 /= (orddo - i + 1);
        }
      }
      denom = pow(h, orddo);
      int ix;
      for(int i = 0; i < n_eval; i++){
        ix = round(Counts[i]);
        if(ix == 0){
          exp_mult = exp((x_eval[i] - x[0]) / h);
          output(i, 0) += betas[orddo] * pow(x[0] - x_eval[i], orddo) / denom * exp_mult * y[0];
          output(i, 1) += tbetas[orddo] * pow(x[0] - x_eval[i], orddo) / denom * exp_mult * y[0];
          for(int j = 0; j <= orddo; j++){
            output(i, 0) += betas[orddo] * coefs[j] * pow(-x_eval[i], orddo - j) * Ry(j, 0) / denom * exp_mult;
            output(i, 1) += tbetas[orddo] * coefs[j] * pow(-x_eval[i], orddo - j) * Ry(j, 0) / denom * exp_mult;
          }
        }
        else{
          exp_mult = exp((x[ix - 1] - x_eval[i]) / h);
          for(int j = 0; j <= orddo; j++){
            output(i, 0) += betas[orddo] * coefs[j] * (pow(x_eval[i], orddo - j) * Ly(j, ix - 1) * exp_mult + pow(-x_eval[i], orddo - j) * Ry(j, ix - 1) / exp_mult) / denom;
            output(i, 1) -= tbetas[orddo] * coefs[j] * (pow(x_eval[i], orddo - j) * Ly(j, ix - 1) * exp_mult - pow(-x_eval[i], orddo - j) * Ry(j, ix - 1) / exp_mult) / denom;
          }
        }
      }
    }
  }
  return output;
}
