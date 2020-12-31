#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

/* d_abs computes absolute value of a double */

double d_abs(double x){
  if(x > 0) return x;
  else return -x;
}


/* fk_md computes the minimum density hyperplane based on the projected data.
 * That is, if x = Xw is the projected data, then the function finds the value
 * of the lowest (penalised) density estimated from x.
 * Arguments:
 * x = vector of projected points . Should have mean zero!
 * y = vector of weights for the points. In place to make binned approximation possible
 * x_eval = initial set of points at which to evaluate the density
 * h = kernel bandwidth
 * betas = kernel coefficients.
 * al = 2*al*sd(x) is the width of interval which doesn't receive penalty. Encourage solution to be close-ish to mean of x
 * C = penalty coefficient. For details on C, al, see Minimum Density Hyperplanes, Pavlidis et al. (2016)
 */



// [[Rcpp::export]]

double fk_md(NumericVector x, NumericVector y, NumericVector x_eval, double h, NumericVector betas, double al, double C){
  /* setup parameters and object to be returned, etc. */
  int n = x.size();
  int n_eval = x_eval.size();
  int ord = betas.size() - 1;

  /* find local minima in the function using its derivative. tbetas are the coefficients required
   * for kernel derivative sums.
   */
  NumericVector tbetas(ord + 1);
  for(int k = 0; k < ord; k++) tbetas[k] = (k + 1) * betas[k + 1] - betas[k];
  tbetas[ord] = -betas[ord];

  /* miny is value of the minimum density hyperplane. To be updated during search. */
  double miny;
  double stdv; /* to be standard deviation of x */
  NumericVector df(n_eval); /* to store the scaled density derivatives at x_eval */
  if(al > 1e-10){
    /* if al is very small, then assumed location of hyperplane at mean(x) = 0.
     * In this case al is not close to zero and so we search for the minimum density location
     * and the value of the minimum
     */
    double denom;
    double exp_mult;

    /* begin by computing standard deviation of the projected data */
    double var = 0;
    for(int i = 0; i < n; i++) var += x[i] * x[i];
    var /= (n - 1);
    stdv = pow(var, .5);

    /* sort projected data and compute recursive sums as in ksum, dksum, etc. */
    std::sort(x.begin(), x.end());
    NumericMatrix Ly(ord + 1, n);
    NumericMatrix Ry(ord + 1, n);
    for(int i = 0; i <= ord; i++) Ly(i, 0) = pow(-x[0], i) * y[0];
    for(int i = 1; i < n; i++){
      for(int j = 0; j <= ord; j++){
        Ly(j, i) = pow(-x[i], j) * y[i] + exp((x[i - 1] - x[i]) / h) * Ly(j, i - 1);
        Ry(j, n - i - 1) = exp((x[n - i - 1] - x[n - i]) / h) * (pow(x[n - i], j) * y[n - i] + Ry(j, n - i));
      }
    }

    /* determine location of the evaluation points among the projected data */
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

    /* loop over terms in the polynomial term in the kernel and incrementally
     * compute (scaled) density derivative at the evaluation points
     */
    for(int orddo = 0; orddo <= ord; orddo++){ /* orddo is the exponent for the current term in the polynomial */
      /* coefs are the binomial coefficients for the expanded polynomial term */
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

      /* add increment to the estimated density derivatives, as in dksum */
      denom = pow(h, orddo);
      int ix;
      for(int i = 0; i < n_eval; i++){
        ix = round(counts[i]);
        if(ix == 0){
          exp_mult = exp((x_eval[i] - x[0]) / h);
          df[i] -= tbetas[orddo] * pow(x[0] - x_eval[i], orddo) / denom * exp_mult;
          for(int j = 0; j <= orddo; j++) df[i] -= tbetas[orddo] * coefs[j] * pow(-x_eval[i], orddo - j) * Ry(j, 0) / denom * exp_mult;
        }
        else{
          exp_mult = exp((x[ix - 1] - x_eval[i]) / h);
          for(int j = 0; j <= orddo; j++) df[i] += tbetas[orddo] * coefs[j] * (pow(x_eval[i], orddo - j) * Ly(j, ix - 1) * exp_mult - pow(-x_eval[i], orddo - j) * Ry(j, ix - 1) / exp_mult) / denom;
        }
      }
    }

    /* add penalty derivative to those points lying outside (-al*sd, al*sd) */
    for(int i = 0; i < n_eval; i++){
      if(x_eval[i] < (-al * stdv)) df[i] -= 2.0 * C * (-al * stdv - x_eval[i]);
      if(x_eval[i] > (al * stdv)) df[i] += 2.0 * C * (-al * stdv + x_eval[i]);
    }

    /* next search over evaluated derivatives, and use binary search to refine the local minima.
     * f_at_min, df_at_min are the values of the penalised density and its derivative
     * at the minimum being considered. The best minimum found so far is stored in miny, and will end as the
     * lowest density hyperplane to be returned once all minima have been considered.
     */
    double f_at_min, df_at_min, lo, hi, mid; /* lo, hi, mid are the lower/upper/middle points for binary search interval */
    miny = 1e10;
    int pos = 0; /* location index over which to search */
    double eps = 1e-10; /* tolerance for derivative at local minimum, and for width of binary search interval */
    while(pos<(n_eval - 1)){
      if(df[pos] < 0 && df[pos + 1] > 0){
        /* a local minimum between x_eval[pos] and x_eval[pos+1] */
        df_at_min = 1.0;
        lo = x_eval[pos];
        hi = x_eval[pos + 1];
        int ix = round(counts[pos]);
        while((hi - lo) > eps && d_abs(df_at_min) > eps){
          /* refine location of minimum by computing the derivative at
           * the mid-point of the interval. Narrow interval until
           * derivative is close enough to zero
           */
          df_at_min = 0.0;
          f_at_min = 0;
          mid = 0.5 * lo + 0.5 * hi;
          /* compute derivative at midpoint of interval, as for initial evaluation points */
          exp_mult = exp((x[ix - 1] - mid) / h);
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
            for(int j = 0; j <= orddo; j++){
              df_at_min += tbetas[orddo] * coefs[j] * (pow(mid, orddo - j) * Ly(j, ix - 1) * exp_mult-pow(-mid, orddo - j) * Ry(j, ix - 1) / exp_mult) / denom;
              f_at_min += betas[orddo] * coefs[j] * (pow(mid, orddo - j) * Ly(j, ix - 1) * exp_mult + pow(-mid, orddo - j) * Ry(j, ix - 1) / exp_mult) / denom;
            }
          }
          if(mid < (-al * stdv)){
            df_at_min -= 2.0 * C * (-al * stdv - mid);
            f_at_min += C * pow(-al * stdv - mid, 2);
          }
          if(mid > (al * stdv)){
            df_at_min += 2.0 * C * (-al * stdv + mid);
            f_at_min += C * pow(-al * stdv + mid, 2);
          }
          /* update upper and lower bound of the search interval */
          if(df_at_min < (-eps)) lo = mid;
          if(df_at_min > eps) hi = mid;
        }
        /* update best local minimum if the current one is lower than the best found so far */
        if(f_at_min < miny){
          miny = f_at_min;
        }
      }
      pos += 1;
    }
  }
  else{
    /* al is very close to zero and we assume the location of the hyperplane is at mean(x) = 0
     * So compute the density directly
     */
    miny = 0;
    for(int i = 1; i < n; i++){
      double add = 0;
      for(int j = 0; j <= ord; j++) add += betas[j] * pow(d_abs(x[i] / h), j);
      miny += add * exp(-d_abs(x[i]) / h);
    }
  }
  return miny;
}






/* fk_md_dp computes the partial derivatives of the value of the penalised density at
 * the minimum with respect to each of the projected data points.
 * That is, if x = Xw is the projected data, then the function finds the partial
 * derivatives, w.r.t. the elements in x, of the value of the lowest (penalised) density estimated from x.
 * Arguments and interpretation of the function are as in fk_md, except where indicated.
 */


// [[Rcpp::export]]

NumericVector fk_md_dp(NumericVector xo, NumericVector y, NumericVector x_eval, double h, NumericVector betas, double al, double C){
  /* argument xo is the vector of projected points, which must be distinct from its sorted counter-part, x, computed below, since
   * ordering of partial derivatives must match indices of xo = Xw, and not the sorted version.
   */
  int n = xo.size();
  int n_eval = x_eval.size();
  int ord = betas.size() - 1;
  NumericVector tbetas(ord + 1);
  for(int k = 0; k < ord; k++) tbetas[k] = (k + 1) * betas[k + 1] - betas[k];
  tbetas[ord] = -betas[ord];
  double miny;
  double minx = 0; /* location of minimum density hyperplane, i.e., f(minx) = miny. Important for derivatives. */
  double stdv;
  NumericVector df(n_eval);
  if(al > 1e-10){
    double var = 0;
    NumericVector x(n);
    for(int i = 0; i < n; i++){
      x[i] = xo[i];
      var += x[i] * x[i];
    }
    var /= (n - 1);
    stdv = pow(var, .5);
    std::sort(x.begin(), x.end());
    double denom;
    double exp_mult;
    NumericMatrix Ly(ord + 1, n);
    NumericMatrix Ry(ord + 1, n);
    for(int i = 0; i <= ord; i++) Ly(i, 0) = pow(-x[0], i) * y[0];
    for(int i = 1; i < n; i++){
      for(int j = 0; j <= ord; j++){
        Ly(j, i) = pow(-x[i], j) * y[i] + exp((x[i - 1] - x[i]) / h) * Ly(j, i - 1);
        Ry(j, n - i - 1) = exp((x[n - i - 1] - x[n - i]) / h) * (pow(x[n - i], j) * y[n - i] + Ry(j, n - i));
      }
    }
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
        if(ix == 0){
          exp_mult = exp((x_eval[i] - x[0]) / h);
          df[i] -= tbetas[orddo] * pow(x[0] - x_eval[i], orddo) / denom * exp_mult;
          for(int j = 0; j <= orddo; j++) df[i] -= tbetas[orddo] * coefs[j] * pow(-x_eval[i], orddo - j) * Ry(j, 0) / denom * exp_mult;
        }
        else{
          exp_mult = exp((x[ix - 1] - x_eval[i]) / h);
          for(int j = 0; j <= orddo; j++) df[i] += tbetas[orddo] * coefs[j] * (pow(x_eval[i], orddo - j) * Ly(j, ix - 1) * exp_mult - pow(-x_eval[i], orddo - j) * Ry(j, ix - 1) / exp_mult) / denom;
        }
      }
    }
    for(int i = 0; i < n_eval; i++){
      if(x_eval[i] < (-al * stdv)) df[i] -= 2.0 * C * (-al * stdv - x_eval[i]);
      if(x_eval[i]>(al * stdv)) df[i] += 2.0 * C * (-al * stdv + x_eval[i]);
    }
    double f_at_min, df_at_min, lo, hi, mid;
    miny = 1e10;
    int pos = 0;
    double eps = 1e-10;
    while(pos < (n_eval - 1)){
      if(df[pos] < 0 && df[pos + 1] > 0){
        df_at_min = 1.0;
        lo = x_eval[pos];
        hi = x_eval[pos + 1];
        int ix = round(counts[pos]);
        while((hi - lo) > eps && d_abs(df_at_min) > eps){
          df_at_min = 0.0;
          f_at_min = 0;
          mid = 0.5 * lo + 0.5 * hi;
          exp_mult = exp((x[ix - 1] - mid) / h);
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
            for(int j = 0; j <= orddo; j++){
              df_at_min += tbetas[orddo] * coefs[j] * (pow(mid, orddo - j) * Ly(j, ix - 1) * exp_mult - pow(-mid, orddo - j) * Ry(j, ix - 1) / exp_mult) / denom;
              f_at_min += betas[orddo] * coefs[j] * (pow(mid, orddo - j) * Ly(j, ix - 1) * exp_mult + pow(-mid, orddo - j) * Ry(j, ix - 1) / exp_mult) / denom;
            }
          }
          if(mid<(-al * stdv)){
            df_at_min -= 2.0 * C * (-al * stdv - mid);
            f_at_min += C * pow(-al * stdv - mid, 2);
          }
          if(mid > (al * stdv)){
            df_at_min += 2.0 * C * (-al * stdv + mid);
            f_at_min += C * pow(-al * stdv + mid, 2);
          }
          if(df_at_min < (-eps)) lo = mid;
          if(df_at_min > eps) hi = mid;
        }
        if(f_at_min < miny){
          miny = f_at_min;
          minx = mid;
        }
      }
      pos += 1;
    }
  }
  else{
    minx = 0; /* if al is close to zero then minx = 0 and miny does not matter */
    stdv = 1; /* stdev given a dummy value to avoid unitialised use in later check. Operation of check when minx = 0 is the same for all stdv > 0 */
  }
  /* dp stores the vector of partial derivatives to be returned */
  NumericVector dp(n);
  double add;
  /* since minx is known, fast kernel computations are not needed, and so partial derivatives
   * are computed directly.
   */
  for(int i = 0; i < n; i++){
    add = 0;
    if(xo[i] > minx){ /* derivative of density has different sign depending on which side the point lies of minx */
      for(int j = 0; j <= ord; j++) add += tbetas[j] * pow((xo[i] - minx) / h, j);
      dp[i] = add * exp((minx - xo[i]) / h);
    }
    else{
      for(int j = 0; j <= ord; j++) add -= tbetas[j] * pow((minx - xo[i]) / h, j);
      dp[i] = add * exp((xo[i] - minx) / h);
    }
  }
  /* if minx lies outside (-al*sd, al*sd) then derivative of penalty must be added */
  if(minx < (-al * stdv)){
    double cnst = 2.0 * al * C / stdv / (n - 1.0) * (minx + al * stdv);
    for(int i = 0; i < n; i++) dp[i] += cnst * xo[i];
  }
  if(minx > (al * stdv)){
    double cnst = 2.0 * al * C / stdv / (n - 1.0) * (al * stdv - minx);
    for(int i = 0; i < n; i++) dp[i] += cnst * xo[i];
  }
  return(dp);
}






/* fk_is_minim_md determines if the current minimum density hyperplane is valid,
 * in the sense that it lies between modes of the projected density.
 * Arguments and interpretation of the function are as in fk_md/fk_md_dp, except where indicated.
 */


// [[Rcpp::export]]

double fk_is_minim_md(NumericVector x, NumericVector y, NumericVector x_eval, double h, NumericVector betas, double al, double C){
  int n = x.size();
  int n_eval = x_eval.size();
  int ord = betas.size() - 1;
  NumericVector tbetas(ord + 1);
  for(int k = 0; k < ord; k++) tbetas[k] = (k + 1) * betas[k + 1] - betas[k];
  tbetas[ord] = -betas[ord];
  double miny;
  double minx = 0;
  double stdv;
  /* mode1 and modef are the leftmost and rightmost modes of the density estimated from x */
  double modef;
  double mode1 = 1e10; /* mode1 is first set to a large number so that it is easy to determine if it has yet been set during search */

  /* derivative of density and of penalty are separated since their sum determines the location of
   * the minimum, minx, but we care if it lies between local maxima of the density only.
   */
  NumericVector df(n_eval);
  NumericVector dfpen(n_eval);
  if(al > 1e-10){
    double var = 0;
    for(int i = 0; i < n; i++) var += x[i] * x[i];
    var /= (n - 1);
    stdv = pow(var, .5);
    std::sort(x.begin(), x.end());
    double denom;
    double exp_mult;
    NumericMatrix Ly(ord + 1, n);
    NumericMatrix Ry(ord + 1, n);
    for(int i = 0; i <= ord; i++) Ly(i, 0) = pow(-x[0], i) * y[0];
    for(int i = 1; i < n; i++){
      for(int j = 0; j <= ord; j++){
        Ly(j, i) = pow(-x[i], j) * y[i] + exp((x[i - 1] - x[i]) / h) * Ly(j, i - 1);
        Ry(j, n - i - 1) = exp((x[n - i - 1] - x[n - i]) / h) * (pow(x[n - i], j) * y[n - i] + Ry(j, n - i));
      }
    }
    int count = 0;
    NumericVector counts(n_eval);
    for(int i = 0; i < n_eval; i++){
      if(x_eval[i] >= x[n - 1]){
        counts[i] = n;
      }
      else{
        while(x[count] <= x_eval[i]){
          count += 1;
          if(count >= n) break;
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
          df[i] -= tbetas[orddo] * pow(x[0] - x_eval[i], orddo) / denom * exp_mult;
          for(int j = 0; j <= orddo; j++) df[i] -= tbetas[orddo] * coefs[j] * pow(-x_eval[i], orddo - j) * Ry(j, 0) / denom * exp_mult;
        }
        else{
          exp_mult = exp((x[ix - 1] - x_eval[i]) / h);
          for(int j = 0; j <= orddo; j++) df[i] += tbetas[orddo] * coefs[j] * (pow(x_eval[i], orddo - j) * Ly(j, ix - 1) * exp_mult - pow(-x_eval[i], orddo - j) * Ry(j, ix - 1) / exp_mult) / denom;
        }
      }
    }
    for(int i = 0; i < n_eval; i++){
      if(x_eval[i] < (-al * stdv)) dfpen[i] -= 2.0 * C * (-al * stdv - x_eval[i]);
      if(x_eval[i] > (al * stdv)) dfpen[i] += 2.0 * C * (-al * stdv + x_eval[i]);
    }
    double f_at_min, df_at_min, lo, hi, mid;
    miny = 1e10;
    int pos = 0;
    double eps = 1e-10;
    /* when looping over x_eval to find the location of minima and maxima, we also
     * store the locations of the first and final mode (local maximum) of the density.
     */
    while(pos < (n_eval - 1)){
      if(df[pos] > 0 && df[pos + 1] < 0){
        /* a mode is found, and if mode1 has not yet been set (i.e, is still a large number)
         * then it is set to the current location. The final mode updates every time a new
         * mode is found.
         */
        if(x_eval[pos] < mode1){
          mode1 = x_eval[pos];
        }
        modef = x_eval[pos];
      }
      if((df[pos] + dfpen[pos]) < 0 && (df[pos + 1] + dfpen[pos + 1]) > 0){
        df_at_min = 1.0;
        lo = x_eval[pos];
        hi = x_eval[pos + 1];
        int ix = round(counts[pos]);
        while((hi - lo) > eps && d_abs(df_at_min) > eps){
          df_at_min = 0.0;
          f_at_min = 0;
          mid = 0.5 * lo + 0.5 * hi;
          exp_mult = exp((x[ix - 1] - mid) / h);
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
            for(int j = 0; j <= orddo; j++){
              df_at_min += tbetas[orddo] * coefs[j] * (pow(mid, orddo - j) * Ly(j, ix - 1) * exp_mult - pow(-mid, orddo - j) * Ry(j, ix - 1) / exp_mult) / denom;
              f_at_min += betas[orddo] * coefs[j] * (pow(mid, orddo - j) * Ly(j, ix - 1) * exp_mult + pow(-mid, orddo - j) * Ry(j, ix - 1) / exp_mult) / denom;
            }
          }
          if(mid < (-al * stdv)){
            df_at_min -= 2.0 * C * (-al * stdv - mid);
            f_at_min += C * pow(-al * stdv - mid, 2);
          }
          if(mid > (al * stdv)){
            df_at_min += 2.0 * C * (-al * stdv + mid);
            f_at_min += C * pow(-al * stdv + mid, 2);
          }
          if(df_at_min < (-eps)) lo = mid;
          if(df_at_min > eps) hi = mid;
        }
        if(f_at_min < miny){
          miny = f_at_min;
          minx = mid;
        }
      }
      pos += 1;
    }
  }
  else minx = 0;
  /* if minx lies between mode1 and modef then return 1, otherwise 0. */
  double ret = 0;
  if(minx > mode1 && minx < modef) ret = 1;
  return(ret);
}




/* fk_md_b determines the location of the minimum density hyperplane, i.e.,
 * the value minx s.t. f(minx) = miny.
 * Arguments and interpretation of the function are as in fk_md, except where indicated.
 */



// [[Rcpp::export]]

double fk_md_b(NumericVector x, NumericVector y, NumericVector x_eval, double h, NumericVector betas, double al, double C){
  int n = x.size();
  int n_eval = x_eval.size();
  int ord = betas.size() - 1;
  NumericVector tbetas(ord + 1);
  for(int k = 0; k < ord; k++) tbetas[k] = (k + 1) * betas[k + 1] - betas[k];
  tbetas[ord] = -betas[ord];
  double miny;
  double minx = 0;
  double stdv;
  NumericVector df(n_eval);
  if(al > 1e-10){
    double var = 0;
    for(int i = 0; i < n; i++) var += x[i] * x[i];
    var /= (n - 1);
    stdv = pow(var, .5);
    std::sort(x.begin(), x.end());
    double denom;
    double exp_mult;
    NumericMatrix Ly(ord + 1, n);
    NumericMatrix Ry(ord + 1, n);
    for(int i = 0; i <= ord; i++) Ly(i, 0) = pow(-x[0], i) * y[0];
    for(int i = 1; i < n; i++){
      for(int j = 0; j <= ord; j++){
        Ly(j, i) = pow(-x[i], j) * y[i] + exp((x[i - 1] - x[i]) / h) * Ly(j, i - 1);
        Ry(j, n - i - 1) = exp((x[n - i - 1] - x[n - i]) / h) * (pow(x[n - i], j) * y[n - i] + Ry(j, n - i));
      }
    }
    int count = 0;
    NumericVector counts(n_eval);
    for(int i = 0; i < n_eval; i++){
      if(x_eval[i] >= x[n - 1]){
        counts[i] = n;
      }
      else{
        while(x[count] <= x_eval[i]){
          count += 1;
          if(count >= n) break;
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
          df[i] -= tbetas[orddo] * pow(x[0] - x_eval[i], orddo) / denom * exp_mult;
          for(int j = 0; j <= orddo; j++) df[i] -= tbetas[orddo] * coefs[j] * pow(-x_eval[i], orddo - j) * Ry(j, 0) / denom * exp_mult;
        }
        else{
          exp_mult = exp((x[ix - 1] - x_eval[i]) / h);
          for(int j = 0; j <= orddo; j++) df[i] += tbetas[orddo] * coefs[j] * (pow(x_eval[i], orddo - j) * Ly(j, ix - 1) * exp_mult - pow(-x_eval[i], orddo - j) * Ry(j, ix - 1) / exp_mult) / denom;
        }
      }
    }
    for(int i = 0; i < n_eval; i++){
      if(x_eval[i] < (-al * stdv)) df[i] -= 2.0 * C * (-al * stdv - x_eval[i]);
      if(x_eval[i] > (al * stdv)) df[i] += 2.0 * C * (-al * stdv + x_eval[i]);
    }
    double f_at_min, df_at_min, lo, hi, mid;
    miny = 1e10;
    int pos = 0;
    double eps = 1e-10;
    while(pos < (n_eval - 1)){
      if(df[pos] < 0 && df[pos + 1] > 0){
        df_at_min = 1.0;
        lo = x_eval[pos];
        hi = x_eval[pos + 1];
        int ix = round(counts[pos]);
        while((hi - lo) > eps && d_abs(df_at_min) > eps){
          df_at_min = 0.0;
          f_at_min = 0;
          mid = 0.5 * lo + 0.5 * hi;
          exp_mult = exp((x[ix - 1] - mid) / h);
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
            for(int j = 0; j <= orddo; j++){
              df_at_min += tbetas[orddo] * coefs[j] * (pow(mid, orddo - j) * Ly(j, ix - 1) * exp_mult - pow(-mid, orddo - j) * Ry(j, ix - 1) / exp_mult) / denom;
              f_at_min += betas[orddo] * coefs[j] * (pow(mid, orddo - j) * Ly(j, ix - 1) * exp_mult + pow(-mid, orddo - j) * Ry(j, ix - 1) / exp_mult) / denom;
            }
          }
          if(mid < (-al * stdv)){
            df_at_min -= 2.0 * C * (-al * stdv - mid);
            f_at_min += C * pow(-al * stdv - mid, 2);
          }
          if(mid > (al * stdv)){
            df_at_min += 2.0 * C * (-al * stdv + mid);
            f_at_min += C * pow(-al * stdv + mid, 2);
          }
          if(df_at_min < (-eps)) lo = mid;
          if(df_at_min > eps) hi = mid;
        }
        if(f_at_min < miny){
          miny = f_at_min;
          minx = mid;
        }
      }
      pos += 1;
    }
  }
  else{
    for(int i = 1; i < n; i++){
      double add = 0;
      for(int j = 0; j <= ord; j++) add += betas[j] * pow(d_abs(x[i] / h), j);
      miny += add * exp(-d_abs(x[i]) / h);
    }
  }
  return minx;
}
