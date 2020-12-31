#include <cmath>
#include <Rcpp.h>
#include <time.h>
#include <algorithm>

using namespace Rcpp;
using namespace std;


// [[Rcpp::export]]
void negative_check(NumericVector dat){

  for(int i=0; i<dat.size(); i++){
    if(dat[i] <= 0 || Rcpp::internal::Rcpp_IsNA(dat[i]) == true){
      stop("Data may not contain gative values, zeros or NAs.");
    }
  }

}


//' Estimating the Shape Parameter by Method of Maximum Likelihood (MLE)
//'
//' This function can be used to estimate the \code{shape} parameter using the 
//' Maximum Likelihood Estimator method (Newman 2005). It can be used to 
//' obtain biased and unbiased estimates of the shape and scale parameters as 
//' well as the confidence interval for the shape parameter for the biased estimates.
//'
//' @param dat vector of observations
//' @param biased TRUE/FALSE to indicate biased or unbiased estimates
//' @param significance level of significance
//'
//' @return A list of the following form:
//' \describe{
//'   \item{shape}{Estimate of the shape parameter of the data}
//'   \item{lower_bound}{Upper error bound of the estimate of shape}
//'   \item{upper_bound}{Lower error bound of the estimate of shape}
//'   \item{scale}{Estimate of the scale parameter of the data (which is taken
//'   to be the minimum of the data)}
//' }
//' 
//' @references
//' Newman MEJ (2005). "Power Laws, Pareto Distributions And Zipf's Law." 
//' Contemporary Physics, 46, 323-351.
//' 
//' @examples
//' x <- generate_pareto(10000, 5, 2)
//' alpha_mle(x, TRUE, 0.05)
//'
//' x <- generate_pareto(10000, 5, 2)
//' alpha_mle(x, FALSE)
//' 
//' @export
// [[Rcpp::export]]
List alpha_mle(NumericVector dat, bool biased = true,
               Rcpp::Nullable<double> significance = R_NilValue){

  negative_check(dat);

  NumericVector::iterator x_min = min_element(dat.begin(), dat.end());
  double r_x_min = *x_min;

  double n = dat.size();
  double log_sum = 0.0;

  for (int i = 0; i < n; i++){
    log_sum += log(dat[i]);
  }

  if (!significance.isNull()){
    double estimate = n / (log_sum - (n * log(*x_min)));

    double z = R::qnorm(1.0 - (Rcpp::as<double>(significance) / 2.0), 0, 1,
                        TRUE, FALSE);
    double err = z * (sqrt(n + 1.0) / (log_sum - n * log(*x_min)));
    double upper_bound = estimate + err;
    double lower_bound = estimate - err;

    return List::create(Named("shape") = estimate,
                        Named("lower_bound") = lower_bound,
                        Named("upper_bound") = upper_bound,
                        Named("scale") = r_x_min);

  }else{
    if(biased == 1){
      return List::create(Named("shape") = (n / (log_sum - (n * log(*x_min)))),
                          Named("scale") = r_x_min);
    }else{
      double shape = n / (log_sum - (n * log(*x_min)));
      return List::create(Named("shape") = ((n - 2.0) / n) * shape,
                          Named("scale") =
                          r_x_min*(1 - (1 / ((n - 1) * shape))));
    }
  }

}



//' Estimating the Shape Parameter by Weighted Least Squares Method (WLS)
//'
//' This function uses the Weighted Least Squares Method (WLS) to estimate the
//' \code{shape} parameter of a given set of data. (Nair et al. 2019)
//'
//' @param dat vector of observations
//'
//' @return A list of the following form:
//' \describe{
//'   \item{shape}{Estimate of the shape parameter of the data}
//'   \item{scale}{Estimate of the scale parameter of the data (which is taken
//'   to be the minimum of the data)}
//' }
//' 
//' @references
//' Nair J, Wierman A, Zwart B (2019). "The Fundamentals Of Heavy Tails: 
//' Properties, Emergence, And Identification." 
//' http://users.cms.caltech.edu/ adamw/heavytails.html.
//' 
//' @examples
//' x <- generate_pareto(10000, 5, 2)
//' alpha_percentile(x)
//' 
//' @export
// [[Rcpp::export]]
List alpha_wls(NumericVector dat){

  negative_check(dat);

  double n = dat.size();

  NumericVector::iterator x_min = min_element(dat.begin(), dat.end());
  double r_x_min = *x_min;

  double log_sum = 0.0;
  double log_numerator = 0.0;

  for(int i = 0; i < n; i++)
  {
    log_sum += log(dat[i]);
  }

  log_sum -= n* log(*x_min);

  NumericVector y(dat.size());

  sort(dat.begin(), dat.end());

  for(int i = 0; i < n; i++)
  {
    if(i!=0){
      if(dat[i]==dat[i-1]){
        NumericVector::iterator location = find(dat.begin(), dat.begin()+i, dat[i]);
        y[i] = n-distance(dat.begin(), location);
      }else{
        y[i] = n-i;
      }
    }else{
      y[i] = n-i;
    }
  }

  for(int i = 0; i < n; i++)
  {
    log_numerator += log(y[i]);
  }

  log_numerator -= n*log(n);

  return List::create(Named("shape")=(-(log_numerator/log_sum)),
                      Named("scale")=r_x_min);
}



//' Estimating the Shape Parameter by Hill's Estimator
//'
//' This function uses the Hill's Estimator to estimate the shape parameter
//' of a given set of data. (Nair et al. 2019; Pokorna 2016; Hill 1975)
//' It is especially useful when the data is known not to follow an exact 
//' Pareto distribution but the tail of the data does. Thus, the specification 
//' of \code{k}, the \code{k}th largest observation, allows to specify the 
//' point from where Pareto-like behavior may be seen. It is also possible to 
//' specify the value at which the tail begins.
//' When \code{k=n}, the Hill's Estimator returns the same estimate as 
//' \code{alpha_mle} with a warning notifying the user.
//' 
//'
//' @param dat vector of observations
//' @param k number of observations / value equal to or greater than to
//' consider for tail
//' @param value (TRUE/FALSE) indicating if the value which is specified
//' in "k" (TRUE)
//'
//' @return A list of the following form:
//' \describe{
//'   \item{shape}{Estimate of the shape parameter of the data}
//'   \item{scale}{Estimate of the scale parameter of the data (which is taken
//'   to be the minimum of the data)}
//' }
//' 
//' @references
//' Nair J, Wierman A, Zwart B (2019). "The Fundamentals Of Heavy Tails: 
//' Properties, Emergence, And Identification."
//' http://users.cms.caltech.edu/ adamw/heavytails.html.
//' 
//' Pokorna M (2016). Estimation and Application of the Tail Index. 
//' Bachelor's thesis, Charles University in Prague, Faculty of Social 
//' Sciences, Institute of Economic Studies.
//' 
//' Hill B (1975). "A Simple General Approach To Inference About The Tail 
//' Of A Distribution."The Annals of Statistics, 3(5), 1163-1174.
//' 
//' @examples
//' x <- generate_pareto(10000, 5, 2)
//' alpha_hills(x, 400)
//' 
//' @export
// [[Rcpp::export]]
List alpha_hills(NumericVector dat, double k, bool value = false){

  negative_check(dat);

  if(k==dat.size()){
    Rcpp::warning("Setting k as the number of observations makes it equivalent to the MLE (alpha_mle function).");
  }else if(k>dat.size() && value==false){
    stop("k cannot be larger than the size of the data.");
  }

  double double_k = (double)k;
  double r_x_min = 0;
  double shape = 0;

  double denominator = 0.0;

  if(value == 1){

    dat = dat[dat >= k];

    if(dat.size()==0){
      stop("There are no values greater than or equal to specified k.");
    }

    for(int i=0; i<dat.size(); i++){
      denominator += log(dat[i]);
    }

    NumericVector::iterator x_min = min_element(dat.begin(), dat.end());
    r_x_min = *x_min;

    denominator -= ((double)dat.size()*log(r_x_min));
    shape = ((double)dat.size()/denominator);

  }else{

    sort(dat.begin(), dat.end(), std::greater<double>());

    for(int i=0; i<k; i++){
      denominator += log(dat[i]);
    }

    denominator -= (double_k*log(dat[k-1]));
    shape = (double_k/denominator);

    r_x_min = dat[k-1];

  }

  return List::create(Named("shape")=shape,
                      Named("scale")=r_x_min);
}



//' Estimating the Shape Parameter by Method of Moments
//'
//' This function uses the Method of Moments to estimate the shape parameter
//' of a given set of data. (Rytgaard 1990) The method of moments is only 
//' accurate if \eqn{\alpha} (\code{shape parameter}) is greater than or equal 
//' to 1 (Brazauskas and Serfling 2000). This function issues a warning if it 
//' detects that \eqn{\alpha} may be less than 1.
//'
//' @param dat vector of observations
//'
//' @return A list of the following form:
//' \describe{
//'   \item{shape}{Estimate of the shape parameter of the data}
//'   \item{scale}{Estimate of the scale parameter of the data (which is taken
//'   to be the minimum of the data)}
//' }
//' 
//' @references
//' Rytgaard M (1990). "Estimation In The Pareto Distribution." ASTIN
//'  Bulletin: The Journal Of The IAA, 20(2), 201-216.
//'  
//'  Brazauskas V, Serfling R (2000). "Robust and Efficient Estimation Of 
//'  The Tail Index Of A Single-Parameter Pareto Distribution." North 
//'  American Actuarial Journal, 4, 12-27.
//' 
//' @examples
//' x <- generate_pareto(10000, 5, 2)
//' alpha_moment(x)
//' @export
// [[Rcpp::export]]
List alpha_moment(NumericVector dat){

  negative_check(dat);

  if(as<double>(alpha_mle(dat)["shape"])<1.0){
    Rcpp::warning("MLE estimates that this data has a shape parameter less than 1. The Moment Estimator is highly incaccurate for such data. Recommend to use another estimator instead.");
  }
  double n = dat.size();

  NumericVector::iterator x_min = min_element(dat.begin(), dat.end());
  double r_x_min = *x_min;

  double sum = 0.0;

  for(int i = 0; i < n; i++){
    sum += dat[i];
  }

  return List::create(Named("shape")=(sum/(sum- *x_min*n)),
                      Named("scale")=r_x_min);
}



// [[Rcpp::export]]
NumericVector get_Percentiles(NumericVector x, NumericVector p){

  sort(x.begin(), x.end());

  double num = p.size();
  double n = x.size();

  NumericVector out(num);

  for(int i=0; i<num; i++){
    double rank = p[i]*(n+1.0);
    out[i] = x[floor(rank)-1] + ((rank - floor(rank))*(x[floor(rank)]-x[floor(rank)-1]));
  }

  return(out);
}



//' Estimating the Shape Parameter by Method of Percentiles
//'
//' This function uses the Method of Percentiles to estimate the \code{shape}
//' parameter of a given set of data. (Bhatti et al. 2018)
//'
//' @param dat vector of observations
//'
//' @return A list of the following form:
//' \describe{
//'   \item{shape}{Estimate of the shape parameter of the data}
//'   \item{scale}{Estimate of the scale parameter of the data (which is taken
//'   to be the minimum of the data)}
//' }
//' 
//' @references
//' Bhatti SH, Hussain S, Ahmad T, Aslam M, Aftab M, Raza MA (2018). "Efficient
//' estimation of Pareto model: Some modified percentile estimators." 
//' PLoS ONE, 13(5), 1-15.
//' 
//' @examples
//' x <- generate_pareto(10000, 5, 2)
//' alpha_percentile(x)
//' @export
// [[Rcpp::export]]
List alpha_percentile(NumericVector dat){

  negative_check(dat);

  NumericVector::iterator x_min = min_element(dat.begin(), dat.end());

  double r_x_min = *x_min;

  NumericVector percentiles = get_Percentiles(dat, NumericVector::create(0.25, 0.75));

  return List::create(Named("shape")=(log(3)/(log(percentiles[1])-log(percentiles[0]))),
                      Named("scale")=r_x_min);
}



//' Estimating the Shape Parameter by Method of Modified Percentiles
//'
//' This function uses the Method of Modified Percentiles to estimate the 
//' \code{shape} parameter of a given set of data. (Bhatti et al. 2018)
//'
//' @param dat vector of observations
//'
//' @return A list of the following form:
//' \describe{
//'   \item{shape}{Estimate of the shape parameter of the data}
//'   \item{scale}{Estimate of the scale parameter of the data (which is taken
//'   to be the minimum of the data)}
//' }
//' 
//' @references
//' Bhatti SH, Hussain S, Ahmad T, Aslam M, Aftab M, Raza MA (2018). "Efficient
//' estimation of Pareto model: Some modified percentile estimators." 
//' PLoS ONE, 13(5), 1-15.
//' 
//' @examples
//' x <- generate_pareto(10000, 5, 2)
//' alpha_modified_percentile(x)
//' @export
// [[Rcpp::export]]
List alpha_modified_percentile(NumericVector dat){

  negative_check(dat);

  NumericVector::iterator x_min = min_element(dat.begin(), dat.end());

  double r_x_min = *x_min;

  NumericVector percentiles = get_Percentiles(dat, NumericVector::create(0.5, 0.75));

  return List::create(Named("shape")=(log(0.25)+log(2.0))/(log(percentiles[0])-log(percentiles[1])),
                      Named("scale")=r_x_min);

}

// [[Rcpp::export]]
double GeoMean(NumericVector x){

  double totMean = 0;

  for(int i = 0; i<x.size(); i++){
      totMean += log(x[i]);
  }

  totMean = totMean/x.size();

  return(exp(totMean));
}


//' Estimating the Shape Parameter by Geometric Method of Percentiles
//'
//' This function uses the Geometric Method of Percentiles to estimate the
//' \code{shape} parameter of a given set of data. (Bhatti et al. 2018)
//'
//' @param dat vector of observations
//'
//' @return A list of the following form:
//' \describe{
//'   \item{shape}{Estimate of the shape parameter of the data}
//'   \item{scale}{Estimate of the scale parameter of the data (which is taken
//'   to be the minimum of the data)}
//' }
//' 
//' @references
//' Bhatti SH, Hussain S, Ahmad T, Aslam M, Aftab M, Raza MA (2018). "Efficient
//' estimation of Pareto model: Some modified percentile estimators." 
//' PLoS ONE, 13(5), 1-15.
//' 
//' @examples
//' x <- generate_pareto(10000, 5, 2)
//' alpha_geometric_percentile(x)
//' 
//' @export
// [[Rcpp::export]]
List alpha_geometric_percentile(NumericVector dat){

  negative_check(dat);

  NumericVector::iterator x_min = min_element(dat.begin(), dat.end());

  double r_x_min = *x_min;

  NumericVector percentiles = get_Percentiles(dat, NumericVector::create(0.75));

  return List::create(Named("shape")=(log(0.25)+1.0)/(log(GeoMean(dat))-log(percentiles[0])),
                      Named("scale")=r_x_min);

}



//' Estimating the Shape Parameter by Method of Least Squares
//'
//' This function uses the Method of Least Squares to estimate the \code{shape}
//' parameter of a given set of data. (Zaher et al. 2014; Nair et al. 2019)
//'
//' @param dat vector of observations
//'
//' @return A list of the following form:
//' \describe{
//'   \item{shape}{Estimate of the shape parameter of the data}
//'   \item{scale}{Estimate of the scale parameter of the data (which is taken
//'   to be the minimum of the data)}
//' }
//' 
//' @references
//' Zaher HM, El-Sheik AA, El-Magd NATA (2014). "Estimation of Pareto 
//' Parameters Using a Fuzzy Least-Squares Method and Other Known Techniques 
//' with a Comparison." British Journal of Mathematics & Computer Science, 
//' 4(14), 2067-2088.
//' 
//' Nair J, Wierman A, Zwart B (2019). "The Fundamentals Of Heavy Tails: 
//' Properties, Emergence, And Identification."
//' http://users.cms.caltech.edu/ adamw/heavytails.html.
//' 
//' @examples
//' x <- generate_pareto(10000, 5, 2)
//' alpha_ls(x)
//' 
//' @export
// [[Rcpp::export]]
List alpha_ls(NumericVector dat){

  negative_check(dat);

  double n = dat.size();

  NumericVector::iterator x_min = min_element(dat.begin(), dat.end());
  double r_x_min = *x_min;

  NumericVector y_arr(n);
  NumericVector t_arr(n);
  double y_bar = 0.0;
  double t_bar = 0.0;

  for(int i = 0; i < n; i++){
    y_arr[i] = log(1-(i/(n+1)));
    y_bar += y_arr[i];
  }

  y_bar = y_bar/n;

  sort(dat.begin(), dat.end());

  for(int i = 0; i < n; i++){
    t_arr[i] = log(dat[i]);
    t_bar += t_arr[i];
  }

  t_bar = t_bar/n;

  double numerator = 0.0;
  double denominator = 0.0;

  for(int i = 0; i < n; i++){
    numerator += (t_arr[i] - t_bar)*(y_arr[i] - y_bar);
    denominator += pow((t_arr[i] - t_bar), 2.0);
  }

  return List::create(Named("shape")=-(numerator/denominator),
                      Named("scale")=r_x_min);
}

//' Goodness of Fit Test for Pareto Distribution
//'
//' The pareto_test function can be used to identify whether the data is 
//' Pareto distributed (Gulati and Shapiro 2008). The test generates a p-value 
//' corresponding to the actual distribution of the data and is tested for 
//' significance. In the case of Pareto data, the p-value should be greater 
//' than the pre-determined significance level (generally taken as 0.05). 
//'
//' @param dat vector of observations
//'
//' @return A list of the following form:
//' \describe{
//'   \item{p-value}{p-value indicating significance of the test}
//' }
//' 
//' @references
//' Gulati S, Shapiro S (2008). "Goodness-of-Fit Tests for Pareto 
//' Distribution." In F Vonta (ed.), Statistical Models and Methods for 
//' Biomedical and Technical Systems, chapter 19, pp. 259-274. Birkhauser 
//' Basel. ISBN 978-0-8176-4619-6. doi:10.1007/978-0-8176-4619-6.
//' 
//' @examples
//' x <- generate_pareto(10000, 5, 2)
//' pareto_test(x)
//' @export
// [[Rcpp::export]]
List pareto_test(NumericVector dat){
  
  sort(dat.begin(), dat.end());
  
  double n = dat.size();
  
  NumericVector T(n);
  for(int i = 0; i < n; i++){
    T[i] = log(dat[i]/dat[0]);
  }
  
  sort(T.begin(), T.end());
  
  NumericVector Y(n);
  NumericVector U((n-1));
  NumericVector iU((n-1));
  double Yn = 0.0;
  
  for(int i = 0; i < n; i++){
    if(i==0){
      Y[i] = (n - i + 1)*(T[i]);
      U[i] = Y[i];
    }else if(i!=(n-1)){
      Y[i] = (n - i + 1)*(T[i]-T[i-1]);
      U[i] = U[i-1] + Y[i];
    }else{
      Y[i] = (n - i + 1)*(T[i]-T[i-1]);
      Yn = U[i-1] + Y[i];
    }
  }
  
  double U_bar = 0.0;
  double iU_bar = 0.0;
  
  for(int i = 0; i< (n-1); i++){
    U[i] = U[i]/Yn;
    iU[i] = ((double)i+1.0)*U[i];
    U_bar += U[i];
    iU_bar += iU[i];
  }
  
  U_bar = U_bar/(n-1.0);
  iU_bar = iU_bar/(n-1.0);
  
  double Z1 = sqrt(12.0*(n - 1.0)) * (U_bar - 0.5);
  double Z2 = sqrt(5.0*(n - 1.0) / ((n + 2.0) * (n - 2.0)))*( n - 2.0 + 6.0*n*U_bar - (12.0 * iU_bar));
  
  double Z0 = pow(Z1, 2) + pow(Z2, 2);
  
  return List::create(Named("p-value")=exp((-Z0/2.0)));
  
}
