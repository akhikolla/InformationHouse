#include <Rcpp.h>
using namespace Rcpp;

//' Fast binning of numeric vector into equidistant bins
//' 
//' Fast binning of numeric vector into equidistant bins
//'
//' Missing values (NA, Inf, NaN) are added at the end of the vector as the last bin returned if missinglast is set to TRUE
//' 
//' @param x A matrix of regressor variables. Must have the same number of rows as the length of y.
//' @param width The width of the bins
//' @param origin The starting point for the bins. Any number smaller than origin will be disregarded
//' @param missinglast Boolean. Should the missing observations be added as a separate element at the end of the returned count vector.
//' @return An list with elements counts (the frequencies), origin (the origin), width (the width), missing (the number of missings), and last_bin_is_missing (boolean) telling whether the missinglast is true or not.
//' @author Hadley Wickham (from SO: https://stackoverflow.com/questions/13661065/superimpose-histogram-fits-in-one-plot-ggplot) - adapted here by Claus Ekstrøm <claus@@rprimer.dk>
//' @examples
//'
//' set.seed(1)
//' x <- sample(10, 20, replace = TRUE)
//' bin(x, 15)
//' 
//' @export
// [[Rcpp::export]]
List bin(NumericVector x, double width, double origin = 0, bool missinglast=false) {
  int bin, nmissing = 0;
  std::vector<int> out;

  if (width<=0)
    stop("width must be positive");
    
  NumericVector::iterator x_it = x.begin();
  for(; x_it != x.end(); ++x_it) {
    double val = *x_it;
    if (ISNAN(val)) {
      ++nmissing;
    } else {      
      if (val<origin)
	continue;
      
      bin = (val - origin) / width;

      //     Rcout << val << "  " << bin << std::endl;
      
      // Make sure there\'s enough space
      if (bin >= out.size()) {
	out.resize(bin + 1);
      }
      ++out[bin];
    }
  }
  
  // Put missing values in the last position
  if (missinglast)
    out.push_back(nmissing);
  //  return out;

  Rcpp::List RVAL =  Rcpp::List::create(Rcpp::Named("counts") = out,
					Rcpp::Named("origin") = origin,
					Rcpp::Named("width") = width,
					Rcpp::Named("missing") = nmissing,
					Rcpp::Named("last_bin_is_missing") = missinglast
					);

  return(RVAL);

}

