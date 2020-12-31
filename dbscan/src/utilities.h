#include <Rcpp.h>
using namespace Rcpp;

// std::to_string is apparently a c++11 only thing that crashes appveyor, so using ostringstream it is!
namespace patch
{
template < typename T > std::string to_string( const T& n )
{
  std::ostringstream stm ;
  stm << n ;
  return stm.str() ;
}
}

template <typename T, typename C> bool contains (const T& container, const C& key)
{
  if (std::find(container.begin(), container.end(), key) != container.end()){
    return true; 
  } else {
    return false; 
  }
}

// [[Rcpp::export]]
IntegerVector lowerTri(IntegerMatrix m){
  int n = m.nrow(); 
  IntegerVector lower_tri = IntegerVector(n * (n - 1) / 2);
  for (int i = 0, c = 0; i < n; ++i){
    for (int j = i + 1; j < n; ++j){
      if (i < j) lower_tri[c++] = m(i, j);
    }
  }
  return(lower_tri);
}

IntegerVector which_cpp( NumericVector x, double value) {
  int nx = x.size();
  std::vector<int> y;
  y.reserve(nx);
  for(int i = 0; i < nx; ++i) { if (x[i] == value) y.push_back(i); }
  return wrap(y);
}

IntegerVector which_cpp( IntegerVector x, int value) {
  int nx = x.size();
  std::vector<int> y;
  y.reserve(nx);
  for(int i = 0; i < nx; ++i) { if (x[i] == value) y.push_back(i); }
  return wrap(y);
}

IntegerVector which_geq( IntegerVector x, int value) {
  int nx = x.size();
  std::vector<int> y;
  y.reserve(nx);
  for(int i = 0; i < nx; ++i) { if (x[i] >= value) y.push_back(i); }
  return wrap(y);
}


// [[Rcpp::export]]
NumericVector combine(const NumericVector& t1, const NumericVector& t2){
  std::size_t n = t1.size() + t2.size();
  NumericVector output = Rcpp::no_init(n);
  std::copy(t1.begin(), t1.end(), output.begin());
  std::copy(t2.begin(), t2.end(), output.begin()+t1.size());
  return output;
}

IntegerVector combine(const IntegerVector& t1, const IntegerVector& t2){
  std::size_t n = t1.size() + t2.size();
  IntegerVector output = Rcpp::no_init(n);
  std::copy(t1.begin(), t1.end(), output.begin());
  std::copy(t2.begin(), t2.end(), output.begin()+t1.size());
  return output;
}

// Faster version of above combine function, assuming you can precompute and store
// the containers needing to be concatenated
// [[Rcpp::export]]
IntegerVector concat_int (List const& container){
  int total_length = 0;
  for (List::const_iterator it = container.begin(); it != container.end(); ++it){
      total_length += as<IntegerVector>(*it).size(); 
  }
  int pos = 0;
  IntegerVector output = Rcpp::no_init(total_length);
  for (List::const_iterator it = container.begin(); it != container.end(); ++it){
    IntegerVector vec = as<IntegerVector>(*it); 
    std::copy(vec.begin(), vec.end(), output.begin() + pos);
    pos += vec.size();
  }
  return(output);
}

