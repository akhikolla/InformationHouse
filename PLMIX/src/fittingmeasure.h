struct Comp{
  Comp(const Rcpp::NumericVector& v ) : _v(v) {}
  bool operator ()(int a, int b) { return _v[a] < _v[b]; }
  const Rcpp::NumericVector& _v;
};

Rcpp::IntegerMatrix top1freq1dim(Rcpp::IntegerMatrix pi_inv) ;
Rcpp::IntegerMatrix tau(Rcpp::IntegerMatrix pi_inv) ;

Rcpp::IntegerMatrix PLMIXsim(int N, int K, int G, Rcpp::NumericMatrix p, Rcpp::NumericMatrix ref_order,Rcpp::NumericVector weights, bool rankingFormat,  Rcpp::IntegerMatrix pi_inv) ;

Rcpp::IntegerMatrix PLMIXsim1dim(int N, int K, int G, Rcpp::NumericMatrix p, Rcpp::NumericMatrix ref_order, Rcpp::NumericVector weights, bool rankingFormat,  Rcpp::IntegerMatrix pi_inv) ;

Rcpp::IntegerMatrix PLMIXsimnotronc1dim(int N, int K, int G, Rcpp::NumericMatrix p, Rcpp::NumericMatrix ref_order, Rcpp::NumericVector weights, bool rankingFormat) ;

Rcpp::IntegerMatrix PLMIXsimnotronc(int N, int K, int G, Rcpp::NumericMatrix p, Rcpp::NumericMatrix ref_order, Rcpp::NumericVector weights, bool rankingFormat) ;


Rcpp::IntegerVector quickintsample(int n, int size, Rcpp::NumericVector prob) ;

Rcpp::IntegerVector quickintsample1dim(int n, int size, Rcpp::NumericVector prob) ; 
