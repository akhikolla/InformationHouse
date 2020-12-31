#include <Rcpp.h>
using namespace Rcpp;
//' Calculate the haversine distance in kilometers given lat/lon pairs
//'
//' @param lat1 A vector of latitudes
//' @param lon1 A vector of longitudes
//' @param lat2 A vector of latitudes
//' @param lon2 A vector of longitudes
//' @return a vector of distances in kilometers
//' @export
//' @examples
//' # simple haversine calculation 
//' lon1 <- runif(-160, -60, n = 10e6)
//' lat1 <- runif(40, 60, n = 10e6)
//' lon2 <- runif(-160, -60, n = 10e6)
//' lat2 <- runif(40, 60, n = 10e6)
//' df <- data.frame(lat1, lon1, lat2, lon2)
//' df$havers <- haversine(df$lat1, df$lon1, df$lat2, df$lon2)
//[[Rcpp::export]]
Rcpp::NumericVector haversine(Rcpp::NumericVector lat1, Rcpp::NumericVector lon1, Rcpp::NumericVector lat2, Rcpp::NumericVector lon2){
    double halfC = M_PI / 180;
    // convert lat and lons to radians
    Rcpp::NumericVector lat1r = lat1 * halfC;
    Rcpp::NumericVector lon1r = lon1 * halfC;
    Rcpp::NumericVector lat2r = lat2 * halfC;
    Rcpp::NumericVector lon2r = lon2 * halfC;
    // haversine formula
    Rcpp::NumericVector dlon = lon2r - lon1r;
    Rcpp::NumericVector dlat = lat2r - lat1r; 
    Rcpp::NumericVector a = pow(sin(dlat/2), 2) + cos(lat1r) * cos(lat2r) * pow(sin(dlon/2), 2);
    Rcpp::NumericVector c = 2 * asin(sqrt(a)); 
    double r = 6371; // radius of earth in km
    Rcpp::NumericVector d = c * r;
    Rcpp::NumericVector haversine = d; // distance in km
    return haversine;
}
