#include <Rcpp.h>
#include <cmath>
#define pi 3.14159265358979323846
#define earthRadiusKm 6378137 // 6371.0
using namespace Rcpp;

//Found much of this source code on stackoverflow

// [[Rcpp::export]] // This function converts decimal degrees to radians
double deg2rad(double deg) {
  return (deg * pi / 180);
}
// [[Rcpp::export]] //  This function converts radians to decimal degrees
double rad2deg(double rad) {
  return (rad * 180 / pi);
}
// distanceEartch calculates Haversine distance
// [[Rcpp::export]] 
double distanceEarth(double lat1d, double lon1d, double lat2d, double lon2d) {
  double lat1r, lon1r, lat2r, lon2r, u, v;
  lat1r = deg2rad(lat1d);
  lon1r = deg2rad(lon1d);
  lat2r = deg2rad(lat2d);
  lon2r = deg2rad(lon2d);
  u = sin((lat2r - lat1r)/2);
  v = sin((lon2r - lon1r)/2);
  return 2.0 * earthRadiusKm * asin(sqrt(u * u + cos(lat1r) * cos(lat2r) * v * v));
}

// [[Rcpp::export]]
int vecminInd(NumericVector x) { // Calculate position of minimum vector; found code online
  // Rcpp supports STL-style iterators
  NumericVector::iterator it = std::min_element(x.begin(), x.end());
  // we want the value so dereference 
  return it - x.begin();
}

// Double Loop to calculate distance for each person in voter file
// [[Rcpp::export]]
NumericVector nearest_dbox(NumericVector lat1d_vec, NumericVector lon1d_vec, 
                           NumericVector lat2d_vec, NumericVector lon2d_vec) {
  //NA Handling
  Rcpp::StringVector na_error = "Error: NA Exists, Omit NAs from vectors and rerun";
    
  if( any( is_na( lat1d_vec ) ) | any( is_na( lon1d_vec ) ) ) {
    
    Rcout << na_error <<"\n";
    return 0;
    
  } else if (any( is_na( lat2d_vec ) ) | any( is_na( lat2d_vec ) ) ) {
    
    Rcout << na_error <<"\n";
    return 0;
    
  } else { // Execute geocode distance calculation
    
    double xsize = lat1d_vec.length();
    double xsize2 = lon1d_vec.length();
    double ysize = lat2d_vec.length();
    double ysize2 = lon2d_vec.length();
    
    //Testing length equality
    Rcpp::StringVector size_varies1 = "Error: lat1d_vec and lon1d_vec are not same length";
    Rcpp::StringVector size_varies2 = "Error: lat2d_vec and lon2d_vec are not same length";
    
    //Report error if x1,x2 or y1,y2 are of different lengths
    if ( xsize != xsize2){
      
      Rcout << size_varies1 <<"\n";
      return 0;
      
    } else if ( ysize != ysize2){
      
      Rcout << size_varies2 <<"\n";
      return 0;
      
    } else { //Distance calculation
      
      NumericVector out(xsize); // final holder for return
      NumericVector hold(ysize); // holder which we then take minimum of at end of each j loop
    
      for (int i = 0; i < xsize; i++) { // Initiate voter loop
    
        for (int j = 0; j < ysize; j++) { // initiate box loop
    
          hold[j] = distanceEarth(lat1d_vec[i], lon1d_vec[i], lat2d_vec[j], lon2d_vec[j]); // Haversine Calc
    
        }
    
        out[i] = hold[vecminInd(hold)]; // select closest box
  
      }
      return out;
    }  
    
  }

}


// [[Rcpp::export]]
List nearest_dbox2(NumericMatrix lat_lon1d, NumericMatrix lat_lon2d) {
  
  //Turn matrices into vectors
  NumericVector lat1d_vec = lat_lon1d(_,0);
  NumericVector lon1d_vec = lat_lon1d(_,1);
  NumericVector lat2d_vec = lat_lon2d(_,0);
  NumericVector lon2d_vec = lat_lon2d(_,1);
  
  //NA Handling
  Rcpp::StringVector na_error = "Error: NA Exists, Omit NAs from vectors and rerun";
  
  if( any( is_na( lat1d_vec ) ) | any( is_na( lon1d_vec ) ) ) {
    
    Rcout << na_error <<"\n";
    return 0;
    
  } else if (any( is_na( lat2d_vec ) ) | any( is_na( lat2d_vec ) ) ) {
    
    Rcout << na_error <<"\n";
    return 0;
    
  } else { // Execute geocode distance calculation
    
    double xsize = lat1d_vec.length();
    double xsize2 = lon1d_vec.length();
    double ysize = lat2d_vec.length();
    double ysize2 = lon2d_vec.length();
    
    //Testing length equality
    Rcpp::StringVector size_varies1 = "Error: lat1d_vec and lon1d_vec are not same length";
    Rcpp::StringVector size_varies2 = "Error: lat2d_vec and lon2d_vec are not same length";
    
    //Report error if x1,x2 or y1,y2 are of different lengths
    if ( xsize != xsize2){
      
      Rcout << size_varies1 <<"\n";
      return 0;
      
    } else if ( ysize != ysize2){
      
      Rcout << size_varies2 <<"\n";
      return 0;
      
    } else { //Distance calculation
      
      NumericVector out(xsize); // final holder for return
      NumericVector hold(ysize); // holder which we then take minimum of at end of each j loop
      NumericVector row_count(xsize);
      
      for (int i = 0; i < xsize; i++) { // Initiate voter loop
        
        for (int j = 0; j < ysize; j++) { // initiate box loop
          
          hold[j] = distanceEarth(lat1d_vec[i], lon1d_vec[i], lat2d_vec[j], lon2d_vec[j]); // Haversine Calc
          
        }
        
        out[i] = hold[vecminInd(hold)]; // select closest box
        row_count[i] = vecminInd(hold)+1; // Need to add one to port back to R
      }
      return List::create(Rcpp::Named("distance")=out,Rcpp::Named("rows")=row_count); //Output is list() type
    }  
    
  }
  
}


