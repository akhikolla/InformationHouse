#include <Rcpp.h>
using namespace Rcpp;

//' Creates C++ loop 
//'
//' Meant to be a helper within date_check(), and subsequently propagate_date()
//
//' @param x Our date vector
//' @param y Our days supply vector
//' @return A new vector to be appended to a users dataframe with adjusted dates

// [[Rcpp::export(name = ".date_checkCpp")]]
NumericVector date_checkCpp(NumericVector x, NumericVector y){
  int data_length = x.size();
  NumericVector out(data_length);
  
  out[0] = x[0];  //Set our out vectors first date to the input vectors first date (will never be adjusted) 
  for(int i = 1; i < data_length; ++i){
    int prior_date = out[i-1] + y[i-1]; //this is basically like saying x[0] + y[0] to begin, but then keeps track of the adjusted date after the first loop iteration by using out[]
    int date_to_check = x[i]; // We need to check the input buffer next date here, not out[i]. out[i] is uninitlazied (to 0 or a random value) until we actually set it.
    
    // if the date prior is greater than the current date... 
    if(prior_date > date_to_check){
      // update the current date if it overlaps with what the next date currently is
      out[i] =  prior_date;
    } else {
      out[i] = x[i];  // supply lasts long enough, just keep the original next date of resupply
    }
  }
  
  return out;
}

 
