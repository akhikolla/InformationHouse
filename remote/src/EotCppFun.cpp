#include <Rcpp.h>
using namespace Rcpp;


///////////////// linear model functions //////////////////////////////////
// [[Rcpp::export]]
double corC(NumericVector x, NumericVector y) {
  int nx = x.size(), ny = y.size();
  
  if (nx != ny) stop("Input vectors must have equal length!");
  
  double sum_x = sum(x), sum_y = sum(y);
  
  NumericVector xy = x * y;
  NumericVector x_squ = x * x, y_squ = y * y;
  
  double sum_xy = sum(xy);
  double sum_x_squ = sum(x_squ), sum_y_squ = sum(y_squ);
  
  double out = ((nx * sum_xy) - (sum_x * sum_y)) / sqrt((nx * sum_x_squ - pow(sum_x, 2.0)) * (nx * sum_y_squ - pow(sum_y, 2.0)));
  
  return out;
}


// [[Rcpp::export]]
List lmC(NumericVector x, NumericVector y) {
  // Check if input vectors have equal length
  int nx = x.size(); // C++ size() is equal to R length()
  int ny = y.size();
  // Throw error message in case of unequal length
  if (nx != ny) stop("Input vectors must have equal length!");

  // Calculate r value based on linear regression formula
  double sum_x = sum(x); // C++ sum() is equal to R sum()
  double sum_y = sum(y);
  
  NumericVector xy = x * y; // Arithmetic operations work the same way as in R
  NumericVector x_squ = x * x;
  NumericVector y_squ = y * y;
  
  double sum_xy = sum(xy);
  double sum_x_squ = sum(x_squ);
  double sum_y_squ = sum(y_squ);
  
  // C++ pow() is equal to R ^
  double lm_r = ((nx * sum_xy) - (sum_x * sum_y)) / sqrt((nx * sum_x_squ - pow(sum_x, 2.0)) * (nx * sum_y_squ - pow(sum_y, 2.0)));

  // Calculate intercept and slope
  double lm_intercept = ((sum_y * sum_x_squ) - (sum_x * sum_xy)) / ((nx * sum_x_squ) - pow(sum_x, 2.0));
  double lm_slope = ((nx * sum_xy) - (sum_x * sum_y)) / ((nx * sum_x_squ) - pow(sum_x, 2.0));

  // Calculate residuals
  NumericVector lm_predicted_values = lm_slope * x + lm_intercept;
  NumericVector lm_residuals = y - lm_predicted_values;
  
  // Calculate test statistics
  int lm_df = ny - 2; // Degrees of freedom 
  
  double lm_stderr = sqrt(sum(pow(lm_residuals, 2.0)) / lm_df) / 
                     sqrt(sum(pow(x - mean(x), 2.0))); // Standard error
                     
  double lm_tscore = lm_slope / lm_stderr; // t score

  
  // List and return single parameters
  return List::create(lm_r, lm_intercept, lm_slope, lm_residuals, lm_tscore, lm_df);
}


// [[Rcpp::export]]
NumericVector predRsquaredSum(NumericMatrix pred_vals, NumericMatrix resp_vals, 
                              bool standardised) {
  // Number of rows of input matrices
  int nrow_pred = pred_vals.nrow(), nrow_resp = resp_vals.nrow();
  
  // Loop through all predictor cells
  NumericVector lm_rsq_sum(nrow_pred);
  for (int i = 0; i < nrow_pred; i++) {
    
    // For the current predictor cell, loop through all response cells
    // and calculate the corresponding R-squared value
    NumericVector lm_rsq(nrow_resp);
    for (int j = 0; j < nrow_resp; j++) {
      lm_rsq[j] = pow(corC(pred_vals(i, _), resp_vals(j, _)), 2.0);
      
      // Perform standardisation (optional)
      if (!standardised) {
        lm_rsq[j] = lm_rsq[j] * var(resp_vals(j, _));
      }
      
      // Replace possible NaN with 0
      if (lm_rsq[j] != lm_rsq[j]) {
        lm_rsq[j] = 0;
      }
    }
    
    // Sum up R-squared values of current predictor cell
    lm_rsq_sum[i] = sum(lm_rsq);
  }

  return lm_rsq_sum;
}


// [[Rcpp::export]]
List respLmParam(NumericMatrix x, NumericMatrix y, int cell) {
  int nrow_y = y.nrow();
  
  List lmC_out(nrow_y);
  for (int i = 0; i < nrow_y; i++) {
    lmC_out[i] = lmC(x(cell, _), y(i, _));
  }
  
  return lmC_out;
}

///////////////// index of agreement functions ////////////////////////////
// [[Rcpp::export]]
NumericVector findudC(NumericVector x) {
  NumericVector v = diff(x);
  NumericVector z = v.size();
  NumericVector mm(v.size(), -1.0);
  NumericVector pp(v.size(), 1.0);
  NumericVector res = ifelse( v > z, pp, mm);
  return res;
}

// [[Rcpp::export]]
double iodaC(NumericVector x, NumericVector y) {
  NumericVector hh(x.size(), 1.0);
  NumericVector mm(x.size(), 0.0);
  NumericVector e = ifelse(findudC(x) == findudC(y), hh, mm);
  double m = mean(e);
  return m;
}

// [[Rcpp::export]]
NumericVector iodaSumC(NumericMatrix pred_vals, NumericMatrix resp_vals) {
  // Number of rows of input matrices
  int nrow_pred = pred_vals.nrow(), nrow_resp = resp_vals.nrow();
  
  // Loop through all predictor cells
  NumericVector ioda_sum(nrow_pred);
  for (int i = 0; i < nrow_pred; i++) {
    
    NumericVector ioda(nrow_resp);
    for (int j = 0; j < nrow_resp; j++) {
      
      ioda[j] = iodaC(pred_vals(i, _), resp_vals(j, _));
      
    }
    
    // Sum up IOA values of current predictor cell
    ioda_sum[i] = sum(ioda);
  }

  return ioda_sum;
}

///////////////// deseasoning functions ////////////////////////////
// [[Rcpp::export]]
NumericMatrix monthlyMeansC(NumericMatrix x, int nCycleWindow) {
  
  // input matrix: number of rows and columns
  int nRows = x.nrow(); 
  int nCols = x.ncol();  
  
  // initialize vector for temporary storage of values per time step
  int nVecLen = nCols / nCycleWindow;
  NumericVector adValues(nVecLen);
  
  // initialize counter
  int n;
  
  // initialize output matrix
  NumericMatrix mdMeans(nRows, nCycleWindow);
  
  // loop over rows (pixels)
  for (int i = 0; i < nRows; i++) {
    
    // per row, loop over cycle window (time steps, e.g. 0-11 for 12 months) 
    for (int j = 0; j < nCycleWindow; j++) {
      
      // reset counter
      n = 0;
      
      // for the current time step, extract values from input matrix
      // (e.g. for january, extract the 1st, 13th, 24th, ... value)
      for (int k = j; k < nCols; k = k + nCycleWindow) {
        
        adValues[n] = x(i, k);
        
        // increment counter by 1        
        n += 1;
        
      }
      
      // for the current time step, insert mean value into output matrix
      mdMeans(i, j) = mean(adValues);
      
    }
  }
  
  return mdMeans;
}

///////////////// denoise functions ////////////////////////////
// [[Rcpp::export]]
NumericMatrix insertReconsC(List lRecons, NumericMatrix mdTemplate) {
  
  int nListLength = lRecons.size();
  // std::cout << nListLength;
  
  NumericVector dListSlot;
  for (int i = 0; i < nListLength; i++) {
    // current slot entries
    dListSlot = lRecons[i];
    // std::cout << dListSlot[30] << "\n";
    
    // insert values into matrix
    mdTemplate(_, i) = dListSlot;
  }
  
  return mdTemplate;
}
