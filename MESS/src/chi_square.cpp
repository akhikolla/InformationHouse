// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
// #include <RcppArmadilloExtensions/sample.h>
using namespace Rcpp;


// Compute chisq test statistic with no correction
double chisq_test(arma::mat x) {
  arma::colvec rowsum = sum(x, 1);
  arma::rowvec colsum = sum(x, 0);

  arma::mat expected = kron(colsum, rowsum)/sum(rowsum);

  double res = accu(square((x - expected))/expected);

  return((res<0.0000000000001) ? 0 : res);
}


// Helper function to convert between arma and Rcpp

template <typename T>
Rcpp::NumericVector arma2vec(const T& x) {
    return Rcpp::NumericVector(x.begin(), x.end());
}



// [[Rcpp::export(.chisq_test_cpp)]]
List chisq_test_cpp(NumericMatrix x, int margin=0, int statistic=1, int B=100000) {

  // Margin: 0 - N, 1 - rows, 2 - both

  arma::mat X(x.begin(), x.nrow(), x.ncol(), false);
  double N = accu(X);
  int nrows = x.nrow();
  int ncols = x.cols();

  // Sanity checks
  if (B<1)
    Rcpp::stop("The number of permutations must be greater than 0");

  // Check that no missing
  // Check that min >= 0
  // Check that integers

  // Compute the test statistic for the original data
  double originaltt = chisq_test(X);

  // Compute the expected values
  arma::colvec rowsum = sum(X, 1);
  arma::rowvec colsum = sum(X, 0);
  arma::mat expected = kron(colsum, rowsum)/sum(rowsum);

  int larger=0;

  if (margin==0) {  // Fixed N
    // Compute sample size and cell percentages
    arma::colvec expvec(expected.memptr(), nrows*ncols, false);
    expvec = normalise(expvec, 1);

    NumericVector cellprobs;
    cellprobs = arma2vec(expvec);

    //    Rcout << "kjhkjH" <<  cellprobs << std::endl;
    IntegerVector ans(cellprobs.size());
    NumericVector ans2(nrows*ncols);

    //    IntegerVector E;
    //    E = arma2vec(expected);

    arma::mat simtable(nrows, ncols);

    for (int i=0; i<B; i++) {
      rmultinom(N, cellprobs.begin(), cellprobs.size(), ans.begin());
      ans2 = as<NumericVector>(ans);  // Convert to numeric.
      // Could maybe get a small speedup here by not having to do this
      // and work directly on the ans vector

      //      sum(Rcpp::pow((ans2-E),2)/E);
      
      arma::mat simtable(ans2.begin(), nrows, ncols, false);
      
      if (chisq_test(simtable) >= originaltt)
	larger +=1;
    }
    
    //    Rcout << "Answer : " << ans << " og " <<  (larger+1.0)/(B+1.0) << std::endl;
    
  } else if (margin==1) {
    // Compute row sizes as they will be fixed
    arma::colvec rowsize = sum(X, 1);

    // Compute column probabilities under H0
    arma::rowvec colprob = sum(X, 0)/N;

    NumericVector colpercent = arma2vec(colprob);
    IntegerVector frame = seq_len(ncols); // The categories to choose from
    
    arma::mat simtable(nrows, ncols);

    for (int i=0; i<B; i++) {
      // Empty the simulated table
      simtable.zeros();    

      // Simulate each row    
      for (int j=0; j<nrows; j++) {
	IntegerVector res = sample(frame, rowsize(j), TRUE, colpercent) ;
	// Fill up the simulated table. Should be possible to do faster
	for (int k = 0; k < rowsize(j); k++) {
	  simtable(j, res(k)-1) +=1;
	}
      }
      if (chisq_test(simtable) >= originaltt)
	larger +=1;
    }
      
    
  } else { // margin = 2

    // Check various conditions

    for (int i=0; i<B; i++) {
      // Empty the simulated table
      //      simtable.zeros();    


      //      R::r2dtable(1, rowsum, colsum);

      
      /*

      
      // Simulate each row    
      for (int j=0; j<nrows; j++) {
	IntegerVector res = sample(frame, rowsize(j), TRUE, colpercent) ;
	// Fill up the simulated table. Should be possible to do faster
	for (int k = 0; k < rowsize(j); k++) {
	  simtable(j, res(k)-1) +=1;
	}
      }

      
      if (chisq_test(simtable) >= originaltt)
	larger +=1;
      */

    }
    

    
  }

  
    /* Slightly faster but more convoluted

    int row=0;
    int idx = 0;
    for (int k = 0; k < N; k++) {
      if (idx == n(row)) {
	row +=1;
	idx=0;
      }
      testtable(row, res(k)-1) +=1;
      idx +=1;
    }

     */


  //  }

  NumericVector stat = NumericVector::create(_["X-squared"] = originaltt) ;

  Rcpp::List RVAL =  Rcpp::List::create(Rcpp::Named("method") = "Two-sided contingency table test with fixed",
					Rcpp::Named("statistic") = stat,
					Rcpp::Named("p.value") = (larger+1.0)/(B+1.0));

  RVAL.attr("class") = "htest";
  
  return(RVAL);
}

