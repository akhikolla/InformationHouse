
#include <Rcpp.h>
#include "typedefs.h"
#include "headers.h"

using namespace Rcpp;

//=============================================================================

// Calculate con and cov of disjunctions
//
//   x   membership scores of antecendens
//   y   membership scores of consequens
// [[Rcpp::export]]
NumericVector C_conCov(const NumericVector x, const NumericVector y, const IntegerVector f){
  int n=x.size();
  NumericVector Sums(3), conCov(2);
  for (int i=0; i<n; i++){
    Sums(0)+=x[i]*f(i);
    Sums(1)+=y[i]*f(i);
    Sums(2)+=std::min(x[i], y[i])*f(i);
  };
  conCov(0) = Sums(2)/Sums(0);
  conCov(1) = Sums(2)/Sums(1);
  return(conCov);
}

//=============================================================================

// C_subsetMin: minimum of a subset of the elements of an integer vector
//   x    vector
//   sub  subset of elements to consider
// [[Rcpp::export]]
double C_subsetMin(const NumericVector x, const IntegerVector sub){
  int len = sub.size();
  double out = x(sub(0)-1);
  if (len == 1) return(out);
  for (int i=1; i<sub.size(); i++){
    out = std::min(x(sub(i)-1), out);
  }
  return(out);
}


// membership scores of a series of conjunctions
//   m   IntegerMatrix, the rows representing the components of the conjunctions
//   x  'scores'-matrix
// [[Rcpp::export]]
NumericMatrix C_conjScore(const NumericMatrix x, const IntegerMatrix m){
  int n=x.nrow(), p=m.nrow();
  NumericMatrix out(n,p);
  for (int i=0; i<n; i++){
    for (int j=0; j<p; j++){
      // Rcpp::Rcout << i << " " << j << " " << C_subsetMin(x(i, _), m(j, _)) << std::endl;
      out(i,j) = C_subsetMin(x(i, _), m(j, _));
    }
  }
  return(out);
}

//==============================================================================

// initialize ii
// [[Rcpp::export]]
IntegerVector C_init_ii(const IntegerVector nn, const LogicalVector st){
  int l=nn.size();
  IntegerVector ini(l);
  ini.fill(0);
  for (int j=0; j<l-1; j++){
    if (st[j]){
      ini[j+1] = ini[j] + 1;
    }
  }
  return(ini);
}

// upper limit of ii
// [[Rcpp::export]]
IntegerVector C_set_lim(const IntegerVector nn, const LogicalVector st){
  int l=nn.size();
  IntegerVector lim(l);
  lim = nn - 1; 
  if (l<=1) return(lim);
  for (int j=l-2; j>=0; j--){
    if (st[j]){
      lim[j] = lim[j+1] - 1;
    }
  }
  return(lim);
}

// C_max_which
// [[Rcpp::export]]
int C_max_which(const LogicalVector x) {
  int n = x.size();
  int i;
  for (i=n-1; i >= 0; i--) {
    if (x[i]) break;
  }
  return(i+1);
}


// increment ii
// [[Rcpp::export]]
IntegerVector C_increment(IntegerVector ii, const IntegerVector nn, 
                          const LogicalVector st, const IntegerVector lim){
  int l = ii.size();
  if (ii[l-1] < lim[l-1]){    // end of loop not yet reached in last position
    ii[l-1] += 1;
    return(ii);
  } 
  // p = first position with increment in index
  int p=C_max_which(ii < lim);
  if (p == 0){
    ii.fill(0);
    return(ii);
  }
  // increment index positions < p
  if (p==1){
    ii[0] += 1; 
  } else {
    ii[Range(0, p-1)] = C_increment(ii[Range(0, p-1)], nn[Range(0, p-1)], st[Range(0, p-2)], lim[Range(0, p-1)]);
  }
  // reset index positions >= p
  for (int j=p; j<l; j++){
    if (st[j-1]){
      ii[j] = ii[j-1] + 1;
    } else {
      ii[j] = 0;
    } 
  }
  return(ii);
}


/*
// [[Rcpp::export]]
int test_increment(IntegerVector nn, LogicalVector st){
  IntegerVector ii=C_init_ii(nn, st);
  IntegerVector lim=C_set_lim(nn, st);
  int count=0;
  do {
    C_increment(ii, nn, st, lim);
    count++;
  } while (as<bool>(any(ii>0)));
  return(count);
} */


/* R
(nn <- rep(2:4, c(2, 1, 3)))
st <- diff(nn) == 0
C_set_lim(nn, st)
C_init_ii(nn, st)

test_increment(c(4L, 4L), T)
test_increment(c(4L, 4L), F)
test_increment(c(3L, 5L, 7L), c(F, F))

get_n <- function(nn){
  choose_args <- unname(unclass(rle(nn)[2:1]))
  prod(do.call(mapply, c(list(choose), choose_args)))
}
mytest <- function(nn){ 
  n <- get_n(nn)
  st <- diff(nn) == 0L
  cat(n, "\n")
  stopifnot(test_increment(nn, st) == n)
}
mytest(c(3L, 5L, 7L))
mytest(c(1L, 2L, 3L, 3L, 4L, 5L, 5L, 5L, 6L, 7L))
mytest(sort(c(1L, 2L, 3L, 3L, 4L, 5L, 5L, 5L, 6L, 7L)))
mytest(sample(c(1L, 2L, 3L, 3L, 4L, 5L, 5L, 5L, 6L, 7L)))

*/

//==============================================================================


// Find asf's of a given structure (=lengths of the conjunctions in the disjunctions)
// [[Rcpp::export]]
IntegerMatrix C_find_asf(const IntegerVector conjlen, const numMatList x, 
                         const NumericVector y, const IntegerVector f,
                         const double con, const double cov, 
                         const int maxSol){
  int n_conj = conjlen.size();

  // st = indicator, =true if conj hast same length as precedent conj
  LogicalVector st(n_conj - 1);
  st = (Rcpp::diff(conjlen) == 0);
  
  // nn = number of conjunctions
  IntegerVector nn(n_conj);
  for (int i=0; i<n_conj; i++){
    nn(i)=as<NumericMatrix>(x[i]).ncol();
  }

  // intialize ii, count, out; define lim
  IntegerVector ii=C_init_ii(nn, st);
  int count=0;
  IntegerMatrix out(maxSol, n_conj);
  IntegerVector lim=C_set_lim(nn, st);
  
  // Main loop over combinations
  do {
    // Calculate membership Scores for disjunction
    //Rcpp::Rcout << ii << " count=" << count << std::endl;
    NumericVector ms=as<NumericMatrix>(x[0])(_, ii[0]);
    //Rcpp::Rcout << /*"min: " << minim << ", */ "max: " << ms << std::endl;
    for (int i=1; i<n_conj; i++){
      ms=pmax(ms, as<NumericMatrix>(x[i])(_, ii[i]));
      //Rcpp::Rcout << /*"min: " << minim << ", */ "max: " << ms << std::endl;
    }
    // calculate con und cov
    NumericVector coco=C_conCov(ms, y, f);
    //Rcpp::Rcout << "coco: " << coco << std::endl;
    
    // Conditionally insert in next row of the matrix
    if (coco(0) >= con && coco(1) >= cov){
      for (int j=0; j<n_conj; j++){
        out(count, j) = ii[j];
      }
      count++;
    }

    // increment:
    C_increment(ii, nn, st, lim);
    
    if (count>=maxSol){
      Rcpp::Rcout << "Not all candidate solutions are recorded (maxSol="<< maxSol << ")" 
                  << std::endl;
      break;
    }
  } while (as<bool>(any(ii>0)));
  
  // if (count > 0){
  //   Rcpp::Rcout << "count="<< count << std::endl;
  // }
  
  // Return matrix with 0 rows if count=0:
  if (count == 0){ 
    IntegerMatrix a(0, n_conj);
    return a;
  }
  
  // Return result matrix
  return out(Range(0, count-1), _);
}

