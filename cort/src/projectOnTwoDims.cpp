#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
Rcpp::List projectOnTwoDims(const NumericMatrix a,
                           const NumericMatrix b,
                           const NumericVector p,
                           const NumericVector f,
                           const NumericVector kern,
                           const NumericVector dims,
                           const NumericMatrix data) {

  // Intialise objects :
  int n_leaves = a.ncol();
  int dim = a.nrow();
  int n_obs = data.nrow();
  double prod;
  int n_edges1,n_edges2,new_n,k,l;
  Rcpp::NumericVector edges1(2*n_leaves);
  Rcpp::NumericVector edges2(2*n_leaves);

  // Compute edges :
  edges1[seq(0,n_leaves-1)] = a(dims(0)-1,_);
  edges1[seq(n_leaves,2*n_leaves-1)] = b(dims(0)-1,_);
  edges1 = sort_unique(edges1);

  edges2[seq(0,n_leaves-1)] = a(dims(1)-1,_);
  edges2[seq(n_leaves,2*n_leaves-1)] = b(dims(1)-1,_);
  edges2 = sort_unique(edges2);

  n_edges1 = edges1.length() - 1;
  n_edges2 = edges2.length() - 1;
  new_n = n_edges1*n_edges2;

  // Setup variables :
  Rcpp::NumericMatrix new_a(new_n,2);
  Rcpp::NumericMatrix new_b(new_n,2);
  Rcpp::NumericVector new_f(new_n);
  Rcpp::NumericVector new_p(new_n);
  Rcpp::NumericVector new_vols(new_n);
  new_p.fill(0.0);
  new_f.fill(0.0);

  for(int i =0; i<new_n;i++){
    // Construct new_a and new_b :
    k = i / n_edges2;
    l = i % n_edges2;
    new_a(i,0) = edges1[k];
    new_a(i,1) = edges2[l];
    new_b(i,0) = edges1[k+1];
    new_b(i,1) = edges2[l+1];

    // Construct vols :
    new_vols(i) = (new_b(i,0)-new_a(i,0))*(new_b(i,1)-new_a(i,1));

    // Construct f :
    for(int n=0; n < n_obs; n++){
      if((new_a(i,0) <= data(n,0))&&(data(n,0) < new_b(i,0))&&(new_a(i,1) <= data(n,1))&&(data(n,1) < new_b(i,1))){
        new_f(i) += 1.0/n_obs;
      }
    }

    // Construct p :
    for(int m=0; m< n_leaves; m++){
      prod = kern(m);
      if(prod != 0){
        for(int d = 0; d < dim; d++){
          if(d  == (dims(0)-1)){
            prod *= std::max(std::min(b(d,m),new_b(i,0))-std::max(a(d,m),new_a(i,0)),0.0);
          } else if(d == (dims(1)-1)){
            prod *= std::max(std::min(b(d,m),new_b(i,1))-std::max(a(d,m),new_a(i,1)),0.0);
          } else {
            prod *= b(d,m)-a(d,m);
          }
          if(prod == 0){
            break;
          }
        }
        new_p(i) += prod;
      }
    }
  }
  Rcpp::List L = List::create(Named("f") = new_f , _["p"] = new_p, _["a"] = new_a, _["b"] = new_b, _["vols"] = new_vols);
  return (L);
}




















