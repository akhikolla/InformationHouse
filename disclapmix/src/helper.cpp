#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
IntegerMatrix rcpp_create_design_matrix(IntegerMatrix x, int clusters) {
  int n = x.nrow();
  int loci = x.ncol();
  
  int N = n * clusters * loci;
  int a = 0;
      
  IntegerMatrix z(N, 2);
  
  for (int i = 0; i < n; i++) {
    for (int c = 0; c < clusters; c++) {
      for (int l = 0; l < loci; l++) {
       //z(a++, _) = IntegerVector::create(c + 1, l + 1);
        z(a, 0) = c + 1;
        z(a, 1) = l + 1;
        a++;
      }
    }
  }

  return(z);
}

// [[Rcpp::export]]
NumericVector rcpp_create_new_weight_vector(NumericMatrix vic, int loci) {
  int n = vic.nrow();
  int clusters = vic.ncol();
  
  int N = n * clusters * loci;
  int a = 0;
      
  NumericVector z(N);
  
  for (int i = 0; i < n; i++) {
    for (int c = 0; c < clusters; c++) {
      for (int l = 0; l < loci; l++) {
        double v = vic(i, c);
        z(a++) = v;
      }
    }
  }

  return(z);
}

// [[Rcpp::export]]
IntegerVector rcpp_create_response_vector(IntegerMatrix x, IntegerMatrix y) {
  int n = x.nrow();
  int loci = x.ncol();
  int clusters = y.nrow();
  
  if (y.ncol() != loci) {
    throw std::range_error("Different number of loci (columns) in x and y");
  }

  int N = n * clusters * loci;
  int a = 0;

  IntegerVector z(N);
  
  for (int i = 0; i < n; i++) {
    IntegerVector xhap = x(i, _);
    
    for (int c = 0; c < clusters; c++) {
      IntegerVector yhap = y(c, _);
      
      for (int l = 0; l < loci; l++) {
        z(a++) = abs(xhap(l) - yhap(l));
      }
    }
  }

  return(z);
}

// [[Rcpp::export]]
NumericMatrix rcpp_calculate_wic(IntegerMatrix x, IntegerMatrix y, NumericMatrix p, NumericVector tau) {
  int n = x.nrow();
  int loci = x.ncol();
  int clusters = y.nrow();
  
  if (y.ncol() != loci) {
    throw std::range_error("Different number of loci (columns) in x and y");
  }
  
  NumericMatrix wic(n, clusters);
  
  for (int c = 0; c < clusters; c++) {
    IntegerVector yhap = y(c, _);
    
    for (int i = 0; i < n; i++) {
      IntegerVector xhap = x(i, _);
      
      double prod = tau(c);
      
      for (int l = 0; l < loci; l++) {
        double pcl = p(c, l);
        prod *= pow(pcl, abs(xhap(l) - yhap(l)))*((1-pcl)/(1+pcl));
      }
      
      wic(i, c) = prod;
    }
  }

  return(wic);
}

// [[Rcpp::export]]
NumericMatrix rcpp_calculate_vic(NumericMatrix wic) {
  int n = wic.nrow();
  int clusters = wic.ncol();
  
  NumericMatrix vic(n, clusters);
  
  for (int i = 0; i < n; i++) {
    NumericVector wi = wic(i, _);
    double wisum = std::accumulate(wi.begin(), wi.end(), 0.0);
    
    for (int c = 0; c < clusters; c++) {
      vic(i, c) = wi(c) / wisum;
    }
  }
  
  return(vic);
}

// [[Rcpp::export]]
NumericVector rcpp_calculate_haplotype_probabilities(IntegerMatrix new_data, IntegerMatrix y, NumericMatrix p, NumericVector tau) {
  int n = new_data.nrow();
  int loci = new_data.ncol();
  int clusters = y.nrow();

  if (y.ncol() != loci) {
    throw std::range_error("Different number of loci (columns) in x and y");
  }
  
  NumericVector happrobs(n);
  
  for (int i = 0; i < n; i++) {
    IntegerVector h = new_data(i, _);
    double hprob = 0.0;
    
    for (int c = 0; c < clusters; c++) {
      IntegerVector yhap = y(c, _);
      double component_prob = tau(c);
      
      for (int l = 0; l < loci; l++) {
        double pcl = p(c, l);
        component_prob *= pow(pcl, abs(h(l) - yhap(l)))*((1-pcl)/(1+pcl));
      }
      
      hprob += component_prob;
    }
    
    happrobs(i) = hprob;
  }
  
  return(happrobs);
}


// [[Rcpp::export]]
NumericVector rcpp_calculate_haplotype_probabilities_se(IntegerMatrix new_data, IntegerMatrix y, NumericVector theta_clusters, NumericVector theta_loci, NumericMatrix vcov, NumericVector tau) {
  return tau;
}

// [[Rcpp::export]]
NumericMatrix rcpp_calculate_haplotype_probabilities_clusterwise(IntegerMatrix new_data, IntegerMatrix y, NumericMatrix p, NumericVector tau) {
  int n = new_data.nrow();
  int loci = new_data.ncol();
  int clusters = y.nrow();

  if (y.ncol() != loci) {
    throw std::range_error("Different number of loci (columns) in x and y");
  }
  
  NumericMatrix happrobs(n, clusters);
  
  for (int i = 0; i < n; i++) {
    IntegerVector h = new_data(i, _);
    
    for (int c = 0; c < clusters; c++) {
      IntegerVector yhap = y(c, _);
      double component_prob = 1.0;
      
      for (int l = 0; l < loci; l++) {
        double pcl = p(c, l);
        component_prob *= pow(pcl, abs(h(l) - yhap(l)))*((1-pcl)/(1+pcl));
      }
      
      happrobs(i, c) = component_prob;
    }    
  }
  
  return(happrobs);
}


int rdisclap_single(double p) {
  if (p < 0.0 || p >= 1.0) {
    stop("0 <= p < 1 expected");
  }
  
  int res = (int)ceil(log(Rf_runif(0.0, 1.0) * (1 + p)/2)/log(p) - 1);
  return (int)(Rf_runif(0.0, 1.0) < 0.5 ? -1 : 1) * res;
}

// [[Rcpp::export]]
IntegerMatrix rcpp_simulate(int nsim, IntegerMatrix y, NumericVector tau_cumsum, NumericMatrix disclap_parameters) {
  int loci = y.ncol();
  int clusters = y.nrow();
  
  IntegerMatrix res(nsim, loci);
  
  for (int i = 0; i < nsim; i++) {
    double clus = Rf_runif(0.0, 1.0);
    int origin = clusters - 1;
    
    for (int j = 0; j < clusters; j++) {
      if (clus <= tau_cumsum[j]) {
        origin = j;
        break;
      }
    }
    
    IntegerVector hap = y(origin, Rcpp::_);
    
    for (int k = 0; k < loci; k++) {
      hap[k] += rdisclap_single(disclap_parameters(origin, k));
    }
    
    res(i, Rcpp::_) = hap;
  }
  
  return res;
}

// [[Rcpp::export]]
int rcpp_find_haplotype_in_matrix(const IntegerMatrix subpop, const IntegerVector h) {
  int nrows = subpop.nrow();
  int ncolumns = subpop.ncol();
  for (int i = 0; i < nrows; i++) {
    bool found = true;
    
    for (int j = 0; j < ncolumns; j++) {
      if (subpop(i, j) != h(j)) {
        found = false;
        break;
      }
    }

    if (found) {
      return i+1; // For R 1-based indexing
    }
  }
  return -1;
}

