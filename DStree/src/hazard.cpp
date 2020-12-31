#include <Rcpp.h>


using namespace Rcpp;

// [[Rcpp::export]] 

List hazard(NumericVector y1, NumericVector y2,NumericVector lev,NumericVector wt){
  
  
  int n = y1.size();
  int k = 0;
  int z = lev.size(); 
  Rcpp::NumericVector ncens(z);  
  Rcpp::NumericVector uncens(z);
  Rcpp::NumericVector ncensk(z-1);
  Rcpp::NumericVector nuncensk(z-1);
  Rcpp::NumericVector r(z);
  Rcpp::NumericVector h(z);
  Rcpp::NumericVector r1(z-1);
  Rcpp::NumericVector cum(z-1);
  Rcpp::NumericVector S(z);  
  Rcpp::NumericVector S1(z+1);
  Rcpp::NumericVector pi(z);  


  
  for(int i=0;i<n;i++){
    if(y2(i)==1){
      k=k+1;
      }
  }
  
  int l = n-k;
  Rcpp::NumericMatrix datcens(k,3);
  Rcpp::NumericMatrix datuncens(l,3);
  
  int g = 0; 
  for(int u=0;u<n;u++){
  if(y2(u)==1){
      datcens(g,0) = y1(u);
      datcens(g,1) = y2(u);
      datcens(g,2) = wt(u);
      g = g+1;

    }
  }
  
  int m = 0;
  for(int a=0;a<n;a++){
    if(y2(a)==2){
      datuncens(m,0) = y1(a);
      datuncens(m,1) = y2(a);
      datuncens(m,2) = wt(a);
      m = m+1;
      }
  }
  
  for(int i=0;i<z;i++){
   for(int j=0;j<k;j++){
    if(datcens(j,0)==lev(i)){
      ncens(i) = ncens(i)+1;
    }
   }
    for(int p=0;p<l;p++){
     if(datuncens(p,0)==lev(i)){
      uncens(i) = uncens(i)+1;
     }
    }
  }
  
  for(int i=0;i<z-1;i++){
  ncensk(i) =  ncens(i);
  nuncensk(i) = uncens(i);
  }
  
  r(0) = sum(wt);
  //cum=cumsum(ncensk+nuncensk);
  
  double acc = 0;
        
  for(int i = 0; i < z-1; i++){
         acc += ncensk(i)+nuncensk(i);
         cum(i) = acc;
    }
  
  for(int i=0;i<z-1;i++){
  r1(i) = r(0)-cum(i);
  }
  
  for(int i=1;i<z;i++){
  r(i) =  r1(i-1);
  }
  
  for(int i=0;i<z;i++){
    if(r(i)==0){
      h(i)=0;
    }
    else{
  h(i)=uncens(i)/r(i);
   }
  }
  
  
  S1(0)=1;
  for (int i=0;i<z;i++){
    S1(i+1) = S1(i)*(1-h(i));
  }
  for (int i=0;i<z;i++){
    S(i) = S1(i+1);
  }
  
  pi(0)=1-S(0);
  for (int i=1;i<z;i++){
    pi(i) = S(i-1)-S(i);
  }
  
  return Rcpp::List::create(Rcpp::Named("n.cens") = ncens,
                          Rcpp::Named("n.uncens") = uncens,
                          Rcpp::Named("pi") = pi,
                          Rcpp::Named("h") = h,
                          Rcpp::Named("S") = S
                          );
  
}

