#include <Rcpp.h>
using namespace Rcpp;

//' rcpp_fastmatclcr
//' 
//' @param I Matrix
//' @param w vector
//' @param MSEmat Matrix
//' @param S Vector
//' @param maxlevel Integer
//' @return Nothing, void
//' @export
// [[Rcpp::export]]
void rcpp_fastmatclcr(NumericMatrix I, NumericVector w, NumericMatrix MSEmat, NumericVector S, int maxlevel){
  
  int ns = S.size();
  int nb = w.size();
  int d = I.ncol();
  
  NumericVector RR(ns);
  
  for(int bindex=0; bindex<=(nb-1); bindex++){
    if(fabs(w[bindex])>0.5){
    RR = MSEmat( _ ,I(bindex,0)-1);
    for(int dim=1;dim<=(d-1); dim++){
      RR = RR*MSEmat( _ ,maxlevel*dim+I(bindex,dim)-1);
    }
    RR = -w[bindex]*RR;
     
     S += RR;
    }
  }
  return;
}


#include <Rcpp.h>
using namespace Rcpp;

//' rcpp_fastmatclcranddclcr
//' 
//' @param I Matrix
//' @param w vector
//' @param MSEmat Matrix
//' @param dMSEmat Matrix
//' @param S Vector
//' @param dS Matrix
//' @param maxlevel Integer
//' @param numpara Integer
//' @return Nothing, void
//' @export
// [[Rcpp::export]]
void rcpp_fastmatclcranddclcr(NumericMatrix I, NumericVector w, NumericMatrix MSEmat,
                              NumericMatrix dMSEmat, NumericVector S, NumericMatrix dS, int maxlevel, int numpara){
  
  int ns = S.size();
  int nb = w.size();
  int d = I.ncol();
  
  NumericVector RR(ns);
  
  for(int bindex=0; bindex<=(nb-1); bindex++){
    if(fabs(w[bindex])>0.5){
    RR = MSEmat( _ ,I(bindex,0)-1);
    for(int dim=1;dim<=(d-1); dim++){
      RR = RR*MSEmat( _ ,maxlevel*dim+I(bindex,dim)-1);
    }
    RR = -w[bindex]*RR;
    
    S += RR;
    
    for(int dim=0;dim<=(d-1); dim++){
      for(int para=0;para<=(numpara-1); para++){
        dS( _ , numpara*dim+para) = dS( _ , numpara*dim+para)+RR*dMSEmat( _ ,maxlevel*numpara*dim+para*maxlevel+I(bindex,dim)-1);
      }
    }
    }
  }
  return;
}
