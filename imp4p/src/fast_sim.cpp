
#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

inline NumericVector fast_si( NumericVector& prot, NumericMatrix& mat) {

  int nr = mat.nrow();

  int n = prot.size();
  LogicalVector out_prot(n);
  LogicalVector out_mat(n);

  double m_p=0;
  double m_m=0;
  double k=0;
  //double v_p=0;
  //double v_m=0;
  //double cov_mp=0;
  double d_m=0;
  //double d_m2=0;
  //double aa=0;
  //double bb=0;
  //double cc=0;
  NumericVector prot2(n);
  NumericVector corr(nr);
  for (int j = 0; j < nr; ++j) {
    prot2=mat(j,_);
    m_p=0;
    m_m=0;
    k=0;
    d_m=0;
    //nb_na_prot2=0;
    //nb_na_prot=0;
    //selection des observations pairwise
    for (int i = 0; i < n; ++i) {
      //prot2[i]=mat(j,i);
      out_prot[i] = NumericVector::is_na(prot[i]);
      out_mat[i] = NumericVector::is_na(prot2[i]);
      if ( out_prot[i] == FALSE ){
        if ( out_mat[i] == FALSE ){
          m_p+=prot[i];
          m_m+=prot2[i];
          d_m+=(prot[i]-prot2[i])*(prot[i]-prot2[i]);
          //nombre de valeurs observ?es aux m?mes endroits
          ++k;
        }//else{nb_na_prot2+=1;}
      }//else{nb_na_prot+=1;}
    }
    // aa=0;
    // bb=0;
    // cc=0;
    // for (int i = 0; i < n; ++i) {
    //   aa+=(prot[i]-m_p)*(prot2[i]-m_m);
    //   bb+=(prot[i]-m_p)*(prot[i]-m_p);
    //   cc+=(prot2[i]-m_m)*(prot2[i]-m_m);
    // }
    // d_m2=(aa)/(sqrt(bb)*sqrt(cc));
    //calcul de k+exp(-x^2/2)
    //nb de valeurs manquantes situ?es aux m?mes endroits divis?es par le nb de valeurs manquantes de prot
    //if (ind=1){
    //    d_na = (n-k)/nb_na_prot;
    //    rd_na = NumericVector::rank(d_na);
    //    if ( k>0 ){
    //      //Si prot a au moins une valeur en commun avec prot2
    //      corr[j]=rd_na+exp(-d_m);
    //    }else{corr[j]=NA_REAL;}
    //}else{
        if ( k>0 ){
          //Si prot a au moins une valeur en commun avec prot2
          if (k>1){
              corr[j]=d_m;
          }else{
            corr[j]=d_m;
          }
        }else{corr[j]=NA_REAL;}
    //}
  }

  return corr;
}

// [[Rcpp::export]]
NumericVector fast_sim(NumericVector prot, NumericMatrix mat) {
  return fast_si(prot, mat);
}
