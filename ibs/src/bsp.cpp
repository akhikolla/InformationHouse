#include<Rcpp.h>
using namespace Rcpp;
double gdiv(double a,double b){if(a==0.0&&b==0.0)return(0.0);else return(a/b);}
double bsp(int i, int ord, double x, int nk, NumericVector kns){
  if(i<0||i>nk-ord-1){
    Rcout<<"illegal i value: i="<<i<<"; nk-ord="<<nk<<"-"<<ord<<"="<<nk-ord<<std::endl;
    return R_NaN;
  }
  if(x<kns[i]||x>kns[i+ord])return(0.0);
  int k=nk-1; while(kns[k]==kns[k-1])k--; k--; 
  if(ord==1){
    if(i!=k){
      return((kns[i]<=x && x<kns[i+1])? 1.0 : 0.0);
    }else{
      if(i==k){
	return((kns[i]<=x && x<=kns[i+1])? 1.0 : 0.0);
      }else return R_NaN;
    }
  }else
    return(gdiv((x-kns[i])*bsp(i,ord-1,x,nk,kns), kns[i+ord-1]-kns[i])+
	   gdiv((kns[i+ord]-x)*bsp(i+1,ord-1,x,nk,kns),kns[i+ord]-kns[i+1])
	   );
}
//[[Rcpp::export]]
NumericVector bsbasesCpp(NumericVector xs, NumericVector kns, int order){
  int nx=xs.size(), nk=kns.size(),i,j;
  int ord=order, nb=nk-ord;
  int ansSize=nb*nx;
  NumericVector ans(ansSize);
  for(i=0;i<nx;i++){
    for(j=0;j<nb;j++){
      ans[j+i*nb]=bsp(j,ord,xs[i],nk,kns);
    }
  }
  return(ans);
}
//[[Rcpp::export]]
NumericVector bsplineCpp(NumericVector xs, int ord, NumericVector kns, NumericVector coef){
  int nx=xs.size(), nk=kns.size(), i, j;
  NumericVector ans(nx,0.0);
  for(i=0;i<nx;i++)
    for(j=0; j<nk-ord; j++)
      ans[i] += coef[j]*bsp(j,ord,xs[i],nk,kns);
  return(ans);
}
//[[Rcpp::export]]
NumericVector ibsCpp(NumericVector xs, int ord, NumericVector kns, NumericVector coef){
  int nx=xs.size(), nk=kns.size(), nc=coef.size(),i, j, k;
  double tmp;
  NumericVector ans(nx,0.0),kns_ext(nk+1),bet(nc);
  for(i=0;i<nk;i++){
    kns_ext[i]=kns[i];
  }
  kns_ext[nk]=kns_ext[nk-1];
  bet[0]=coef[0]*(kns_ext[ord]-kns[0])/ord;
  for(i=1;i<nc;i++){
    bet[i]=bet[i-1]+coef[i]*(kns_ext[i+ord]-kns_ext[i])/ord;
  }
  ans=bsplineCpp(xs,ord+1,kns_ext,bet);
  return(ans);
}
