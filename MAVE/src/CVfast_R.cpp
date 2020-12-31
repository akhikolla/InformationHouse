#include <vector>
#include <map>
#include <algorithm>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace std;

arma::mat tiedrank(arma::mat x){
  int nrow=x.n_rows,ncol=x.n_cols;
  arma::mat x_sorted=sort(x);
  arma::mat t_rank;t_rank.zeros(nrow,ncol);
  arma::mat rank;
  rank.zeros(x.n_rows,x.n_cols);
  for(int i=0;i<ncol;++i){

    map<double, double> val2t_rank;
    double flag=x_sorted(0,i);
    int ppos=0;
    for(int j=1;j<nrow;++j){
      if((abs(x_sorted(j,i)-flag)>1e-12)){
        for(int k=ppos;k<j;++k){
          rank(k,i)=(ppos+j-1.0)/2;
          val2t_rank.insert(make_pair(x_sorted(k,i),(ppos+j-1.0)/2));
        }
        flag=x_sorted(j,i);
        ppos=j;
      }
    }
    int j = x.n_rows;
    for(int k=ppos;k<j;++k){
      rank(k,i)=(ppos+j-1.0)/2;
      val2t_rank.insert(make_pair(x_sorted(k,i),(ppos+j-1.0)/2));
    }

    for(int j=0;j<nrow;++j){
      t_rank(j,i)=val2t_rank[x(j,i)];
    }
  }

  return t_rank;
}


bool cmp(pair<double,int>first,pair<double,int>second){
  if(first.first<second.first)
    return false;
  else
    return true;
}


arma::uvec max_num(arma::colvec x,int num) {
  int n=x.size();
  vector<pair<double,int> >win;
  for(int i=0;i<num;++i){
    win.push_back(make_pair(x[i],i));
  }

  make_heap(win.begin(),win.end(),cmp);

  for(int i=num;i<n;++i){

    if(x[i]>win[0].first){
      pop_heap(win.begin(),win.end(),cmp);
      win.pop_back();
      win.push_back(make_pair(x[i],i));
      push_heap(win.begin(),win.end(),cmp);
    }
  }
  arma::uvec maxIdx;
  maxIdx.resize(num);
  for(int i=0;i<num;++i){
    maxIdx[i]=win[i].second;
  }
  return maxIdx;
}


//[[Rcpp::export]]
double CVfastCpp(const arma::mat& x,const arma::mat& ky){
  int n = x.n_rows, p=x.n_cols;
  double cv = 0;
  arma::mat nx = tiedrank(x)/(n+1);
  double h=arma::mean(arma::mean(arma::stddev(nx,0,0)))/(pow(n,1/(p+4.0)));

  for(int i=0;i<n;++i){
    arma::colvec dxi = arma::sum(abs(nx-repmat(nx.row(i),n,1)),1)/h;
    arma::colvec k=1/pow(1+dxi,4);
    k[i]=0;
    k=k/(sum(k)+1e-6);
    arma::rowvec ye=k.t()*ky;
    cv+=arma::mean(abs(ye-ky.row(i)))/n;
  }

  return cv;
}
