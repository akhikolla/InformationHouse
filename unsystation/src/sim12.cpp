#include <RcppArmadilloExtensions/sample.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

//' @keywords internal
// [[Rcpp::export]]
NumericMatrix funcSimX(NumericVector coef, NumericMatrix buffMat){

  int i, j, t;
	int p = coef.size();
	int N = buffMat.ncol();
	int T = buffMat.nrow()-p;
	NumericMatrix simX(T, N);

	for(i=0; i<N; i++){
		for(t=0; t<T; t++){
			simX(t, i) = buffMat(t+p, i);
			for(j=0; j<p; j++){
				if(t-j-1<0){
					simX(t, i) += coef(j)*buffMat(t+p-j-1, i);
				} else{
					simX(t, i) += coef(j)*simX(t-j-1, i);
				}
			}
		}
	}

	return(wrap(simX));

}

//' @keywords internal
// [[Rcpp::export]]
NumericVector func_coef(NumericVector z, int scale){
	int len = z.size();
	int lenw = pow(2.0, -scale);
	int t, j;
	NumericVector coef(len);
	NumericVector wave(lenw);
	NumericVector pad(2*len);

	for(j=0; j<lenw/2; j++){
		wave(j) = sqrt(pow(2.0, scale));
		wave(j+lenw/2) = -wave(j);
	}

	for(t=0; t<len; t++) pad(t) = pad(len+t) = z(t);

	for(t=0; t<len; t++){
		for(j=0; j<lenw; j++) coef(t) += pad(t+j)*wave(j);
	}

 	return(coef);
}

//' @keywords internal
// [[Rcpp::export]]
NumericMatrix funcSEMat(int T, int minLength, NumericVector probVec, IntegerVector topCand, IntegerVector bottomCand){

  	int i, j, br, s, e, t, iter;
  	int M = topCand.size();
  	NumericMatrix seMat(6, M);
	IntegerVector se(2);
	IntegerVector frame = Rcpp::Range(1, T);

	for(i=0; i<M; i++){
		iter = 0;
		t = topCand(i);
		while(TRUE){
			se = RcppArmadillo::sample(frame, se.size(), FALSE, probVec);
			s = min(se); e = max(se);
			if(s>t || e<=t) continue;
			if(e-s+1 <= minLength) continue;
			iter += 1;
			if(i==0) break;
			br = 0;
			for(j=0; j<i; j++){
				if(s==seMat(0, j) && e==seMat(2, j)){
					br = 1;
				}
			}
			if(iter==10 || br==0) break;
		}
		seMat(0, i) = s; seMat(1, i) = t; seMat(2, i) = e;

		t = bottomCand(i);
		while(TRUE){
			se = RcppArmadillo::sample(frame, se.size(), FALSE, probVec);
			s = min(se); e = max(se);
			if(s>t || e<=t) continue;
			if(e-s+1 <= minLength) continue;
			iter += 1;
			if(i==0) break;
			br = 0;
			for(j=0; j<i; j++){
				if(s==seMat(0, j) && e==seMat(2, j)){
					br = 1;
				}
			}
			if(iter==10 || br==0) break;
		}
		seMat(3+0, i) = s; seMat(3+1, i) = t; seMat(3+2, i) = e;
	}

	return(wrap(seMat));

}

//' @keywords internal
// [[Rcpp::export]]
List funcRes(NumericMatrix yMat, int M, int minLength, NumericVector probVec, IntegerVector topCand, IntegerVector bottomCand, NumericVector var){

  int i, j, k, l, t, s1, t1, e1, s2, t2, e2, ns1, ne1, ns2, ne2, ts1, tt1, te1, ts2, tt2, te2, n1, n2;
  int maxLevel = yMat.ncol();
  int T = yMat.nrow();
  int Msq = M*M;
  double mean1, mean2;
  double tmp = 0;
  NumericMatrix resMat(Msq, 4+maxLevel+1);
  NumericMatrix seMat;

  seMat = funcSEMat(T, minLength, probVec, topCand, bottomCand);

  for(i=0; i<M; i++){
    s1 = seMat(0, i); t1 = seMat(1, i); e1 = seMat(2, i);
    for(j=0; j<M; j++){
      l = M*i+j;

      ts1 = s1; tt1 = t1; te1  = e1;
      ts2 = s2 = seMat(3+0, j); tt2 = t2 = seMat(3+1, j); te2 = e2 = seMat(3+2, j);
      if(t2 < t1){
        tt1 = t2; tt2 = t1;
      }
      if(s2 < s1){
        ts1 = s2; ts2 = s1;
      }
      if(e2 < e1){
        te1 = e2; te2 = e1;
      }
      ns1 = ts1; ne2 = te2;
      if(te1 < ts2){
        ne1 = te1; ns2 = ts2;
      } else{
        if(ts2 > tt1){
          if(tt2 > te1){
            ne1 = ts2; ns2 = max(IntegerVector::create(ts2+1, te1));
          } else{
            ne1 = ts2; ns2 = ts2+1;
          }
        } else{
          if(tt2 > te1){
            ne1 = te1; ns2 = te1+1;
          } else{
            ne1 = (tt1+tt2)/2; ns2 = ne1+1;
          }
        }
      }
      n1 = ne1-ns1+1; n2 = ne2-ns2+1;
      resMat(l, 0) = ns1; resMat(l, 1) = ne1; resMat(l, 2) = ns2; resMat(l, 3) = ne2;
      if(n1>minLength && n2>minLength){
        for(k=0; k<maxLevel; k++){
          mean1 = mean2 = 0;
          for(t=ns1-1; t<ne1; t++){
            mean1 += yMat(t, k);
          }
          mean1 = mean1/n1;
          for(t=ns2-1; t<ne2; t++){
            mean2 += yMat(t, k);
          }
          mean2 = mean2/n2;
          resMat(l, 4+k) = (mean1-mean2)/var(k);
          if(fabs(resMat(l, 4+k)) > tmp) tmp = fabs(resMat(l, 4+k));
        }
        resMat(l, 4+maxLevel) = fabs(tmp);
      }
    }

  }
  return(Rcpp::List::create(Rcpp::Named("res")=resMat, Rcpp::Named("se")=seMat, Rcpp::Named("var")=var));

}

//' @keywords internal
// [[Rcpp::export]]
NumericMatrix funcResVar(NumericMatrix yMat, NumericMatrix seMat, NumericVector tmpVar){

	int i, k, t, ns1, ne1, ns2, ne2, n1, n2;
	int maxLevel = yMat.ncol();
	int Msq = seMat.nrow();
	double mean1, mean2;
	NumericMatrix res(Msq, maxLevel);

	for(i=0; i<Msq; i++){
		ns1 = seMat(i, 0); ne1 = seMat(i, 1); ns2 = seMat(i, 2); ne2 = seMat(i, 3);
		n1 = ne1-ns1+1; n2 = ne2-ns2+1;
		for(k=0; k<maxLevel; k++){
			mean1 = mean2 = 0;
			for(t=ns1-1; t<ne1; t++) mean1 = mean1 + yMat(t, k);
			mean1 = mean1/n1;
			for(t=ns2-1; t<ne2; t++) mean2 = mean2 + yMat(t, k);
			mean2 = mean2/n2;
			res(i, k) = (mean1-mean2)/tmpVar(k);
		}
	}
	return(res);

}

//' @keywords internal
// [[Rcpp::export]]
NumericMatrix funcApplyVar(NumericMatrix nullStatMat, int maxLevel, int Msq){

  int i, j, k;
  NumericMatrix bootVar(Msq, maxLevel);

  j = k = 0;
  for(i=0; i<Msq*maxLevel; i++){
    bootVar(j, k) = var(nullStatMat(_, i));
    bootVar(j, k) = sqrt(bootVar(j, k));
    j += 1;
    if(i==Msq*(k+1)-1){
      j = 0; k += 1;
    }
  }

  return(bootVar);
}


