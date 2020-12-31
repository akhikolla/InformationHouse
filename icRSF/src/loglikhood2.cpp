#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double loglikC(NumericVector parm, NumericMatrix Dmat, NumericVector x) {
  int nsub = Dmat.nrow(), J = Dmat.ncol() - 1, i, j;
  double result = 0, temp;
  for (i = 0; i < nsub; ++i) {
    temp = Dmat(i, 0);
  	for (j = 1; j < J+1; ++j) {
  		temp += Dmat(i, j)*exp(-exp(x[i]*parm[0] + parm[j]));
  	}
  	result += log(temp);
  }
  result = -result;
  return result;   
}

// [[Rcpp::export]]
NumericVector gradlikC(NumericVector parm, NumericMatrix Dmat, NumericVector x) {
  int nsub = Dmat.nrow(), J = Dmat.ncol() - 1, i, j;
  NumericVector Deriv(J + 1);
  double denom, temp;
  for (i = 0; i < nsub; ++i) {
  	denom = Dmat(i, 0);
  	for (j = 1; j < J+1; ++j) {
  	  denom += Dmat(i, j)*exp(-exp(x[i]*parm[0] + parm[j]));
  	} 
  	for (j = 1; j < J+1; ++j) {
  	  temp = Dmat(i, j)*exp(-exp(x[i]*parm[0] + parm[j]))*exp(x[i]*parm[0] + parm[j]);
  	  Deriv[0] += temp*x[i]/denom;
  	  Deriv[j] += temp/denom;	
  	}
  }
  return Deriv;
}



// [[Rcpp::export]]
double loglikCD(NumericVector parm, NumericMatrix Dmat, NumericVector x) {
  int nsub = Dmat.nrow(), J = Dmat.ncol() - 1, i, j;
  double result = 0, temp;
  NumericVector parm1(J+1);
  parm1[0] = parm[0];
  parm1[1] = parm[1];
  for (i = 2; i < J +1; i++) {
    parm1[i] = parm1[i-1] + parm[i];
  }
  for (i = 0; i < nsub; ++i) {
    temp = Dmat(i, 0);
  	for (j = 1; j < J+1; ++j) {
  		temp += Dmat(i, j)*exp(-exp(x[i]*parm[0] + parm1[j]));
  	}
  	result += log(temp);
  }
  result = -result;
  return result;   
}

// [[Rcpp::export]]
NumericVector gradlikCD(NumericVector parm, NumericMatrix Dmat, NumericVector x) {
  int nsub = Dmat.nrow(), J = Dmat.ncol() - 1, i, j;
  NumericVector Deriv(J + 1), parm1(J+1);
  double denom, temp;
  parm1[0] = parm[0];
  parm1[1] = parm[1];
  for (i = 2; i < J +1; i++) {
    parm1[i] = parm1[i-1] + parm[i];
  }
  for (i = 0; i < nsub; ++i) {
  	denom = Dmat(i, 0);
  	for (j = 1; j < J+1; ++j) {
  	  denom += Dmat(i, j)*exp(-exp(x[i]*parm1[0] + parm1[j]));
  	} 
  	for (j = 1; j < J+1; ++j) {
  	  temp = Dmat(i, j)*exp(-exp(x[i]*parm1[0] + parm1[j]))*exp(x[i]*parm1[0] + parm1[j]);
  	  Deriv[0] += temp*x[i]/denom;
  	  Deriv[j] += temp/denom;	
  	}
  }
  for (j=J-1; j>0; j--) {
    Deriv[j] = Deriv[j] + Deriv[j+1];
  }
  return Deriv;
}

// [[Rcpp::export]]
double loglikCD0(NumericVector parm, NumericMatrix Dmat) {
  int nsub = Dmat.nrow(), J = Dmat.ncol() - 1, i, j;
  double result = 0, temp;
  NumericVector parm1(J);
  parm1[0] = parm[0];
  for (i = 1; i < J; i++) {
    parm1[i] = parm1[i-1] + parm[i];
  }
  for (i = 0; i < nsub; ++i) {
    temp = Dmat(i, 0);
  	for (j = 0; j < J; ++j) {
  		temp += Dmat(i, j+1)*exp(-exp(parm1[j]));
  	}
  	result += log(temp);
  }
  result = -result;
  return result;   
}


// [[Rcpp::export]]
double loglikC0(NumericVector parm, NumericMatrix Dmat) {
  int nsub = Dmat.nrow(), J = Dmat.ncol() - 1, i, j;
  double result = 0, temp;
  for (i = 0; i < nsub; ++i) {
    temp = Dmat(i, 0);
  	for (j = 0; j < J; ++j) {
  		temp += Dmat(i, j+1)*exp(-exp(parm[j]));
  	}
  	result += log(temp);
  }
  result = -result;
  return result;   
}

// [[Rcpp::export]]
NumericVector gradlikCD0(NumericVector parm, NumericMatrix Dmat) {
  int nsub = Dmat.nrow(), J = Dmat.ncol() - 1, i, j;
  NumericVector Deriv(J), parm1(J);
  double denom, temp;
  parm1[0] = parm[0];
  for (i = 1; i < J; i++) {
    parm1[i] = parm1[i-1] + parm[i];
  }
  for (i = 0; i < nsub; ++i) {
  	denom = Dmat(i, 0);
  	for (j = 0; j < J; ++j) {
  	  denom += Dmat(i, j+1)*exp(-exp(parm1[j]));
  	} 
  	for (j = 0; j < J; ++j) {
  	  temp = Dmat(i, j+1)*exp(-exp(parm1[j]))*exp(parm1[j]);
  	  Deriv[j] += temp/denom;	
  	}
  }
  for (j=J-2; j>=0; j--) {
    Deriv[j] = Deriv[j] + Deriv[j+1];
  }
  return Deriv;
}

// [[Rcpp::export]]
NumericVector gradlikC0(NumericVector parm, NumericMatrix Dmat) {
  int nsub = Dmat.nrow(), J = Dmat.ncol() - 1, i, j;
  NumericVector Deriv(J);
  double denom, temp;
  for (i = 0; i < nsub; ++i) {
  	denom = Dmat(i, 0);
  	for (j = 0; j < J; ++j) {
  	  denom += Dmat(i, j+1)*exp(-exp(parm[j]));
  	} 
  	for (j = 0; j < J; ++j) {
  	  temp = Dmat(i, j+1)*exp(-exp(parm[j]))*exp(parm[j]);
  	  Deriv[j] += temp/denom;	
  	}
  }
  return Deriv;
}


// [[Rcpp::export]]
NumericVector splitpointC(NumericMatrix Dm, NumericVector x, Function f) {
  int J = Dm.ncol() - 1;
  NumericVector ux = unique(x).sort(), chisq(J+2), chisqtemp(J+3, -1.0);
  int nux = ux.size();
  for (int i=0; i <= nux-2; i++) {
    chisq = f(Dm, x >= (ux[i]+ux[i+1])/2);
    if (chisq[0] >= chisqtemp[1]) {
       chisqtemp[0] = (ux[i]+ux[i+1])/2;
       for (int j=1; j<J+3; j++) {
         chisqtemp[j] = chisq[j-1];
       }
    }  
  }
  return chisqtemp; 
}
// [[Rcpp::export]]
NumericVector splitpointCD(NumericMatrix Dm, NumericVector x, Function f) {
  int J = Dm.ncol() - 1;
  NumericVector chisq(J+2), chisqtemp(J+3, -1.0);
  double ux = 0.5;
 
    chisq = f(Dm, x >= ux);
       chisqtemp[0] = ux;
       for (int j=1; j<J+3; j++) {
         chisqtemp[j] = chisq[j-1];
       }

  return chisqtemp; 
}
// [[Rcpp::export]]
List bestsplitC(NumericMatrix Dm, NumericMatrix Xmat, Function f){
  int J = Dm.ncol() - 1, nsub = Dm.nrow(), nbeta = Xmat.ncol();
  NumericVector result(J+3), temp(J+3), xtemp(nsub); int id=-1;
  result[1]=-100;
  for (int i=0; i <nbeta; i++) {
    xtemp = Xmat(_, i);
    temp = splitpointC(Dm, xtemp, f);
    if(temp[1]>=result[1]){
      result[0] = temp[0];
      result[1] = temp[1];
      id = i;
      for(int j=2; j<= J+2; j++)
      result[j] = temp[j];
    }
  }
  return List::create(id + 1, result);
}
// [[Rcpp::export]]
List bestsplitCD(NumericMatrix Dm, NumericMatrix Xmat, Function f) {
  int J = Dm.ncol() - 1, nsub = Dm.nrow(), nbeta = Xmat.ncol();
  NumericVector result(J+2), temp(J+2), xtemp(nsub); int id=-1;
  result[0]=-100;
  for (int i=0; i <nbeta; i++) {
    xtemp = Xmat(_, i);
    temp = f(Dm, xtemp);
    if(temp[0]>=result[0]){
      result[0] = temp[0];
      id = i;
      for(int j=1; j<= J+1; j++)
      result[j] = temp[j];
    }
  }
  return List::create(id + 1, result);
}


// [[Rcpp::export]]
double scorefun(double beta, NumericVector x, NumericVector parm, NumericMatrix Dmat) {
  int J = Dmat.ncol() - 1, nsub = Dmat.nrow(), i, j, j1;
  double score = 0.0, denom, info = 0.0;
  NumericVector y(J), dy(J), ly(J), jacob(J), dfdy(J);
  for (i = 0; i < nsub; i++) {
    denom = Dmat(i, 0);
    y = exp(-exp(parm + beta*x[i]));
    ly = exp(parm + beta*x[i]);
    jacob = -y*ly*x[i];
    for (j=0; j< J; j++) {
      denom += Dmat(i, j+1)*y[j];
    }
    for (j = 0; j < J; j++) {
      dfdy[j] = Dmat(i, j+1)/denom;
      score += jacob[j]*dfdy[j];
    }
    for (j = 0; j < J; j++) {
      for (j1 = 0; j1 <J; j1++) {
        info += -jacob[j]*jacob[j1]*dfdy[j]*dfdy[j1];
      }
    }
    info += sum(y*ly*ly*dfdy - y*ly*dfdy)*x[i]*x[i];
  }
  return -score*score/info;
} 

// [[Rcpp::export]]
double scorefun0(NumericVector x, NumericVector parm, NumericMatrix Dmat) {
  int J = Dmat.ncol() - 1, nsub = Dmat.nrow(), i, j, j1;
  double score = 0.0, denom, info = 0.0;
  NumericVector y(J), dy(J), ly(J), jacob(J), dfdy(J);
  for (i = 0; i < nsub; i++) {
    denom = Dmat(i, 0);
    y = exp(-exp(parm));
    ly = exp(parm);
    jacob = -y*ly*x[i];
    for (j=0; j< J; j++) {
      denom += Dmat(i, j+1)*y[j];
    }
    for (j = 0; j < J; j++) {
      dfdy[j] = Dmat(i, j+1)/denom;
      score += jacob[j]*dfdy[j];
    }
    for (j = 0; j < J; j++) {
      for (j1 = 0; j1 <J; j1++) {
        info += -jacob[j]*jacob[j1]*dfdy[j]*dfdy[j1];
      }
    }
    info += sum(y*ly*ly*dfdy - y*ly*dfdy)*x[i]*x[i];
  }
  return -score*score/info;
}

// [[Rcpp::export]]
NumericVector splitpt(NumericMatrix Dmat, NumericVector x, NumericVector parm) {
  int nsub = Dmat.nrow(), j;
  NumericVector ux = unique(x).sort(), result(2, -1.0), xtemp(nsub);
  double chisq;
  int nux = ux.size();
  for (int i=0; i < nux-1; i++) {
    for (j=0; j<nsub; j++) {
      if (x[j] > ux[i]) {
        xtemp[j] = 1.0;
      } else {
        xtemp[j] = -1.0;
      }
    }
    chisq = scorefun0(xtemp, parm, Dmat);
    if (chisq >= result[1]) {
       result[0] = ux[i];
	   result[1] = chisq;
    }  
  }
  return result;   
}

// [[Rcpp::export]]
List bsplitC(NumericMatrix Dmat, NumericMatrix Xmat, NumericVector parm) {
  int nsub = Dmat.nrow(), nbeta = Xmat.ncol(), id = -1;
  NumericVector result(2, -9999999999), temp(2), xtemp(nsub);
  for (int i=0; i < nbeta; i++) {
    xtemp = Xmat(_, i);
    temp = splitpt(Dmat, xtemp, parm);
    if (temp[1] >= result[1]) {
      id = i;
      result[0] = temp[0];
      result[1] = temp[1];
    }
  }
  return List::create(id + 1, result);
}


