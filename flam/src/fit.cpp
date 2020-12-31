#include <Rcpp.h>


using namespace Rcpp;

// Dynamic programming algorithm for the 1d fused lasso problem
// (Ryan's implementation of Nick Johnson's algorithm)


// [[Rcpp::export]]
NumericVector tf_dp(int n, NumericVector y, double lam) {
  
    NumericVector beta (n);
    
    // Take care of a few trivial cases
    if (n==0) return(beta);
    if (n==1 || lam==0) {
        for (int i=0; i<n; i++) beta[i] = y[i];
        return(beta);
    }
    
    // These are used to store the derivative of the
    // piecewise quadratic function of interest
    double afirst, alast, bfirst, blast;
    double *x = (double*)malloc(2*n*sizeof(double));
    double *a = (double*)malloc(2*n*sizeof(double));
    double *b = (double*)malloc(2*n*sizeof(double));
    int l,r;
    
    // These are the knots of the back-pointers
    double *tm = (double*)malloc((n-1)*sizeof(double));
    double *tp = (double*)malloc((n-1)*sizeof(double));
    
    // We step through the first iteration manually
    tm[0] = -lam+y[0];
    tp[0] = lam+y[0];
    l = n-1;
    r = n;
    x[l] = tm[0];
    x[r] = tp[0];
    a[l] = 1;
    b[l] = -y[0]+lam;
    a[r] = -1;
    b[r] = y[0]+lam;
    afirst = 1;
    bfirst = -lam-y[1];
    alast = -1;
    blast = -lam+y[1];
    
    // Now iterations 2 through n-1
    int lo, hi;
    double alo, blo, ahi, bhi;
    for (int k=1; k<n-1; k++) {
        // Compute lo: step up from l until the
        // derivative is greater than -lam
        alo = afirst;
        blo = bfirst;
        for (lo=l; lo<=r; lo++) {
            if (alo*x[lo]+blo > -lam) break;
            alo += a[lo];
            blo += b[lo];
        }
        
        // Compute the negative knot
        tm[k] = (-lam-blo)/alo;
        l = lo-1;
        x[l] = tm[k];
        
        // Compute hi: step down from r until the
        // derivative is less than lam
        ahi = alast;
        bhi = blast;
        for (hi=r; hi>=l; hi--) {
            if (-ahi*x[hi]-bhi < lam) break;
            ahi += a[hi];
            bhi += b[hi];
        }
        
        // Compute the positive knot
        tp[k] = (lam+bhi)/(-ahi);
        r = hi+1;
        x[r] = tp[k];
        
        // Update a and b
        a[l] = alo;
        b[l] = blo+lam;
        a[r] = ahi;
        b[r] = bhi+lam;
        afirst = 1;
        bfirst = -lam-y[k+1];
        alast = -1;
        blast = -lam+y[k+1];
    }
    
    // Compute the last coefficient: this is where
    // the function has zero derivative
    
    alo = afirst;
    blo = bfirst;
    for (lo=l; lo<=r; lo++) {
        if (alo*x[lo]+blo > 0) break;
        alo += a[lo];
        blo += b[lo];
    }
    beta[n-1] = -blo/alo;
    
    // Compute the rest of the coefficients, by the
    // back-pointers
    for (int k=n-2; k>=0; k--) {
        if (beta[k+1]>tp[k]) beta[k] = tp[k];
        else if (beta[k+1]<tm[k]) beta[k] = tm[k];
        else beta[k] = beta[k+1];
    }
    
    // Done! Free up memory
    free(x);
    free(a);
    free(b);
    free(tm);
    free(tp);
    
    return(beta);
}

// Soft-thresholding
void soft_thresh(int n, double *y, double lam, double *beta) {
    for (int i=0; i<n; i++) {
        if (y[i]>lam) beta[i] = y[i]-lam;
        else if (y[i]<-lam) beta[i] = y[i]+lam;
        else beta[i]=0;
    }
}

// [[Rcpp::export]]
void updatecolumn(NumericMatrix thetamat, NumericVector column, int colnum, int nrow) {
   
   for (int k=0; k<nrow; k++) {
        thetamat(k,colnum-1) = column[k];
   }
}

// [[Rcpp::export]]
void updatevector(NumericVector vector, NumericVector update, int length) {
   
   for (int k=0; k<length; k++) {
        vector[k] = update[k];
   }
}

// [[Rcpp::export]]
void updateresidual(NumericVector resid, NumericVector y, NumericMatrix thetamat, int j, int n, int p) {
   
   for (int k=0; k<n; k++) {
        resid[k] = y[k];
        if (p>1) {
        for (int l=0; l<p; l++) {
          resid[k] -= thetamat(k,l);
        }
        resid[k] += thetamat(k,j-1);
        }
   }
}

// [[Rcpp::export]]
double calcsum(NumericMatrix thetamat, int n, int p) {
   double sum=0;
   
   for (int k=1; k<n; k++) {
      for (int l=0; l<p; l++) {
          if (thetamat(k,l)-thetamat(k-1,l) > 0) sum += thetamat(k,l)-thetamat(k-1,l); else sum -= thetamat(k,l)-thetamat(k-1,l);
      }
   }
   return(sum);
}

// [[Rcpp::export]]
IntegerVector order_(NumericVector x) {
  NumericVector sorted = clone(x).sort();
  return match(sorted, x);
}

// [[Rcpp::export]]
double maxLambda_a1_C_single(NumericVector ycenter, NumericVector xcol, int n) {
   double lam = 0;
   double cumsum=0;
   IntegerVector index = order_(xcol);
    
   for (int k=0; k<(n-1); k++) {
          cumsum += ycenter[index[k]-1];
          if (lam < cumsum) lam = cumsum;
          if (lam < -cumsum) lam = -cumsum;
   }
   return(lam);
}

// [[Rcpp::export]]
NumericVector ordertheta(NumericMatrix thetamat, IntegerMatrix orderx, int i, int n) {

   NumericVector newcol (n);
    
   for (int k=0; k<n; k++) {
          newcol[k] = thetamat(orderx(k,i-1)-1,i-1);
    }
   return(newcol);
}

// [[Rcpp::export]]
NumericMatrix orderthetamat(NumericMatrix thetamat, IntegerMatrix orderx, int n, int p) {
   NumericMatrix newthetamat (n,p);
    
   for (int k=0; k<n; k++) {
     for (int l=0; l<p; l++) {
          newthetamat(k,l) = thetamat(orderx(k,l)-1,l);
     }
    }
   return(newthetamat);
}


// [[Rcpp::export]]
List flamstep(NumericMatrix initialthetamat, NumericVector y, double lambda, double alpha, int n, int p, IntegerMatrix orderx, IntegerMatrix rankx, double tolerance) {
   NumericMatrix thetamat = initialthetamat;
  double oldobj=0; double newobj=0; int niter=0; int converge=0; double dev=0; double sqerr=0; double sum1=0; double sum2=0; double err=0; double colnorm=0; double beta0=0;
  NumericVector r (n); NumericVector betasol (n); double sumsqest=0;
  NumericVector est (n); NumericVector rr (n); NumericVector update (n);
  double sumest=0; 
  
  while (converge==0 && niter<1000) {
    ++niter;
    
    for (int j=1; j<p+1; j++) {
    beta0 = 0;
			updateresidual(r, y, thetamat, j, n, p);
			
    for (int k=0; k<n; k++) {
        rr[k] = r[orderx(k,j-1)-1];
      }
      
			betasol = tf_dp(n, rr, alpha*lambda);
    for (int k=0; k<n; k++) {
        est[k] = betasol[rankx(k,j-1)-1];
      }
			
      sumest = 0; sumsqest = 0;
      for (int k=0; k<n; k++) {
        sumest += est[k];
      }
			beta0 += sumest/n;
			for (int k=0; k<n; k++) {
        est[k] -= sumest/n;
			}
     for (int k=0; k<n; k++) {
        sumsqest += pow(est[k],2.0);
      }
			if (sumsqest!=0) {
        double scalef = 1 - (1 - alpha) * lambda/pow(sumsqest,0.5);
        
        if (scalef < 0) {
        for (int k=0; k<n; k++) {
          update[k] = 0;}
        } else {
        for (int k=0; k<n; k++) {
        update[k] = scalef * est[k];
        }
        }
			  updatecolumn(thetamat,update,j,n);
      } else {
        updatecolumn(thetamat,est,j,n); 
			}
    }
    
    NumericMatrix thetamatreord = orderthetamat(thetamat, orderx, n, p);
    
    sum1 = calcsum(thetamatreord, n, p);
		
    sum2 = 0;
      for (int l=0; l<p; l++) { 
       colnorm = 0;
   for (int k=0; k<n; k++) {
          colnorm += pow(thetamat(k,l),2.0);
     }
     sum2 += pow(colnorm,0.5);
    }
		
    sqerr = 0;
     for (int k=0; k<n; k++) {
       err = y[k] - beta0;
     for (int l=0; l<p; l++) {
          err -= thetamat(k,l);
     }
     sqerr += pow(err,2.0);
    }
 		oldobj = newobj;
    newobj = 0.5 * sqerr + lambda * alpha * sum1 + lambda * (1-alpha) * sum2;

		if (niter>1) {
			dev = (oldobj - newobj) / oldobj;
			if (dev < tolerance) converge = 1;
		} 

  }
   return List::create(
     _["thetamat"] = thetamat,
     _["beta0"] = beta0
   );
}

// [[Rcpp::export]]
List flamsteplogistic(NumericMatrix initialthetamat, NumericVector y, double lambda, double alpha, int n, int p, IntegerMatrix orderx, IntegerMatrix rankx, double tolerance) {

  NumericMatrix thetamat = initialthetamat;
  double oldobj=0; double newobj=0; int niter=0; int converge=0; double dev=0; 
  double sum1=0; double sum2=0; double colnorm=0;  double beta0=0;
  NumericVector r (n); NumericVector betasol (n); double sumsqest=0;
  NumericVector rr (n); NumericVector update (n);
  NumericVector yhat (n);
  double ftot=0;
  
  double L = (p + 1)/4;
    for (int k=0; k<n; k++) {
    yhat[k] = beta0;
    for (int l=0; l<p; l++) {
      yhat[k] += thetamat(k,l);
    }
  }
  
  while (converge==0 && niter<1000) {

    ++niter;
    
    for (int k=0; k<n; k++) {
      beta0 += (y[k] - exp(yhat[k])/(1 + exp(yhat[k])))/(n * L);
    }
  
   for (int j=1; j<p+1; j++) {
    
    for (int k=0; k<n; k++) {  
       r[k] = thetamat(k,j-1) + (y[k] - exp(yhat[k])/(1 + exp(yhat[k])))/L;
    }
    for (int k=0; k<n; k++) {
        rr[k] = r[orderx(k,j-1)-1];
      }

  		betasol = tf_dp(n, rr, alpha*lambda/L);
    
    for (int k=0; k<n; k++) {
        thetamat(k,j-1) = betasol[rankx(k,j-1)-1];
      }
		}
   for (int j=1; j<p+1; j++) {
     sumsqest = 0;
       for (int k=0; k<n; k++) {
        sumsqest += pow(thetamat(k,j-1),2.0);
      }
  		if (sumsqest!=0) {
        double scalef = 1 - (1 - alpha) * lambda/(L * pow(sumsqest,0.5));
        
        if (scalef < 0) {
        for (int k=0; k<n; k++) {
          thetamat(k,j-1) = 0;}
        } else {
        for (int k=0; k<n; k++) {
        thetamat(k,j-1) *= scalef;
        }
        }
      }  
  }
    
    NumericMatrix thetamatreord = orderthetamat(thetamat, orderx, n, p);

    sum1 = calcsum(thetamatreord, n, p);	
    sum2 = 0;
      for (int l=0; l<p; l++) { 
       colnorm = 0;
   for (int k=0; k<n; k++) {
          colnorm += pow(thetamat(k,l),2.0);
     }
     sum2 += pow(colnorm,0.5);
    }


  for (int k=0; k<n; k++) {
    yhat[k] = beta0;
    for (int l=0; l<p; l++) {
      yhat[k] += thetamat(k,l);
    }
  }
  
  ftot = 0;
  for (int k=0; k<n; k++) {
    ftot += -y[k] * yhat[k] + log(1 + exp(yhat[k]));
  }
 		oldobj = newobj;
    newobj = ftot + lambda * alpha * sum1 + lambda * (1-alpha) * sum2;

		if (niter>1) {
			dev = (oldobj - newobj) / oldobj;
			if (dev < tolerance) converge = 1;
		} 

  }
  return List::create(
    _["thetamat"] = thetamat,
    _["beta0"] = beta0
  );
}

