#include "fdrseg.h"

double sign(double a) {return (a == 0) ? 0 : (a < 0 ? -1 : 1); }

// [[Rcpp::export(".scalar_quantile")]]
double scalar_quantile(double p, int n, int r)
{
    static List simul; 
    IntegerVector narray;
    IntegerVector rarray;
    List darray;
        
    int isdone = 0;
    int i, j, k;
    if (simul.size() > 0) {
        narray = simul[0];
        rarray = simul[1];
        darray = simul[2];
        for (i = 0; i < narray.size(); ++i)
            if (n == narray[i] && r <= rarray[i]) {
                isdone = 1;
                r = rarray[i];
                break;
            }
    }
    NumericVector data(r);

    if (isdone) {
        data = darray[i];
    } else { // simulation
        double sm;
        NumericVector Y(n);
        NumericVector S(n+1);
        GetRNGstate();
        for (i = 0; i < r; ++i) {
            Y = rnorm(n, 0, 1);
            S[0] = 0;
            for (j = 0; j < n; ++j) {
                S[j+1] = S[j] + Y[j];
            }
            data[i] = -INFINITY;
            for (j = 0; j < n; ++j) {
                sm = -INFINITY;
                for (k = 0; k < n-j; ++k) {
                    sm = fmax(fabs(S[k+j+1] - S[k]), sm);
                }
                data[i] = fmax(data[i], sm/sqrt(j+1.)-sqrt(2+2*log(n/(j+1.))));
            }
        }
        PutRNGstate();
        sort(data.begin(), data.end());
        
        narray = Language("c", narray, n).eval();
        rarray = Language("c", rarray, r).eval();
        darray = Language("c", darray, List::create(data)).eval();
        simul = List::create(narray, rarray, darray);
    }
    return data[round(fmin(fmax(p*r-1,0),r-1))];
}

// [[Rcpp::export(".vector_quantile")]]
NumericVector vector_quantile(double p, int n, int r)
{
    NumericVector q(n); // output
    
    // keep simple, the list will not increase
    static List simul;  // n, r, data
        
    int i, j, k;
    int isdone = 0;
    if (simul.size() > 0) {
        int n0  = simul[0];
        int r0  = simul[1];
        if (n <= n0 && r <= r0) {
            n = n0;
            r = r0;
            isdone = 1;
        }
    } 
    
    NumericVector data(r * n);
    NumericVector qntl(r);
    if (isdone) {
        data = simul[2];
    } else { // simulation
        NumericVector Y(n);
        NumericVector S(n+1);
        
        double mY, mYo;
        NumericVector sma(n);
        NumericVector smb(n);
        GetRNGstate();
        for (i = 0; i < r; ++i) {
            Y = rnorm(n, 0, 1);
            S[0] = 0;
            mYo  = 0;
            for (j = 0; j < n; ++j) {
                S[j+1] = S[j] + Y[j];
                mY     = S[j+1] / (j+1.);
                sma[j]      = -INFINITY;
                smb[j]      = -INFINITY;
                data[i*n+j] = -INFINITY;
                for (k = 0; k <= j; ++k) {
                    sma[k] = fmax( (S[j+1] - S[j-k] - (k+1.)*mY)/sqrt(k+1.), sqrt(k+1.)*(mYo-mY) + sma[k]);
                    smb[k] = fmax(-(S[j+1] - S[j-k] - (k+1.)*mY)/sqrt(k+1.), sqrt(k+1.)*(mY-mYo) + smb[k]);
                    data[i*n+j] = fmax(data[i*n+j], fmax(sma[k], smb[k]) - sqrt(2+2*log((j+1.)/(k+1.))));
                }
                mYo = mY;
            }
        }
        PutRNGstate();
        // store simulated results
        simul = List::create(n, r, data);
    }
    
    for (i = 0; i < q.size(); ++i) { // note: n might > q.size()
        for (j = 0; j < r; ++j) {
            qntl[j] = data[j*n+i];
        }
        sort(qntl.begin(), qntl.end());
        q[i] = qntl[round(fmin(fmax(p*r-1,0),r-1))];
    }
    
    return q;
}

// multiresolution statistics up to scale Y.size()
// [[Rcpp::export(".mrstatvec_cpp")]]
NumericVector mrstatvec_cpp(NumericVector Y) // int id, double param, int penId, int fmId
{
    int j, k;
    int n = Y.size();
    NumericVector mr(n);              // output
    NumericVector S(n+1);             // cumulative sums
    
    double mY, mYo;
    NumericVector sma(n);
    NumericVector smb(n);
    S[0] = 0;
    mYo  = 0;
    for (j = 0; j < n; ++j) {
      S[j+1] = S[j] + Y[j];
      mY     = S[j+1] / (j+1.);
      sma[j] = -INFINITY;
      smb[j] = -INFINITY;
      mr[j]  = -INFINITY;
      for (k = 0; k <= j; ++k) {
        sma[k] = fmax( (S[j+1] - S[j-k] - (k+1.)*mY)/sqrt(k+1.), sqrt(k+1.)*(mYo-mY) + sma[k]);
        smb[k] = fmax(-(S[j+1] - S[j-k] - (k+1.)*mY)/sqrt(k+1.), sqrt(k+1.)*(mY-mYo) + smb[k]);
        mr[j]  = fmax(mr[j], fmax(sma[k], smb[k]) - sqrt(2. + 2*log((j+1.)/double(k+1.))));
      }
      mYo = mY;
    }
    
    return mr;
}

double pens(int n, int lens, double q)
{
    return fmax(1e-3, q + sqrt(2 + 2*log(n/double(lens))));
}

double penfs(int n, int len, NumericVector q)
{
    return fmax(1e-3, q[n-1] + sqrt(2. + 2*log(n/double(len))));
}

double eta(int llen, int ulen, int n, double sd, NumericVector q)
{
    return sd * (penfs(n, llen, q)/sqrt(llen) + penfs(n, ulen, q)/sqrt(ulen));
}

