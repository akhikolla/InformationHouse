
#include "fdrseg.h"

// [[Rcpp::export(".smuce_cpp")]]
List smuce_cpp(NumericVector Y, double q, double sd)
{
    int i, j, cnt;
    int nseg;
    double rad, aux_mu, aux_c, aux_m, lstar;
    
    int n = Y.size();
    
    NumericVector S(n+1);   // cumulative sum
    S[0] = 0;
    for (i = 0; i < n; ++i)
        S[i+1] = S[i] + Y[i];
    
    IntegerVector leE(n);   // leftmost end points
    NumericVector mu(n);    // optimal value of step function
    NumericVector C(n+1);   // optimal cost
    IntegerVector L(n+1);   // left bound of jump location
    IntegerVector R(n+1);   // right bound of jump location
    NumericVector loB(n+1); // auxilary memory for lower bounds
    NumericVector upB(n+1); // auxilary memory for upper bounds
    
    // initialization
    nseg = 0;
    R[0] = 0;
    R[1] = 0;
    L[0] = 0;
    for (i = 0; i < n+1; ++i)
        loB[i] = -INFINITY;
    for (i = 0; i < n+1; ++i)
        upB[i] = INFINITY;
    for (i = 0; i < n; ++i)
        C[i+1] = INFINITY;
    C[0] = 0.;
    
    // search rightward
    while (isinf(C[n])) {
        ++nseg;
        lstar = R[nseg-1];
        for (i = R[nseg]; i < n; ++i) {
            for (j = i; j > R[nseg]; --j) {
                aux_m  = (S[i+1] - S[j]) / double(i-j+1);
                rad    = pens(n, i-j+1, q) * sd / sqrt(i-j+1);
                loB[j] = fmax(fmax(aux_m - rad, loB[j]), loB[j+1]);
                upB[j] = fmin(fmin(aux_m + rad, upB[j]), upB[j+1]);
            }
            for (j = R[nseg]; j >= lstar; --j) {
                aux_m  = (S[i+1] - S[j]) / double(i-j+1);
                rad    = pens(n, i-j+1, q) * sd / sqrt(i-j+1);
                loB[j] = fmax(fmax(aux_m - rad, loB[j]), loB[j+1]);
                upB[j] = fmin(fmin(aux_m + rad, upB[j]), upB[j+1]);
                if (loB[j] < upB[j]) {
                    if (aux_m > upB[j])
                        aux_mu = upB[j];
                    else if (aux_m < loB[j])
                        aux_mu = loB[j];
                    else
                        aux_mu = aux_m;
                    aux_c = (aux_mu*aux_mu - 2*aux_mu*aux_m)*(i-j+1) + C[j];
                    if (aux_c < C[i+1]) {
                        C[i+1] = aux_c;
                        mu[i] = aux_mu;
                        leE[i] = j + 1;
                    }
                } else {
                    lstar = j + 1;
                    break;
                }
            }
            if (i == R[nseg])
                L[nseg] = j + 1 + (loB[j] >= upB[j]);
            if (isinf(C[i+1]))
                break;
        }
        R[nseg+1] = i;
        L[nseg+1] = R[nseg];
    }
    NumericVector value(nseg);
    IntegerVector left(nseg);
    
    cnt = n-1;
    for (i = nseg-1; i >= 0; --i) {
        left[i] = leE[cnt];
        value[i] = mu[cnt];
        cnt = leE[cnt] - 2;
    }
    
    return List::create(Named("value") = value, Named("left") = left, Named("n") = n);
}
