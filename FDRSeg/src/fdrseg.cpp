
#include "fdrseg.h"

// [[Rcpp::export(".fdrseg_cpp")]]
List fdrseg_cpp(NumericVector Y, NumericVector q, NumericVector qm, double sd)
{
    int i, j, k, cnt;
    int isstop, njmp, jmin;
    double lb, ub, rad, aux_m, aux_mu, aux_c;
    int lbi, ubi;
    int imin, aux_imin, imax, aux_imax;
    int n = Y.size();
    NumericVector S(n+1); // cumulative sum
    NumericVector lm(n);  // lm[i]: maximal mean on intervals of length (i+1) within some large interval --> lower bound
    NumericVector um(n);  // um[i]: minimal mean on intervals of length (i+1) within some large interval --> upper bound
    IntegerVector L(n);   // leftmost end points
    NumericVector mu(n);  // optimal value of step function
    NumericVector C(n);   // optimal cost
    NumericVector J(n);   // number of jumps
    
    // initialization (separate loop favors memory continuity)
    S[0] = 0;
    for (i = 0; i < n; ++i)
        S[i+1] = S[i] + Y[i];
    for (i = 0; i < n; ++i)
        C[i] = INFINITY;
    for (i = 0; i < n; ++i)
        lm[i] = -INFINITY;
    for (i = 0; i < n; ++i)
        um[i] = INFINITY;
    for (i = 0; i < n; ++i)
        J[i] = n;
    
    imax = 0;
    // first search
    for (i = 0; i < n; ++i) {
        lb = -INFINITY;
        ub = INFINITY;
        for (j = 0; j <= i; ++j) {
            aux_m = (S[i+1] - S[i-j]) / double(j+1);
            // maximal mean on intervals of length (j+1) within [0, i]
            lm[j] = fmax(aux_m, lm[j]);
            // minimal mean on intervals of length (j+1) within [0, i]
            um[j] = fmin(aux_m, um[j]);
            rad = penfs(i+1, j+1, q) * sd / sqrt(j+1);
            if (lm[j] - rad > lb) {
                lbi = j;
                lb = lm[j] - rad;
            }
            if (um[j] + rad < ub) {
                ubi = j;
                ub = um[j] + rad;
            }
        }
        if (lb < ub) {
            if (aux_m > ub)
                mu[i] = ub;
            else if (aux_m < lb)
                mu[i] = lb;
            else
                mu[i] = aux_m;
            if (i > imax)
                imax = i;
            C[i] = (mu[i] * mu[i] - 2 * mu[i] * aux_m)*(i+1.);
            L[i] = 1;
            J[i] = 0;
        } else if (lm[lbi] - um[ubi] > eta(lbi+1, ubi+1, n, sd, qm))
            break;
    }
    
    // later searches
    imin = 0;
    njmp = 0;
    while (isinf(C[n-1]) && njmp < n ) {
        jmin = n;
        aux_imin = n;
        aux_imax = 0;
        for (i = imin; i < n; ++i) {
            if (J[i] == njmp) {
                if (i < jmin)
                    jmin = i;
            } else if (isinf(C[i])) {
                isstop = 0;
                for (j = 0; j < i-jmin; ++j)
                    lm[j] = -INFINITY;
                for (j = 0; j < i-jmin; ++j)
                    um[j] = INFINITY;
                for (j = i-1; j >= jmin; --j) {
                    for (k = 0; k < i-j; ++k) {
                        aux_m = (S[j+2+k]- S[j+1]) / double(k+1.);
                        // maximal mean on intervals of length (k+1) within [j+1, i]
                        lm[k] = fmax(aux_m, lm[k]);
                        // minimal mean on intervals of length (k+1) within [j+1, i]
                        um[k] = fmin(aux_m, um[k]);
                    }
                    if (J[j] == njmp) {
                        lb = -INFINITY;
                        ub = INFINITY;
                        for (k = 0; k < i-j; ++k) {
                            rad = penfs(i-j, k+1, q) * sd / sqrt(k+1.);
                            if (lm[k] - rad > lb) {
                                lbi = k;
                                lb = lm[k] - rad;
                            }
                            if (um[k] + rad < ub) {
                                ubi = k;
                                ub = um[k] + rad;
                            }
                        }
                        if (lb < ub) {
                            if (aux_m > ub)
                                aux_mu = ub;
                            else if (aux_m < lb)
                                aux_mu = lb;
                            else
                                aux_mu = aux_m;
                            if (i < aux_imin)
                                aux_imin = i;
                            if (i > aux_imax)
                                aux_imax = i;
                            aux_c = (aux_mu*aux_mu - 2*aux_mu*aux_m)*(i-j) + C[j];
                            if (aux_c < C[i]) {
                                C[i]  = aux_c;
                                mu[i] = aux_mu;
                                L[i]  = j+2;
                                J[i]  = njmp + 1;
                            }
                        } else {
                            if (lm[lbi] - um[ubi] > eta(lbi+1, ubi+1, n-jmin, sd, qm) && j >= imax)
                                isstop = 1;
                            if (lm[lbi] - um[ubi] > eta(lbi+1, ubi+1, i-jmin, sd, qm))
                                break;
                        }
                    }
                }
                if (isinf(C[i]) && isstop)
                    break;
            }
        }
        imin = aux_imin;
        imax = aux_imax;
        
        ++njmp;
    }
    
    NumericVector value(njmp+1);
    IntegerVector left(njmp+1);
    cnt = n-1;
    for (i = njmp; i >= 0; --i) {
        left[i] = L[cnt];
        value[i] = mu[cnt];
        cnt = L[cnt] - 2;
    }
    
    return List::create(Named("value") = value, Named("left") = left, Named("n") = n);
}
