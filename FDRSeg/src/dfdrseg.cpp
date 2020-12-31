
#include "fdrseg.h"

// cost is computed on whole interval

// [[Rcpp::export(".dfdrseg_cpp")]]
List dfdrseg_cpp(NumericVector Y, NumericVector q, NumericVector qm, double sd, int lag) // int penId, int isconf
{    
    int i, j, k, cnt;
    int isstop, njmp, jmin;
    double lb, ub, rad, aux_m, aux_mu, aux_c;
    int lbi, ubi;
    int imin, aux_imin, imax, aux_imax;
    int n = Y.size();
    NumericVector S(n+1);  // cumulative sum
    NumericVector Sq(n+1); // cumulative sum of squares
    NumericVector lm(n);   // lm[i]: maximal mean on intervals of length (i+1) within some large interval --> lower bound
    NumericVector um(n);   // um[i]: minimal mean on intervals of length (i+1) within some large interval --> upper bound
    IntegerVector L(n);    // leftmost end points (number in R style)
    NumericVector mu(n);   // optimal value of step function
    NumericVector C(n);    // optimal cost
    NumericVector J(n);    // number of jumps
    
    if (n <= lag) {
        Rprintf("Error in JSMURF-FDR: the length of data should > the length of kernel! \n");
        return -1;
    }
    
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
    Sq[0] = 0;
    for (i = 0; i < n; ++i) {
        Sq[i+1] = Sq[i] + (Y[i] * Y[i]);
    }
    
    imax = 0;
    // first search
    for (i = 0; i < n; ++i) {
        lb = -INFINITY;
        ub = INFINITY;
        for (j = 0; j <= i-lag+1; ++j) { // on interval [i-j, i] the filtered signal is constant
            aux_m = (S[i+1] - S[i-j]) / double(j + 1.);
            // maximal mean on intervals of length (j+1) within [0, i]
            lm[j] = fmax(aux_m, lm[j]);
            // minimal mean on intervals of length (j+1) within [0, i]
            um[j] = fmin(aux_m, um[j]);
            rad   = penfs(i-lag+2, j+1, q) * sd / sqrt(j+1.);
            if (lm[j] - rad > lb) {
                lbi = j;
                lb  = lm[j] - rad;
            }
            if (um[j] + rad < ub) {
                ubi = j;
                ub  = um[j] + rad;
            }
        }
        if (lb < ub) {
            aux_m = S[i+1] / double(i+1.);
            if (aux_m > ub)
                mu[i] = ub;
            else if (aux_m < lb)
                mu[i] = lb;
            else
                mu[i] = aux_m;
            if (i > imax)
                imax = i;
            C[i] = (mu[i] * mu[i] - 2 * mu[i] * aux_m) * (i+1.) + Sq[i+1];
            L[i] = 1;
            J[i] = 0;
        } else if (lm[lbi] - um[ubi] > eta(lbi+1, ubi+1, n-lag+1, sd, qm))
            break;
    }
    
    // later searches
    imin = 0; // lag - 1
    njmp = 0;
    while (isinf(C[n-1]) && njmp < n ) {
        jmin     = n;
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
                    for (k = 0; k < i-j-lag+1; ++k) {
                        aux_m = (S[j+lag+k+1]- S[j+lag]) / double(k + 1.);
                        // maximal mean on intervals of length (k+1) within [j+lag, i]
                        lm[k] = fmax(aux_m, lm[k]);
                        // minimal mean on intervals of length (k+1) within [j+lag, i]
                        um[k] = fmin(aux_m, um[k]);
                    }
                    if (J[j] == njmp) {
                        lb = -INFINITY;
                        ub = INFINITY;
                        for (k = 0; k < i-j-lag+1; ++k) {
                            rad = penfs(i-j-lag+1, k+1, q) * sd / sqrt(k+1.);
                            if (lm[k] - rad > lb) {
                                lbi = k;
                                lb  = lm[k] - rad;
                            }
                            if (um[k] + rad < ub) {
                                ubi = k;
                                ub  = um[k] + rad;
                            }
                        }
                        if (lb < ub) {
                            aux_m = (S[i+1] - S[j+1]) / double(i-j);
                            if (aux_m > ub)
                                aux_mu = ub;
                            else if (aux_m < lb)
                                aux_mu = lb;
                            else
                                aux_mu = aux_m;
                            aux_c = (aux_mu*aux_mu - 2*aux_mu*aux_m)*(i-j) + (Sq[i+1] - Sq[j+1]) + C[j];
                            if (i < aux_imin)
                                aux_imin = i;
                            if (i > aux_imax)
                                aux_imax = i;
                            if (aux_c < C[i]) {
                                C[i]  = aux_c;
                                mu[i] = aux_mu;
                                L[i]  = j + 2;
                                J[i]  = njmp + 1;
                            }
                        } else {
                            if (lm[lbi] - um[ubi] > eta(lbi+1, ubi+1, n-jmin-lag+1, sd, qm) && j >= imax)
                                isstop = 1;
                            if (lm[lbi] - um[ubi] > eta(lbi+1, ubi+1, i-jmin-lag+1, sd, qm))
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
