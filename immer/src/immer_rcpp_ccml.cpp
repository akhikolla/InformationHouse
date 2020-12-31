//// File Name: immer_rcpp_ccml.cpp
//// File Version: 0.720


// [[Rcpp::depends(RcppArmadillo)]]

// #include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;
// using namespace arma;





///********************************************************************
///** immer_ccml_proc_freq_item_pair
// [[Rcpp::export]]
Rcpp::NumericMatrix immer_ccml_proc_freq_item_pair( Rcpp::IntegerMatrix dat,
    Rcpp::IntegerMatrix dat_resp, int K, Rcpp::NumericVector weights, int ii, int jj )
{
    Rcpp::NumericMatrix freq(K+1,K+1);
    int N = dat.nrow();
    int ii1 = ii - 1;
    int jj1 = jj - 1;
    for (int nn=0; nn<N; nn++){
        if ( ( dat_resp(nn,ii1) == 1 ) & ( dat_resp(nn,jj1) == 1 ) ){
            freq( dat(nn,ii1), dat(nn,jj1) ) += weights[nn];
        }
    }
    //-- output
    return freq;
}
///********************************************************************


///********************************************************************
///** immer_ccml_proc_freq
// [[Rcpp::export]]
Rcpp::NumericMatrix immer_ccml_proc_freq( Rcpp::IntegerMatrix dat,
    Rcpp::IntegerMatrix dat_resp, int K, Rcpp::NumericVector weights)
{

    // ll_index: 0
    // item1: 1
    // item2: 2
    // score: 3
    // cat1: 4
    // cat2: 5
    // n: 6

    int I = dat.ncol();
    int K1 = K+1;
    int N_pair = I * (I-1) / 2;
    int NR = N_pair * K1 * K1;
    Rcpp::NumericMatrix dfr(NR,7);
    Rcpp::NumericMatrix freq_ii_jj;
    int ww=0;
    int lli=1;
    int ss=0;

    for (int ii=1; ii<I+1;ii++){
        for (int jj=ii+1; jj<I+1; jj++){
            freq_ii_jj = immer_ccml_proc_freq_item_pair( dat, dat_resp, K, weights, ii, jj );
            for (int rr=1; rr<2*K; rr++){
                for (int uu=0; uu<rr+1; uu++){
                    ww=rr-uu;
                    if ( (uu<K1) & (ww<K1) ){
                        dfr(ss,0) = lli; // ll_index
                        dfr(ss,1) = ii;  // item 1
                        dfr(ss,2) = jj;  // item 2
                        dfr(ss,3) = rr; // score
                        dfr(ss,4) = uu; // cat1
                        dfr(ss,5) = ww; // cat2
                        dfr(ss,6) = freq_ii_jj(uu,ww); // n
                        ss ++;
                    }
                }
                lli++;
            }
        }
    }

    //-- output
    return dfr;
}
///********************************************************************

///********************************************************************
///** immer_ccml_calc_item_intercepts
// [[Rcpp::export]]
Rcpp::NumericMatrix immer_ccml_calc_item_intercepts(Rcpp::NumericMatrix b_fixed,
        Rcpp::NumericVector A_, Rcpp::NumericVector par)
{
    int I = b_fixed.nrow();
    int K = b_fixed.ncol() - 1;
    int NX = par.size();
    double A_val=0;

    Rcpp::NumericMatrix b(I,K+1);
    for (int kk=0; kk<K; kk++){
        b(_,kk) = b_fixed(_,kk);
    }

    for (int ii=0; ii<I; ii++){
        for (int kk=0; kk<K; kk++){
            for (int pp=0; pp<NX; pp++){
                A_val = A_[ ii + kk*I + pp*I*K ];
                if ( A_val != 0 ){
                    b(ii,kk+1) += A_val * par[pp];
                }
            }
        }
    }

    //-- output
    return b;
}
///********************************************************************


///********************************************************************
///** immer_ccml_probs
// [[Rcpp::export]]
Rcpp::NumericVector immer_ccml_probs(Rcpp::NumericMatrix b, Rcpp::NumericVector ll_index1,
    Rcpp::NumericVector item10, Rcpp::NumericVector item20, Rcpp::NumericVector cat1,
    Rcpp::NumericVector cat2, int max_ll_index )
{
    int NL = ll_index1.size();
    Rcpp::NumericVector zaehl(NL);
    Rcpp::NumericVector sz(max_ll_index);
    Rcpp::NumericVector opt_fct(NL);
    // dfr$zaehl <- exp( - b[ cbind( dfr$item1, dfr$cat1 + 1) ] - b[ cbind( dfr$item2, dfr$cat2 + 1) ] )
    // sz <- rowsum( dfr$zaehl, dfr$ll_index )
    // dfr$nenn <- sz[dfr$ll_index]
    // dfr$opt_fct <- dfr$zaehl / dfr$nenn
    // compute denominator
    for (int ll=0; ll<NL; ll++){
        zaehl[ll] = std::exp( - b( item10[ll], cat1[ll] ) - b( item20[ll], cat2[ll] ) );
        sz[ ll_index1[ll] ] += zaehl[ll];
    }
    // compute conditional probabilities
    for (int ll=0; ll<NL; ll++){
        opt_fct[ll] = zaehl[ll] / sz[ ll_index1[ll] ];
    }
    //-- output
    return opt_fct;
}
///********************************************************************

///********************************************************************
///** immer_ccml_probs
// [[Rcpp::export]]
Rcpp::NumericVector immer_ccml_probs_from_par(Rcpp::NumericMatrix b_fixed,
    Rcpp::NumericVector A_, Rcpp::NumericVector par, Rcpp::NumericVector ll_index1,
    Rcpp::NumericVector item10, Rcpp::NumericVector item20, Rcpp::NumericVector cat1,
    Rcpp::NumericVector cat2, int max_ll_index, int pp1, int pp2, double h1, double h2 )
{
    int NX = par.size();
    Rcpp::NumericVector par1(NX);
    for (int pp=0; pp<NX; pp++){
        par1[pp] = par[pp];
    }
    if (pp1 >= 0){
        par1[pp1] = par[pp1] + h1;
    }
    if (pp2 >= 0){
        par1[pp2] = par[pp2] + h2;
    }
    //-- compute item intercepts
    Rcpp::NumericMatrix b = immer_ccml_calc_item_intercepts(b_fixed, A_, par1);
    //-- compute conditional probabilities
    Rcpp::NumericVector opt_fct = immer_ccml_probs(b, ll_index1, item10, item20, cat1, cat2, max_ll_index );
    //-- output
    return opt_fct;
}
///********************************************************************

///********************************************************************
///** log_eps
double log_eps(double x, double eps)
{
    // double val = x;
    // if (val < eps){ val = eps; }
    double val=x+eps;
    val = std::log(val);
    return val;
}
///********************************************************************

///********************************************************************
///** immer_ccml_opt_function
// [[Rcpp::export]]
double immer_ccml_opt_function(Rcpp::NumericMatrix b, Rcpp::NumericVector ll_index1,
    Rcpp::NumericVector item10, Rcpp::NumericVector item20, Rcpp::NumericVector cat1,
    Rcpp::NumericVector cat2, Rcpp::NumericVector n, Rcpp::NumericVector ntot,
    int max_ll_index )
{
    // int I = b.nrow();
    // int K = b.ncol();
    int NL = ll_index1.size();
    double val=0;
    double eps = 1e-50;
    // dfr$log_opt_fct <- log( dfr$opt_fct )
    // sum( dfr$log_opt_fct * dfr$n )

    // compute conditional probabilities
    Rcpp::NumericVector opt_fct = immer_ccml_probs(b, ll_index1, item10, item20, cat1, cat2, max_ll_index );

    // compute optimization function
    for (int ll=0; ll<NL; ll++){
        val += n[ll] * log_eps( opt_fct[ll], eps );
    }
    val = - val;

    //-- output
    return val;
}
///********************************************************************

///********************************************************************
///** immer_ccml_opt_function_par
// [[Rcpp::export]]
double immer_ccml_opt_function_par( Rcpp::NumericMatrix b_fixed,
    Rcpp::NumericVector A_, Rcpp::NumericVector par, Rcpp::NumericVector ll_index1,
    Rcpp::NumericVector item10, Rcpp::NumericVector item20, Rcpp::NumericVector cat1,
    Rcpp::NumericVector cat2, Rcpp::NumericVector n, Rcpp::NumericVector ntot,
    int max_ll_index )
{
    //-- compute item intercepts
    Rcpp::NumericMatrix b = immer_ccml_calc_item_intercepts(b_fixed, A_, par);
    //-- compute optimization function
    double val = immer_ccml_opt_function(b, ll_index1, item10, item20, cat1,
                    cat2, n, ntot, max_ll_index );
    //-- output
    return val;
}
///********************************************************************

///********************************************************************
///** immer_ccml_gradient
// [[Rcpp::export]]
Rcpp::NumericMatrix immer_ccml_gradient(Rcpp::NumericMatrix b, Rcpp::NumericVector ll_index1,
    Rcpp::NumericVector item10, Rcpp::NumericVector item20, Rcpp::NumericVector cat1,
    Rcpp::NumericVector cat2, Rcpp::NumericVector n, Rcpp::NumericVector ntot,
    int max_ll_index )
{
    int I = b.nrow();
    int K = b.ncol();
    int NL = ll_index1.size();
    Rcpp::NumericMatrix grad(I, K-1);
    grad.fill(0);

    // dfr$opt_fct <- dfr$zaehl / dfr$nenn
    // ind <- c( which( ( dfr$item1 == ii0 ) & ( dfr$cat1 == cat0 ) ), which( ( dfr$item2 == ii0 ) &
    //                ( dfr$cat2 == cat0 ) ) )
    // - sum( dfr$n[ind] - dfr$ntot[ind] * dfr$opt_fct[ind] )

    // compute conditional probabilities
    Rcpp::NumericVector opt_fct = immer_ccml_probs(b, ll_index1, item10, item20, cat1, cat2, max_ll_index );
    // compute derivatives
    for (int ll=0; ll<NL; ll++){
        if (cat1[ll]>0){
            grad( item10[ll], cat1[ll]-1 ) += n[ll] - ntot[ll] * opt_fct[ll];
        }
        if (cat2[ll]>0){
            grad( item20[ll], cat2[ll]-1 ) += n[ll] - ntot[ll] * opt_fct[ll];
        }
    }

    //-- output
    return grad;
}
///********************************************************************

///********************************************************************
///** immer_ccml_gradient_par
// [[Rcpp::export]]
Rcpp::NumericVector immer_ccml_gradient_par( Rcpp::NumericMatrix b_fixed,
    Rcpp::NumericVector A_, Rcpp::NumericVector par, Rcpp::NumericVector ll_index1,
    Rcpp::NumericVector item10, Rcpp::NumericVector item20, Rcpp::NumericVector cat1,
    Rcpp::NumericVector cat2, Rcpp::NumericVector n, Rcpp::NumericVector ntot,
    int max_ll_index )
{
    //-- compute item intercepts
    Rcpp::NumericMatrix b = immer_ccml_calc_item_intercepts(b_fixed, A_, par);
    //-- compute gradient with respect zo b
    Rcpp::NumericMatrix grad = immer_ccml_gradient(b, ll_index1, item10, item20, cat1,
                    cat2, n, ntot, max_ll_index );

    // compute derivative with respect to par
    int NX = par.size();
    int K = b_fixed.ncol()-1;
    int I = b_fixed.nrow();
    Rcpp::NumericVector grad_par(NX);
    double Aval=0;

    // d grad / d par = d grad / d b * d b / d par

    for (int pp=0; pp<NX; pp++){
        for (int ii=0; ii<I; ii++){
            for (int kk=0; kk<K; kk++){
                Aval = A_[ ii + kk*I + pp*I*K ];
                if (Aval != 0 ){
                    grad_par[pp] += Aval * grad(ii,kk);
                }
            }
        }
    }

    //-- output
    return grad_par;
}
///********************************************************************


///********************************************************************
// cross-product information
///** immer_ccml_se
// [[Rcpp::export]]
Rcpp::List immer_ccml_se( Rcpp::NumericMatrix b_fixed,
    Rcpp::NumericVector A_, Rcpp::NumericVector par, Rcpp::NumericVector ll_index1,
    Rcpp::NumericVector item10, Rcpp::NumericVector item20, Rcpp::NumericVector cat1,
    Rcpp::NumericVector cat2, Rcpp::NumericVector n, Rcpp::NumericVector ntot,
    int max_ll_index, double h )
{
    //-- compute item intercepts
    Rcpp::NumericMatrix b = immer_ccml_calc_item_intercepts(b_fixed, A_, par);
    int NX = par.size();
    int NL = ll_index1.size();
    Rcpp::NumericMatrix xpd_mat(NX,NX);
    Rcpp::NumericMatrix obs_mat(NX,NX);
    Rcpp::NumericMatrix der1_mat(NL,NX);
    double eps=1e-50;

    int N=0;
    for (int ll=0; ll<NL; ll++){
        N += n[ll];
    }

    //-----------------------------------------------
    //-- compute cross-product information
    Rcpp::NumericVector opt_p, opt_m, opt_0;
    Rcpp::NumericMatrix opt_p0(NL,NX);
    Rcpp::NumericMatrix opt_m0(NL,NX);
    Rcpp::NumericVector par_temp(NX);
    for (int pp=0; pp<NX; pp++){
        par_temp[pp] = par[pp];
    }

    for (int pp=0; pp<NX; pp++){
        // f(x+h)
        opt_p = immer_ccml_probs_from_par(b_fixed, A_, par, ll_index1, item10, item20, cat1, cat2,
                            max_ll_index, pp, -1, h, 0 );
        opt_p0(_,pp)=opt_p;
        // f(x-h)
        opt_m = immer_ccml_probs_from_par(b_fixed, A_, par, ll_index1, item10, item20, cat1, cat2,
                            max_ll_index, pp, -1, -h, 0 );
        opt_m0(_,pp)=opt_m;
        for (int ll=0; ll<NL; ll++){
            der1_mat(ll,pp) = - ( ( log_eps( opt_p[ll], eps ) - log_eps( opt_m[ll], eps ) ) ) / (2*h);
        }
        par_temp[pp] = par[pp];
    }

    for (int pp1=0; pp1<NX; pp1++){
        for (int pp2=pp1; pp2<NX; pp2++){
            for (int ll=0; ll<NL; ll++){
                xpd_mat(pp1,pp2) += n[ll] * der1_mat(ll,pp1) * der1_mat(ll,pp2);
            }
            xpd_mat(pp1,pp2) = xpd_mat(pp1,pp2) / N;
            if (pp1 < pp2){
                xpd_mat(pp2,pp1) = xpd_mat(pp1,pp2);
            }
        }
    }

    //-----------------------------------------------
    //-- compute observed information
    b = immer_ccml_calc_item_intercepts(b_fixed, A_, par);
    opt_0 = immer_ccml_probs(b, ll_index1, item10, item20, cat1, cat2, max_ll_index );
    double val=0;
    // observed information: diagonal elements
    for (int pp=0; pp<NX; pp++){
        for (int ll=0; ll<NL; ll++){
            val = log_eps( opt_p0(ll,pp), eps ) - 2*log_eps( opt_0[ll], eps ) + log_eps( opt_m0(ll,pp), eps );
            obs_mat(pp,pp) += - n[ll] * val / std::pow(h, 2.0);
        }
        obs_mat(pp,pp) = obs_mat(pp,pp) / N;
    }
    // observed information non-diagonal elements
    Rcpp::NumericVector opt_11;
    Rcpp::NumericVector opt_10;
    Rcpp::NumericVector opt_01;
    Rcpp::NumericVector opt_00;

    for (int pp1=0; pp1<NX-1; pp1++){
        for (int pp2=pp1+1; pp2<NX; pp2++){

            // f(x+h, y+h)
            opt_11 = immer_ccml_probs_from_par(b_fixed, A_, par, ll_index1, item10, item20, cat1, cat2,
                                    max_ll_index, pp1, pp2, h, h );
            // f(x+h, y-h)
            opt_10 = immer_ccml_probs_from_par(b_fixed, A_, par, ll_index1, item10, item20, cat1, cat2,
                                    max_ll_index, pp1, pp2, h, -h );
            // f(x-h, y+h)
            opt_01 = immer_ccml_probs_from_par(b_fixed, A_, par, ll_index1, item10, item20, cat1, cat2,
                                    max_ll_index, pp1, pp2, -h, h );
            // f(x-h, y-h)
            opt_00 = immer_ccml_probs_from_par(b_fixed, A_, par, ll_index1, item10, item20, cat1, cat2,
                                    max_ll_index, pp1, pp2, -h, -h );
            // f1,1 - f1,-1 - f-1,1 + f-1,-1 / 4*h^2
            for (int ll=0; ll<NL; ll++){
                val = log_eps( opt_11[ll], eps) - log_eps( opt_10[ll], eps) - log_eps( opt_01[ll], eps) +
                            log_eps( opt_00[ll], eps);
                obs_mat(pp1,pp2) += - n[ll] * val / ( 4 * h * h );
            }
            obs_mat(pp1,pp2) = obs_mat(pp1,pp2) / N;
            obs_mat(pp2,pp1) = obs_mat(pp1,pp2);
        }
    }

    //-- output
    return Rcpp::List::create(
            Rcpp::Named("xpd_mat") = xpd_mat,
            Rcpp::Named("obs_mat") = obs_mat,
            Rcpp::Named("der1_mat") = der1_mat,
            Rcpp::Named("N") = N
        );
}
///********************************************************************

