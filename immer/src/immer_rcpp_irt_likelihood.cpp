//// File Name: immer_rcpp_irt_likelihood.cpp
//// File Version: 0.812


// [[Rcpp::depends(RcppArmadillo)]]

// #include <RcppArmadillo.h>
#include <Rcpp.h>
#include <Rmath.h>

using namespace Rcpp;
// using namespace arma;


///********************************************************************
///** immer_gpcm_prob_one_item_one_person
// [[Rcpp::export]]
Rcpp::NumericVector immer_gpcm_prob_one_item_one_person(double theta1,
        Rcpp::NumericVector b_ii, double a)
{
    int K = b_ii.size();
    Rcpp::NumericVector probs_ii(K);
    double temp=0;
    double tot=0;
    for (int kk=0; kk<K; kk++){
        temp = std::exp( a * kk * theta1 - b_ii[kk] );
        tot += temp;
        probs_ii[kk] = temp;
    }
    for (int kk=0; kk<K; kk++){
        probs_ii[kk] = probs_ii[kk] / tot;
    }
    //-- output
    return probs_ii;
}
///********************************************************************

///********************************************************************
///** immer_gpcm_calc_probs
// [[Rcpp::export]]
Rcpp::NumericVector immer_gpcm_calc_probs(Rcpp::NumericVector theta,
    Rcpp::NumericMatrix b, Rcpp::NumericVector a)
{
    int TP = theta.size();
    int K = b.ncol();
    int I = a.size();
    Rcpp::NumericVector probs(I*K*TP);
    Rcpp::NumericVector probs_ii;
    // calculate probabilities
    for (int ii=0; ii<I; ii++){
        for (int tt=0; tt<TP; tt++){
            probs_ii = immer_gpcm_prob_one_item_one_person( theta[tt], b(ii,_), a[ii] );
            for (int kk=0; kk<K; kk++){
                probs[ii + kk*I + tt*I*K] = probs_ii[kk];
            }
        }
    }
    //-- output
    return probs;
}
///********************************************************************

///********************************************************************
///** immer_irt_likelihood_gpcm
// [[Rcpp::export]]
Rcpp::NumericMatrix immer_irt_likelihood_gpcm(Rcpp::NumericVector probs,
    Rcpp::IntegerMatrix dat, Rcpp::IntegerMatrix dat_resp, int TP, int K)
{
    int I = dat.ncol();
    int N = dat.nrow();
    Rcpp::NumericMatrix like(N,TP);
    like.fill(1);
    for (int nn=0; nn<N; nn++){
        for (int tt=0; tt<TP; tt++){
            for (int ii=0; ii<I; ii++){
                if (dat_resp(nn,ii)==1){
                    like(nn,tt) = like(nn,tt) * probs[ ii + dat(nn,ii)*I + tt*I*K];
                }
            }
        }
    }
    //-- output
    return like;
}
///********************************************************************
