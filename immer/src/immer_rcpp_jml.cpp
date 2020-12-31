//// File Name: immer_rcpp_jml.cpp
//// File Version: 0.924


// [[Rcpp::depends(RcppArmadillo)]]

// #include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;
// using namespace arma;

///********************************************************************
///** immer_jml_prob_one_item_one_person
Rcpp::NumericVector immer_jml_prob_one_item_one_person(double theta1, Rcpp::NumericVector b_ii)
{
    int K = b_ii.size();
    Rcpp::NumericVector probs_ii(K);
    double temp=0;
    double tot=0;
    for (int kk=0; kk<K; kk++){
        temp = std::exp( kk * theta1 - b_ii[kk] );
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
///** immer_rcpp_trim_increment
double immer_rcpp_trim_increment(double incr, double max_incr)
{
    double incr1 = incr;
    bool do_cut = TRUE;
    while ( do_cut ){
        do_cut = FALSE;
        if (incr1 > max_incr){
            incr1 = incr1 / 2;
            do_cut = TRUE;
        }
        if (incr1 < - max_incr){
            incr1 = incr1 / 2;
            do_cut = TRUE;
        }
    }
    //-- output
    return incr1;
}
///********************************************************************

///********************************************************************
///** immer_jml_update_item_derivatives
// [[Rcpp::export]]
Rcpp::List immer_jml_update_item_derivatives(Rcpp::NumericVector theta,
        Rcpp::NumericMatrix score_items, int N, int K, int I, Rcpp::IntegerMatrix dat_resp,
        Rcpp::NumericMatrix b, Rcpp::NumericVector A_, Rcpp::NumericVector xsi,
        Rcpp::NumericVector max_incr, Rcpp::NumericMatrix b_fixed,
        Rcpp::NumericVector ItemScore, Rcpp::NumericVector update,
        Rcpp::NumericVector update_weights )
{
    Rcpp::NumericMatrix der1_b(I,K);
    Rcpp::NumericMatrix der2_b(I,K);
    Rcpp::NumericVector b_ii(K+1);
    Rcpp::NumericVector probs_ii(K+1);
    Rcpp::NumericMatrix r(I,K);
    Rcpp::NumericVector rr(I*K*K);

    //        for (ii in 1:I){
    //            sc_ii <- score_items[ii,]
    //            b_ii <- c(0, b[ii,] )
    //            probs_ii <- probs_pcm_one_item(theta=theta, b_ii=b_ii)
    //            der1_b[ii,] <- - sc_ii + colSums(probs_ii * dat_resp[,ii] )[-1]
    //            der2_b[ii,] <- colSums( probs_ii*(1-probs_ii)* dat_resp[,ii] )[-1]
    //        }

    for (int ii=0; ii<I; ii++){
        for (int kk=0;kk<K;kk++){
            b_ii[kk+1] = b(ii,kk);
            r(ii,kk) = 0;
            for (int hh=0;hh<K;hh++){
                rr[ ii + I*kk + I*K*hh ]=0;
            }
        }
        for (int nn=0;nn<N; nn++){
            if ( dat_resp(nn,ii) == 1 ){
                if (update[nn] == 1){
                    probs_ii = immer_jml_prob_one_item_one_person(theta[nn], b_ii);
                    for (int kk=0; kk<K; kk++){
                        r(ii,kk) += update_weights[nn]*probs_ii[kk+1];
                        for (int hh=0;hh<K;hh++){
                            rr[ ii + I*kk + I*K*hh ] += update_weights[nn]*probs_ii[kk+1]*probs_ii[hh+1];
                        }
                    }
                }
            }
        } // end nn
    }  // end ii

    //        der1_xsi <- rep(0,NX)
    //        der2_xsi <- rep(0,NX)
    //        for (kk in 1:K){
    //            der1_xsi <- der1_xsi + colSums( der1_b[,kk] * A[,kk,] )
    //            der2_xsi <- der2_xsi + colSums( der2_b[,kk] * A[,kk,] )
    //        }

    int NX = xsi.size();
    Rcpp::NumericMatrix A_bari(NX,I);
    Rcpp::NumericMatrix AA_bari(NX,I);
    Rcpp::NumericMatrix A_Sq(NX,I);
    Rcpp::NumericVector expected(NX);
    Rcpp::NumericVector err(NX);
    Rcpp::NumericVector err_inv(NX);
    Rcpp::NumericVector scores(NX);
    Rcpp::NumericVector incr(NX);

    //    A_Sq <- AA_bari <- A_bari <- matrix( 0, nrow=NX, ncol=I )
    //    for (kk in 1:K){
    //        A_bari <- A_bari + t( A[, kk, ] * r[, kk ] )
    //        AA_bari <- AA_bari + t( A[, kk, ]^2 * r[, kk ] )
    //    }
    //    for (kk1 in 1:K){
    //        for (kk2 in 1:K){
    //            A_Sq <- A_Sq + t( A[,kk1,] * A[,kk2,] * rr[, kk1, kk2 ] )
    //        }
    //    }
    double Aval=0;
    double Aval2=0;
    double eps=1e-20;

    for (int pp=0; pp<NX; pp++){
        for (int ii=0; ii<I; ii++){
            for (int kk=0; kk<K; kk++){
                Aval = A_[ ii + kk*I + pp*I*K ];
                if ( Aval != 0 ){
                    A_bari(pp,ii) += Aval*r(ii,kk);
                    AA_bari(pp,ii) += Aval*Aval*r(ii,kk);
                    for (int hh=0; hh<K; hh++){
                        Aval2 = A_[ ii + hh*I + pp*I*K ];
                        if (Aval2 != 0){
                            A_Sq(pp,ii) += Aval*Aval2*rr[ii + kk*I + hh*I*K];
                        }
                    }
                }
            }
        }
    }

    //    expected <- rowSums(A_bari, na.rm=TRUE) # sum over items
    //    err <- rowSums(AA_bari - A_Sq, na.rm=TRUE)   #sum over the items
    //    err_inv <- abs(1/( abs(err) + eps ))
    //    scores <- ItemScore - expected
    //    incr <- - err_inv*scores
    //    incr <- immer_trim_increment(incr=incr, max_incr=max_incr)

    for (int pp=0; pp<NX; pp++){
        for (int ii=0; ii<I;ii++){
            expected[pp] += A_bari(pp,ii);
            err[pp] += AA_bari(pp,ii) - A_Sq(pp,ii);
        }
        err_inv[pp] = 1 / ( std::abs(err[pp]) + eps );
        scores[pp] = ItemScore[pp] - expected[pp];
        incr[pp] = - err_inv[pp] * scores[pp];
    }

    Rcpp::NumericVector xsi_new(NX);

    //-- compute increments
    for (int hh=0; hh<NX; hh++){
        incr[hh] = immer_rcpp_trim_increment( incr[hh], max_incr[hh]);
        xsi_new[hh] = xsi[hh] + incr[hh];
    }

    Rcpp::NumericMatrix b_new(I,K);
    for (int ii=0;ii<I;ii++){
        for (int kk=0;kk<K;kk++){
            b_new(ii,kk) = b_fixed(ii,kk);
            for (int hh=0; hh<NX; hh++){
                Aval = A_[ ii + kk*I + hh*I*K ];
                if ( Aval != 0 ){
                    b_new(ii,kk) += Aval * xsi_new[hh];
                }
            }
        }
    }

    //-- output
    return Rcpp::List::create(
                Rcpp::Named("incr") = incr,
                Rcpp::Named("xsi") = xsi_new,
                Rcpp::Named("b") = b_new,
                Rcpp::Named("der2_xsi") = err
            );
}
///********************************************************************

///********************************************************************
///** immer_jml_update_theta_derivatives
// [[Rcpp::export]]
Rcpp::List immer_jml_update_theta_derivatives(Rcpp::NumericVector theta,
        Rcpp::NumericVector score_pers, int N, int K, int I,
        Rcpp::NumericMatrix b, double max_incr, Rcpp::IntegerMatrix dat_resp,
        Rcpp::NumericVector update)
{
    Rcpp::NumericVector theta_new(N);
    int NP=N*I*(K+1);
    Rcpp::NumericVector probs(NP);

    //        der1 <- score_pers
    //        der2 <- rep(0,N)
    //        for (ii in 1:I){
    //            b_ii <- c(0,b[ii,])
    //            probs_ii <- probs_pcm_one_item(theta=theta, b_ii=b_ii)
    //            probs[,ii,] <- probs_ii
    //            M_ii <- rowSums( KM * probs_ii * dat_resp[,ii] )
    //            Var_ii <- rowSums( KM^2 * probs_ii * dat_resp[,ii] ) - M_ii^2
    //            der1 <- der1 - M_ii
    //            der2 <- der2 - Var_ii
    //        }

    Rcpp::NumericVector der1(N);
    Rcpp::NumericVector der2(N);
    Rcpp::NumericVector der2_out(N);
    Rcpp::NumericVector b_ii(K+1);
    Rcpp::NumericVector probs_ii(K+1);
    Rcpp::NumericMatrix bM(I,K+1);

    //-- compute increments
    Rcpp::NumericVector incr(N);
    double eps=1e-20;
    double theta_temp=0;

    for (int ii=0; ii<I; ii++){
        bM(ii,0) = 0;
        for (int kk=0;kk<K;kk++){
            bM(ii,kk+1) = b(ii,kk);
        }
    }

    double M_ii=0;
    double Var_ii=0;
    for (int nn=0; nn<N; nn++){
        der1[nn] = score_pers[nn];
    }

    for (int ii=0; ii<I; ii++){
        b_ii = bM(ii,_);
        for (int nn=0; nn<N; nn++){
            if (dat_resp(nn,ii)==1){
                if (update[nn]==1){
                    probs_ii = immer_jml_prob_one_item_one_person(theta[nn], b_ii);
                    M_ii=0;
                    Var_ii=0;
                    for (int kk=0;kk<K;kk++){
                        M_ii += (kk+1) * probs_ii[kk+1];
                        Var_ii += std::pow(kk+1, 2.0) * probs_ii[kk+1];
                    }
                    Var_ii = Var_ii - std::pow(M_ii, 2.0);
                }
                for (int kk=0; kk<K+1; kk++){
                    probs[nn + ii*N + kk*N*I] = probs_ii[kk];
                }
                der1[nn] += - M_ii;
                der2[nn] += - Var_ii;
            }
        }
    }

    for (int nn=0; nn<N; nn++){
        if (update[nn]==1){
            incr[nn] = der1[nn] / ( std::abs( der2[nn] ) + eps );
            incr[nn] = immer_rcpp_trim_increment( incr[nn], max_incr);
            theta_temp = theta[nn] + incr[nn];
        }
        der2_out[nn] = der2[nn];
        theta_new[nn] = theta_temp;
    }

    //-- output
    return Rcpp::List::create(
                Rcpp::Named("theta") = theta_new,
                Rcpp::Named("der1") = der1,
                Rcpp::Named("der2") = der2_out,
                Rcpp::Named("probs") = probs
            );

}
///********************************************************************
