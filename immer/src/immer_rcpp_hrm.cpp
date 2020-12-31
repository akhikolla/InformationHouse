//// File Name: immer_rcpp_hrm.cpp
//// File Version: 0.464


// [[Rcpp::depends(RcppArmadillo)]]
// [[RcppNOplugins(unwindProtect)]]

// includes from the plugin
#include <RcppArmadillo.h>
// #include <Rcpp.h>


using namespace Rcpp;
// using namespace arma;





//**********************************************************************
// sample from multinomial distribution
// subimmer_sample_prob_index
// [[Rcpp::export]]
Rcpp::NumericVector subimmer_sample_prob_index( Rcpp::NumericMatrix probs, Rcpp::NumericVector rn )
{
    int N = rn.size();
    int K = probs.ncol();
    Rcpp::NumericVector xi(N);
    double t1=0;
    for (int nn=0;nn<N;nn++){
        t1=0;
        for (int kk=0;kk<K;kk++){
            t1 = t1 + probs(nn,kk);
            if ( rn[nn] < t1 ){
                xi[nn] = kk;
                break;
            }
        }
    }
    return xi;
}
//**********************************************************************

//**********************************************************************
//*********** probabilities GPCM ***************************************
//
// subimmer_probs_gpcm_rcpp
// [[Rcpp::export]]
Rcpp::NumericVector subimmer_probs_gpcm_rcpp( Rcpp::NumericVector x,
    Rcpp::NumericVector theta,  Rcpp::NumericVector b,
    Rcpp::NumericVector a, int K, Rcpp::NumericVector x_ind )
{

    int N = x.size();
    Rcpp::NumericVector l1(K+1);
    //    N <- length(theta)
    //    KM <- matrix( 0:K, nrow = N, ncol= K+1, byrow=TRUE)
    //    b0 <- c( 0, b[1:K] )
    Rcpp::NumericVector b0(K+1);
    b0[0] = 0;
    for (int kk=1;kk<K+1;kk++){
        b0[kk] = b[kk-1];
    }

    //    bM <- matrix( b0, nrow = N, ncol= K+1, byrow=TRUE)
    //    probs <- exp( a * KM *  theta - bM )
    //    probs <- probs / rowSums(probs, na.rm=TRUE)
    Rcpp::NumericVector probs(N);
    double t1 = 0;
    for (int nn=0;nn<N;nn++){
        if ( x_ind[nn] > 0 ){
            t1 = 0;
            for (int kk=0; kk<K+1;kk++){
                l1[kk] = std::exp( a[0] * kk * theta[nn] - b0[kk] );
                t1 = t1 + l1[kk];
            }
            probs[nn] = l1[ x[nn]  ] / t1;
        } else {
            probs[nn] = 1;
        }
    }
    // OUTPUT
    return probs;
}
//**********************************************************************

//**********************************************************************
//****** probabilities HRM ************************************
//
// subimmer_probs_hrm_rcpp
// [[Rcpp::export]]
Rcpp::NumericVector subimmer_probs_hrm_rcpp( Rcpp::NumericVector x,
    Rcpp::NumericVector xi, Rcpp::NumericVector phi,
    Rcpp::NumericVector psi, int K, Rcpp::NumericVector x_ind )
{
    int N = x.size();
    Rcpp::NumericVector l1(K+1);
    Rcpp::NumericVector probs(N);
    double t1 = 0;
    for (int nn=0; nn < N; nn++){
        //    KM <- matrix( 0:K, nrow=N, ncol=K+1, byrow=TRUE )
        //    p1 <- exp( - ( KM - ( xi + phi ) )^2 / ( 2 * psi ) )
        //    probs <- p1 / rowSums(p1, na.rm=TRUE)
        if ( x_ind[nn] > 0 ){
            t1 = 0;
            for (int kk=0;kk<K+1;kk++){
                l1[kk] = std::exp( - std::pow( kk - xi[nn] - phi[nn], 2.0 ) / 2 / psi[nn] );
                t1 += l1[kk];
            }
            probs[nn] = l1[ x[nn] ] / t1;
        } else {
            probs[nn] = 1;
        }
    }
    return probs;
}
//**********************************************************************

//**********************************************************************
//*********** probabilities GPCM testlet *******************************
//
// subimmer_probs_gpcm_testlet_rcpp
// [[Rcpp::export]]
Rcpp::NumericVector subimmer_probs_gpcm_testlet_rcpp( Rcpp::NumericVector x,
    Rcpp::NumericVector theta,  Rcpp::NumericVector u, Rcpp::NumericVector b,
    Rcpp::NumericVector a, int K, Rcpp::NumericVector x_ind )
{
    int N = x.size();
    Rcpp::NumericVector l1(K+1);

    //    N <- length(theta)
    //    KM <- matrix( 0:K, nrow = N, ncol= K+1, byrow=TRUE)
    //    b0 <- c( 0, b[1:K] )
    Rcpp::NumericVector b0(K+1);
    b0[0] = 0;
    for (int kk=1;kk<K+1;kk++){
        b0[kk] = b[kk-1];
    }

    //    bM <- matrix( b0, nrow = N, ncol= K+1, byrow=TRUE)
    //    probs <- exp( a * KM *  theta - bM )
    //    probs <- probs / rowSums(probs, na.rm=TRUE)
    Rcpp::NumericVector probs(N);

    double t1 = 0;
    for (int nn=0;nn<N;nn++){
        if ( x_ind[nn] > 0 ){
            t1 = 0;
            for (int kk=0; kk<K+1;kk++){
                l1[kk] = std::exp( a[0] * kk * theta[nn] + kk*u[nn] - b0[kk] );
                t1 = t1 + l1[kk];
            }
            probs[nn] = l1[ x[nn]  ] / t1;
        } else {
            probs[nn] = 1;
        }
    }
    // OUTPUT
    return probs;
}
//**********************************************************************

//**********************************************************************
/// immer_sampling_xi
// [[Rcpp::export]]
Rcpp::NumericVector immer_sampling_xi( Rcpp::NumericVector x,
    Rcpp::NumericVector theta, Rcpp::NumericVector b, Rcpp::NumericVector a,
    int K, Rcpp::NumericVector x_ind, Rcpp::NumericVector phi,
    Rcpp::NumericVector psi, double eps, Rcpp::NumericVector pid,
    Rcpp::NumericVector rater, int N )
{
    Rcpp::NumericMatrix probs(N,K+1);
    int ND = x.size();
    // vector for cumulated probabilities
    Rcpp::NumericVector p3(N);
    Rcpp::NumericVector x_kk(N);
    // vector of phi for computation of HRM probabilities
    Rcpp::NumericVector phi_rater(ND);
    Rcpp::NumericVector psi_rater(ND);
    Rcpp::NumericVector xi_kk(ND);

    for( int kk=0; kk < K+1; kk++){ // category kk
        for (int nn=0;nn<N;nn++){
            x_kk[nn] = kk;
        }
        // probabilities GPCM
        Rcpp::NumericVector p1 = subimmer_probs_gpcm_rcpp( x_kk, theta, b, a, K, x_ind );
        // probabilities HRM
        for (int nn=0;nn<ND;nn++){
            phi_rater[nn] = phi[ rater[nn] - 1 ];
            psi_rater[nn] = psi[ rater[nn] - 1 ];
            xi_kk[nn] = kk;
        }
        Rcpp::NumericVector p2 = subimmer_probs_hrm_rcpp( x, xi_kk, phi_rater, psi_rater, K, x_ind );
        for (int nn=0;nn<N;nn++){
            p3[nn] = 0;
        }
        for (int nn=0;nn<ND;nn++){
            p2[nn] = std::log( p2[nn] + eps );
            p3[ pid[nn]-1 ] += p2[nn];
        }
        for (int nn=0; nn<N; nn++){
            probs(nn,kk) = std::exp( p3[nn] ) * p1[nn];
        }
    } // end category kk

    // standardize probabilities
    double t1=0;
    for (int nn=0;nn<N;nn++){
        t1 = 0;
        for (int kk=0;kk<K+1;kk++){
            t1 += probs(nn,kk);
        }
        for (int kk=0;kk<K+1;kk++){
            probs(nn,kk) = probs(nn,kk) / t1;
        }
    }

    // draw a random number
    Rcpp::NumericVector rn(N);
    for (int nn=0;nn<N;nn++){
        rn[nn] = ::Rf_runif( 0.0, 1.0 );
    }

    // draw xi vector
    Rcpp::NumericVector xi_samp = subimmer_sample_prob_index( probs, rn);

    // OUTPUT
    return xi_samp;
}
//**********************************************************************

//**********************************************************************
/// probs_gpcm_rcpp
// [[Rcpp::export]]
Rcpp::NumericVector probs_gpcm_rcpp( Rcpp::NumericVector x,
        Rcpp::NumericVector theta, Rcpp::NumericVector b, Rcpp::NumericVector a,
        int K, Rcpp::NumericVector x_ind )
{
    int N = x.size();
    Rcpp::NumericVector l1(K+1);
    Rcpp::NumericVector b0(K+1);
    b0[0] = 0;
    for (int kk=1;kk<K+1;kk++){
        b0[kk] = b[kk-1];
    }
    Rcpp::NumericVector probs(N);
    double t1=0;
    for (int nn=0;nn<N;nn++){
        if ( x_ind[nn] > 0 ){
            t1 = 0;
            for (int kk=0; kk<K+1;kk++){
                l1[kk] = std::exp( a[0] * kk * theta[nn] - b0[kk] );
                t1 = t1 + l1[kk];
            }
            probs[nn] = l1[ x[nn] ] / t1;
        } else {
                probs[nn] = 1;
        }
    }
    return probs;
}
//**********************************************************************

//**********************************************************************
/// probs_hrm_rcpp
// [[Rcpp::export]]
Rcpp::NumericVector probs_hrm_rcpp( Rcpp::NumericVector x, Rcpp::NumericVector xi,
        Rcpp::NumericVector phi, Rcpp::NumericVector psi, int K,  Rcpp::NumericVector x_ind )
{
    int N = x.size();
    Rcpp::NumericVector l1(K+1);
    Rcpp::NumericVector probs(N);
    double t1 = 0;
    for (int nn=0; nn < N; nn++){
        if ( x_ind[nn] > 0 ){
            t1 = 0;
            for (int kk=0;kk<K+1;kk++){
                l1[kk] = std::exp( - std::pow( kk - xi[nn] - phi[nn], 2.0 ) / 2 / psi[nn] );
                t1 += l1[kk];
            }
            probs[nn] = l1[ x[nn] ] / t1;
        } else {
            probs[nn] = 1;
        }
    }
    // OUTPUT
    return probs;
}
//**********************************************************************

//**********************************************************************
/// sample_prob_index
// [[Rcpp::export]]
Rcpp::NumericVector sample_prob_index( Rcpp::NumericMatrix probs, Rcpp::NumericVector rn)
{
    int N = rn.size();
    int K = probs.ncol();
    Rcpp::NumericVector xi(N);
    double t1=0;
    for (int nn=0;nn<N;nn++){
        t1=0;
        for (int kk=0;kk<K;kk++){
            t1 = t1 + probs(nn,kk);
            if ( rn[nn] < t1 ){
                xi[nn] = kk;
                break;
            }
        }
    }
    // OUTPUT
    return xi;
}
//**********************************************************************

//**********************************************************************
/// probs_gpcm_testlet_rcpp
// [[Rcpp::export]]
Rcpp::NumericVector probs_gpcm_testlet_rcpp( Rcpp::NumericVector x,
        Rcpp::NumericVector theta, Rcpp::NumericVector u,
        Rcpp::NumericVector b, Rcpp::NumericVector a, int K,
        Rcpp::NumericVector x_ind )
{
    int N = x.size();
    Rcpp::NumericVector l1(K+1);
    Rcpp::NumericVector b0(K+1);
    b0[0] = 0;
    for (int kk=1;kk<K+1;kk++){
        b0[kk] = b[kk-1];
    }

    Rcpp::NumericVector probs(N);
    double t1 = 0;
    for (int nn=0;nn<N;nn++){
        if ( x_ind[nn] > 0 ){
            t1 = 0;
            for (int kk=0; kk<K+1;kk++){
                l1[kk] = std::exp( a[0] * kk * theta[nn] + kk*u[nn] - b0[kk] );
                t1 = t1 + l1[kk];
            }
            probs[nn] = l1[ x[nn]  ] / t1;
        } else {
            probs[nn] = 1;
        }
    }
    // OUTPUT
    return probs;
}
//**********************************************************************
