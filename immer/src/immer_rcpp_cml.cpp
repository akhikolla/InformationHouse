//// File Name: immer_rcpp_cml.cpp
//// File Version: 0.543



// [[Rcpp::depends(RcppArmadillo)]]


#include <Rcpp.h>

using namespace Rcpp;

///********************************************************************
///** immer_cml_splitvec
// [[Rcpp::export]]
Rcpp::List immer_cml_splitvec( Rcpp::NumericVector esf_par1,
                Rcpp::IntegerVector splitvec_len_pp)
{
    int NSPL = splitvec_len_pp.size();
    Rcpp::List esf_par(NSPL);
    int nl=0;
    int uu=0;
    for (int ss=0; ss<NSPL; ss++){
        nl = splitvec_len_pp[ss];
        Rcpp::NumericVector v1(nl);
        for (int hh=0; hh<nl; hh++){
            v1[hh] = esf_par1[uu];
            uu ++;
        }
        esf_par[ss] = v1;
    }
    //-- output
    return esf_par;
}
///********************************************************************


///********************************************************************
///** immer_cml_extract_parmindex
// [[Rcpp::export]]
Rcpp::NumericVector immer_cml_extract_parmindex( Rcpp::NumericMatrix esf_par0,
    Rcpp::IntegerVector parm_index_pp )
{
    int N0 = parm_index_pp.size();
    Rcpp::NumericVector esf_par1(N0);
    for (int ii=0; ii<N0; ii++){
        esf_par1[ii] = esf_par0( parm_index_pp[ii], 0);
    }
    //-- output
    return esf_par1;
}
///********************************************************************


///********************************************************************
///** psychotools_esf
// [[Rcpp::export]]
Rcpp::List psychotools_esf( Rcpp::List esf_par, int order,  bool diff )
{
    // attach function
    Rcpp::Environment pkg = Environment::namespace_env("psychotools");
    Rcpp::Function psychotools_esf_intern = pkg["elementary_symmetric_functions"];
    // apply psychotools function
    Rcpp::List res = psychotools_esf_intern( Rcpp::Named("par") = esf_par,
                        Rcpp::Named("order") = order, Rcpp::Named("diff") = diff );
    //-- output
    return res;
}
///********************************************************************

///********************************************************************
///** immer_cml_cloglik_helper
// [[Rcpp::export]]
double immer_cml_cloglik_helper( Rcpp::NumericMatrix esf_par0,
    Rcpp::List parm_index, Rcpp::List splitvec_len, Rcpp::List suffstat,
    Rcpp::List score_freq, bool diff, int NP )
{
    int order=0;
    double val=0;
    double temp=0;
    double eps = 1e-20;

    for (int pp=0; pp<NP; pp++){

        Rcpp::IntegerVector parm_index_pp = parm_index[pp];
        Rcpp::IntegerVector splitvec_len_pp = splitvec_len[pp];
        Rcpp::NumericVector suffstat_pp = suffstat[pp];
        Rcpp::NumericVector score_freq_pp = score_freq[pp];

        Rcpp::NumericVector esf_par1 = immer_cml_extract_parmindex( esf_par0, parm_index_pp );
        Rcpp::List esf_par = immer_cml_splitvec( esf_par1, splitvec_len_pp);

        // apply psychotools function
        Rcpp::List res0 = psychotools_esf( esf_par, order, diff );

        // res <- - sum(suffstat[[pp]] * esf_par ) - sum(score_freq[[pp]] * log(esf))
        Rcpp::NumericVector esf_res = res0[0];

        int NP1 = esf_par1.size();
        for (int uu=0; uu<NP1; uu++){
            val = val - suffstat_pp[uu] * esf_par1[uu];
        }
        int NE = esf_res.size();
        for (int ee=0; ee<NE; ee++){
            temp = esf_res[ee];
            if (esf_res[ee] < eps){ temp = eps; }
            val = val - score_freq_pp[ee] * log(temp);
        }
    }
    //-- output
    return val;
}
///********************************************************************

///********************************************************************
///** immer_cml_agrad_helper
// [[Rcpp::export]]
Rcpp::NumericVector immer_cml_agrad_helper( Rcpp::NumericMatrix esf_par0,
    Rcpp::List parm_index, Rcpp::List splitvec_len, Rcpp::List suffstat,
    Rcpp::List score_freq, bool diff, int NP, Rcpp::NumericMatrix W,
    Rcpp::LogicalMatrix W_logical, Rcpp::List gr1_list )
{
    int order=1;
    int npar = W.ncol();
    Rcpp::NumericVector agrad(npar);
    double eps = 1e-20;

    for (int pp=0; pp<NP; pp++){

        Rcpp::IntegerVector parm_index_pp = parm_index[pp];
        Rcpp::IntegerVector splitvec_len_pp = splitvec_len[pp];
        Rcpp::NumericVector suffstat_pp = suffstat[pp];
        Rcpp::NumericVector score_freq_pp = score_freq[pp];

        Rcpp::NumericVector esf_par1 = immer_cml_extract_parmindex( esf_par0, parm_index_pp );
        Rcpp::List esf_par = immer_cml_splitvec( esf_par1, splitvec_len_pp);

        // apply psychotools function
        Rcpp::List res0 = psychotools_esf( esf_par, order, diff );

        //    gamma0 <- esf[[1]]
        //    gamma1 <- esf[[2]]
        Rcpp::NumericVector gamma0 = res0[0];
        Rcpp::NumericMatrix gamma1 = res0[1];

        //    W1 <- W[ parm_index[[pp]],, drop=FALSE ]
        int npp = parm_index_pp.size();
        int NW = W.ncol();
        Rcpp::NumericMatrix W1(npp,NW);
        Rcpp::LogicalMatrix W1_logical(npp,NW);
        for (int ii=0; ii<npp; ii++){
            W1(ii,_) = W( parm_index_pp[ii],_);
            W1_logical(ii,_) = W_logical( parm_index_pp[ii],_);
        }
        //    gr1 <- suffstat[[pp]] %*% W1
        Rcpp::NumericVector gr1 = gr1_list[pp];

        //    gr2 <- - colSums( ( score_freq[[pp]] * (gamma1 / gamma0))  %*% W1 )
        int NV = gamma0.size();
        Rcpp::NumericMatrix V(NV,NW);
        Rcpp::NumericVector gr2(NW);
        Rcpp::NumericVector gr(NW);
        double temp=0;
        for (int rr=0; rr<NV; rr++){
            if (gamma0[rr] < eps){
                gamma0[rr] = eps;
            }
            temp = score_freq_pp[rr] / gamma0[rr];
            for (int cc=0; cc<NW; cc++){
                for (int hh=0; hh<npp; hh++){
                    if (W1_logical(hh,cc) ){
                        V(rr,cc) += gamma1(rr,hh) * W1(hh,cc);
                    }
                }
                V(rr,cc) = V(rr,cc) * temp;
            }
        }
        for (int cc=0; cc<NW; cc++){
            for (int rr=0; rr<NV; rr++){
                gr2[cc] += - V(rr,cc);
            }
            gr[cc] = gr1[cc] + gr2[cc];
            agrad[cc] += gr[cc];
        }
    }
    //-- output
    return agrad;
}
///********************************************************************

