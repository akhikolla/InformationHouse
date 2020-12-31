//// File Name: cdm_rcpp_din_deterministic_devcrit.cpp
//// File Version: 3.21


// #include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;


///********************************************************************
///** cdm_rcpp_din_deterministic_devcrit
// [[Rcpp::export]]
Rcpp::List cdm_rcpp_din_deterministic_devcrit( Rcpp::NumericMatrix DAT, Rcpp::NumericMatrix DATRESP,
        Rcpp::NumericMatrix LATRESP, Rcpp::NumericVector GUESS, Rcpp::NumericVector SLIP )
{
    // define row and column numbers
    int N=DAT.nrow();
    int I=DAT.ncol();
    int AP=LATRESP.nrow();

    // create output ghat matrix
    Rcpp::NumericMatrix devcrit (N,AP);
    Rcpp::NumericVector mincrit (N);
    Rcpp::NumericVector rn (1);
    Rcpp::NumericVector indexcrit (N);
    mincrit.fill(10000);

    for (int aa=0;aa<AP;aa++){ // begin attributes
        for (int nn=0;nn<N;nn++){ // begin cases
            for (int ii=0;ii<I;ii++){ // begin item loop
                if ( DATRESP(nn,ii)==1 ){  // begin if datresp == 1
                    if ( (LATRESP(aa,ii) == 1 ) && ( DAT(nn,ii)==0 ) ){
                        devcrit(nn,aa) += SLIP[ii];
                    }
                    if ( (LATRESP(aa,ii) == 0 ) && ( DAT(nn,ii)==1 ) ){
                        devcrit(nn,aa) += GUESS[ii];
                    }
                } // end if datresp == 1
            } // end loop over items ii
            if (mincrit[nn]>devcrit(nn,aa)){
                mincrit[nn]=devcrit(nn,aa);
                indexcrit[nn]=aa+1;
            }
            // handle ties
            if (mincrit[nn]==devcrit(nn,aa)){
                rn(0)=R::runif(0,1);
                if (rn(0)>0.5){
                    mincrit[nn]=devcrit(nn,aa);
                    indexcrit[nn]=aa+1;
                }
            }
        } // end cases
    } // end attributes

    ///////////////////////////////////////
    /// OUTPUT
    return Rcpp::List::create(
                Rcpp::Named("devcrit")=devcrit,
                Rcpp::Named("mincrit")=mincrit,
                Rcpp::Named("indexcrit")=indexcrit
            );
}
///********************************************************************

///********************************************************************
///** cdm_rcpp_din_jml_devcrit
// [[Rcpp::export]]
Rcpp::List cdm_rcpp_din_jml_devcrit( Rcpp::NumericMatrix DAT, Rcpp::NumericMatrix DATRESP,
        Rcpp::NumericMatrix LATRESP, Rcpp::NumericVector GUESS, Rcpp::NumericVector SLIP )
{

    // define row and column numbers
    int N=DAT.nrow();
    int I=DAT.ncol();
    int AP=LATRESP.nrow();

    // create output ghat matrix
    Rcpp::NumericMatrix devcrit (N,AP);
    Rcpp::NumericVector mincrit (N);
    Rcpp::NumericVector indexcrit (N);
    Rcpp::NumericVector rn (1);
    mincrit.fill(-1);
    devcrit.fill(1);

    for (int aa=0;aa<AP;aa++){ // begin attributes
        for (int nn=0;nn<N;nn++){ // begin cases
            for (int ii=0;ii<I;ii++){ // begin item loop
                if ( DATRESP(nn,ii)==1 ){  // begin if datresp == 1
                    if ( (LATRESP(aa,ii) == 1 ) && ( DAT(nn,ii)==0 ) ){
                        devcrit(nn,aa) *= SLIP[ii];
                    }
                    if ( (LATRESP(aa,ii) == 1 ) && ( DAT(nn,ii)==1 ) ){
                        devcrit(nn,aa) *= 1-SLIP[ii];
                    }
                    if ( (LATRESP(aa,ii) == 0 ) && ( DAT(nn,ii)==1 ) ){
                        devcrit(nn,aa) *= GUESS[ii];
                    }
                    if ( (LATRESP(aa,ii) == 0 ) && ( DAT(nn,ii)==0 ) ){
                        devcrit(nn,aa) *= 1-GUESS[ii];
                    }
                } // end if datresp == 1
            } // end loop over items ii
            if (mincrit[nn]<devcrit(nn,aa)){
                mincrit[nn]=devcrit(nn,aa);
                indexcrit[nn]=aa+1;
            }
            // handle ties
            if (mincrit[nn]==devcrit(nn,aa)){
                rn(0)=R::runif(0,1);
                if (rn(0)>0.5){
                    mincrit[nn]=devcrit(nn,aa);
                    indexcrit[nn]=aa+1;
                }
            }
        } // end cases
    } // end attributes

    ///////////////////////////////////////
    /// OUTPUT
    return Rcpp::List::create(
                Rcpp::Named("devcrit")=devcrit,
                Rcpp::Named("mincrit")= mincrit,
                Rcpp::Named("indexcrit")=indexcrit
            );
}
///********************************************************************

