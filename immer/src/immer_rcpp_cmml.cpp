//// File Name: immer_rcpp_cmml.cpp
//// File Version: 0.8904



// [[Rcpp::depends(RcppArmadillo)]]

// #include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;
// using namespace arma;

// [include_header_file]
#include "immer_rcpp_dnorm_pbivnorm_drezner.h"



///********************************************************************
///** immer_cmml_calc_probs
// [[Rcpp::export]]
Rcpp::NumericVector immer_cmml_calc_probs(Rcpp::NumericMatrix tau,
        Rcpp::NumericMatrix rho, Rcpp::IntegerMatrix item_table)
{
    int IP = item_table.nrow();
    Rcpp::NumericVector probs(IP);
    int item1, item2, cat1, cat2, cat1p, cat2p;
    double F11, F10, F01, F00, rho_temp, tau1a, tau0a, tau1b, tau0b, val;

    for (int pp=0; pp<IP; pp++){
        item1 = item_table(pp,0);
        item2 = item_table(pp,1);
        cat1 = item_table(pp,2);
        cat2 = item_table(pp,3);
        cat1p = cat1 + 1;
        cat2p = cat2 + 1;
        rho_temp = rho( item1, item2);
        tau1a = tau( item1, cat1p);
        tau0a = tau( item1, cat1);
        tau1b = tau( item2, cat2p);
        tau0b = tau( item2, cat2);
        F11 = pbivnorm_drezner_numeric_arguments( tau1a, tau1b, rho_temp );
        F10 = pbivnorm_drezner_numeric_arguments( tau1a, tau0b, rho_temp );
        F01 = pbivnorm_drezner_numeric_arguments( tau0a, tau1b, rho_temp );
        F00 = pbivnorm_drezner_numeric_arguments( tau0a, tau0b, rho_temp );
        val = F11 - F10 - F01 + F00;
        if (val<0){ val=0.0;}
        probs[pp] = val;
    }
    return probs;
}
///********************************************************************

///********************************************************************
///** immer_cmml_probs_derivatives_tau_rho
// item_table ... item i, item j, category h, category k
// calculation of derivatives of probabilities p_ijhk with respect to
// tau and rho
// [[Rcpp::export]]
Rcpp::List immer_cmml_probs_derivatives_tau_rho(Rcpp::NumericMatrix tau,
        Rcpp::NumericMatrix rho, Rcpp::IntegerMatrix item_table,
        Rcpp::IntegerMatrix rho_index, int N_rho,
        Rcpp::IntegerMatrix tau_index, int N_tau )
{
    int IP = item_table.nrow();
    Rcpp::NumericMatrix rho_der(IP, N_rho);
    Rcpp::NumericMatrix tau_der(IP, N_tau);
    rho_der.fill(0);
    tau_der.fill(0);
    int item1, item2, cat1, cat2, cat1p, cat2p;
    double F11, F10, F01, F00, rho_temp, tau1a, tau0a, tau1b, tau0b, val1, val2;

    for (int pp=0; pp<IP; pp++){
        item1 = item_table(pp,0);
        item2 = item_table(pp,1);
        cat1 = item_table(pp,2);
        cat2 = item_table(pp,3);
        cat1p = cat1 + 1;
        cat2p = cat2 + 1;
        rho_temp = rho( item1, item2);
        tau1a = tau( item1, cat1p);
        tau0a = tau( item1, cat1);
        tau1b = tau( item2, cat2p);
        tau0b = tau( item2, cat2);

        // derivative with respect to rho: dF / d(rho)
        F11 = pbivnorm_drezner_derivative_rho_numeric( tau1a, tau1b, rho_temp );
        F10 = pbivnorm_drezner_derivative_rho_numeric( tau1a, tau0b, rho_temp );
        F01 = pbivnorm_drezner_derivative_rho_numeric( tau0a, tau1b, rho_temp );
        F00 = pbivnorm_drezner_derivative_rho_numeric( tau0a, tau0b, rho_temp );

        val1 = F11 - F10 - F01 + F00;
        rho_der(pp, rho_index(item1, item2) ) = val1;

        // derivative with respect to tau: dF/d(tau) first item
        F11 = pbivnorm_drezner_derivative_a_numeric( tau1a, tau1b, rho_temp );
        F10 = pbivnorm_drezner_derivative_a_numeric( tau1a, tau0b, rho_temp );
        F01 = pbivnorm_drezner_derivative_a_numeric( tau0a, tau1b, rho_temp );
        F00 = pbivnorm_drezner_derivative_a_numeric( tau0a, tau0b, rho_temp );

        val1 = F11 - F10;
        val2 = -F01 + F00;
        tau_der(pp, tau_index(item1, cat1p) ) = val1;
        tau_der(pp, tau_index(item1, cat1) ) = val2;

        // derivative with respect to tau: dF/d(tau) second item
        F11 = pbivnorm_drezner_derivative_b_numeric( tau1a, tau1b, rho_temp );
        F10 = pbivnorm_drezner_derivative_b_numeric( tau1a, tau0b, rho_temp );
        F01 = pbivnorm_drezner_derivative_b_numeric( tau0a, tau1b, rho_temp );
        F00 = pbivnorm_drezner_derivative_b_numeric( tau0a, tau0b, rho_temp );

        val1 = F11 - F01;
        val2 = -F10 + F00;
        tau_der(pp, tau_index(item2, cat2p) ) = val1;
        tau_der(pp, tau_index(item2, cat2) ) = val2;
    }

    return Rcpp::List::create(
                Rcpp::Named("probs_der_rho") = rho_der,
                Rcpp::Named("probs_der_tau") = tau_der
        );

}
///********************************************************************


///********************************************************************
///** immer_bilinear_form
///  t(x) %*% A %*% y
double immer_bilinear_form(Rcpp::NumericVector x, Rcpp::NumericMatrix A,
    Rcpp::NumericVector y )
{
    int D=A.ncol();
    double res=0.0;
    for (int dd=0; dd<D; dd++){
        for (int ee=0; ee<D; ee++){
            res += x[dd]*y[ee]*A(dd,ee);
        }
    }
    return res;
}
///********************************************************************

///********************************************************************
///** immer_derivative_correlation
// [[Rcpp::export]]
double immer_derivative_correlation( double cov12, double var1, double var2,
    double cov12_der, double var1_der, double var2_der )
{
    double val=0;
    double sd1 = std::sqrt(var1);
    double sd2 = std::sqrt(var2);
    double temp1 = cov12_der / sd1 / sd2;
    val = val + temp1;
    temp1 = -0.5 * var1_der * cov12 / std::pow( sd1, 3.0) / sd2;
    val = val + temp1;
    temp1 = -0.5 * var2_der * cov12 / std::pow( sd2, 3.0) / sd1;
    val = val + temp1;
    return val;
}
///********************************************************************

///********************************************************************
///** immer_cmml_trafo_irt_parameters
// [[Rcpp::export]]
Rcpp::List immer_cmml_trafo_irt_parameters(Rcpp::NumericMatrix LAM,
        Rcpp::NumericMatrix PHI, Rcpp::NumericMatrix PSI, Rcpp::NumericMatrix GAM,
        Rcpp::IntegerMatrix item_pairs, Rcpp::IntegerMatrix LAM_index,
        int N_LAM, Rcpp::IntegerMatrix PHI_index, int N_PHI, Rcpp::IntegerMatrix PSI_index,
        int N_PSI, Rcpp::IntegerMatrix tau_index, int N_tau, Rcpp::IntegerMatrix GAM_index,
        int N_GAM, Rcpp::IntegerMatrix rho_index, int N_rho )
{
    int D = PHI.ncol();
    int I = LAM.nrow();
    int IP = item_pairs.nrow();
    int K = GAM.ncol();

    Rcpp::NumericVector var_item(I);
    Rcpp::NumericVector sd_item(I);
    Rcpp::NumericVector cov_item_pair(IP);
    Rcpp::NumericVector cor_item_pair(IP);
    Rcpp::NumericMatrix rho(I,I);
    Rcpp::NumericMatrix tau(I,K);
    Rcpp::NumericVector x(D);
    Rcpp::NumericVector y(D);
    int item1, item2;

    // item-wise variances and standard deviations
    for (int ii=0; ii<I; ii++){
        x = LAM(ii,_);
        var_item[ii] = immer_bilinear_form( x, PHI, x) + 1;
        sd_item[ii] = std::sqrt( var_item[ii] );
    }

    // item-wise thresholds
    for (int ii=0; ii<I; ii++){
        tau(ii,0) = -999;
        for (int hh=1; hh<K-1; hh++){
            tau(ii,hh) = GAM(ii,hh) / sd_item[ii];
        }
        tau(ii,K-1)=999;
    }

    // item-pair-wise covariances and correlations
    for (int pp=0; pp<IP; pp++){
        item1 = item_pairs(pp,0);
        item2 = item_pairs(pp,1);
        x = LAM(item1,_);
        y = LAM(item2,_);
        cov_item_pair[pp] = immer_bilinear_form( x, PHI, y) + PSI(item1,item2);
        cor_item_pair[pp] = cov_item_pair[pp] / sd_item[item1] / sd_item[item2];
        rho(item1, item2) = cor_item_pair[pp];
        rho(item2, item1) = cor_item_pair[pp];
    }

    //--- output
    return Rcpp::List::create(
            Rcpp::Named("var_item") = var_item,   // 1
            Rcpp::Named("sd_item") = sd_item,     // 2
            Rcpp::Named("cov_item_pair") = cov_item_pair,  // 3
            Rcpp::Named("cor_item_pair") = cor_item_pair,   // 4
            Rcpp::Named("rho") = rho,    // 5
            Rcpp::Named("tau") = tau     // 6
        );
}
///********************************************************************

///********************************************************************
///** immer_cmml_trafo_variances_covariances
// [[Rcpp::export]]
Rcpp::List immer_cmml_trafo_variances_covariances(Rcpp::NumericMatrix LAM,
        Rcpp::NumericMatrix PHI, Rcpp::NumericMatrix PSI, Rcpp::NumericMatrix GAM,
        Rcpp::IntegerMatrix item_pairs, Rcpp::IntegerMatrix LAM_index,
        int N_LAM, Rcpp::IntegerMatrix PHI_index, int N_PHI, Rcpp::IntegerMatrix PSI_index,
        int N_PSI, Rcpp::IntegerMatrix tau_index, int N_tau, Rcpp::IntegerMatrix GAM_index,
        int N_GAM, Rcpp::IntegerMatrix rho_index, int N_rho )
{
    int D = PHI.ncol();
    int I = LAM.nrow();
    int IP = item_pairs.nrow();
    Rcpp::NumericVector var_item(I);
    Rcpp::NumericVector sd_item(I);
    Rcpp::NumericVector cov_item_pair(IP);
    Rcpp::NumericVector cor_item_pair(IP);
    Rcpp::NumericVector x(D);
    Rcpp::NumericVector y(D);
    int item1, item2;

    // item-wise variances and standard deviations
    for (int ii=0; ii<I; ii++){
        x = LAM(ii,_);
        var_item[ii] = immer_bilinear_form( x, PHI, x) + 1;
        sd_item[ii] = std::sqrt( var_item[ii] );
    }

    // item-pair-wise covariances and correlations
    for (int pp=0; pp<IP; pp++){
        item1 = item_pairs(pp,0);
        item2 = item_pairs(pp,1);
        x = LAM(item1,_);
        y = LAM(item2,_);
        cov_item_pair[pp] = immer_bilinear_form( x, PHI, y) + PSI(item1,item2);
        cor_item_pair[pp] = cov_item_pair[pp] / sd_item[item1] / sd_item[item2];
    }

    //--- output
    return Rcpp::List::create(
            Rcpp::Named("var_item") = var_item,   // 1
            Rcpp::Named("sd_item") = sd_item,     // 2
            Rcpp::Named("cov_item_pair") = cov_item_pair,  // 3
            Rcpp::Named("cor_item_pair") = cor_item_pair   // 4
        );

}
///********************************************************************

///********************************************************************
///** immer_cmml_trafo_irt_parameters_derivative
// [[Rcpp::export]]
Rcpp::List immer_cmml_trafo_irt_parameters_derivative(Rcpp::NumericMatrix LAM,
        Rcpp::NumericMatrix PHI, Rcpp::NumericMatrix PSI, Rcpp::NumericMatrix GAM,
        Rcpp::IntegerMatrix item_pairs, Rcpp::IntegerMatrix LAM_index,
        int N_LAM, Rcpp::IntegerMatrix PHI_index, int N_PHI, Rcpp::IntegerMatrix PSI_index,
        int N_PSI, Rcpp::IntegerMatrix tau_index, int N_tau, Rcpp::IntegerMatrix GAM_index,
        int N_GAM, Rcpp::IntegerMatrix rho_index, int N_rho )
{
    int D = PHI.ncol();
    int I = LAM.nrow();
    int IP = item_pairs.nrow();
    int K = GAM.ncol();

    Rcpp::NumericVector var_item(I);
    Rcpp::NumericVector sd_item(I);
    Rcpp::NumericVector cov_item_pair(IP);
    Rcpp::NumericVector cor_item_pair(IP);
    Rcpp::NumericMatrix rho(I,I);
    Rcpp::NumericMatrix tau(I,K);
    Rcpp::NumericMatrix cov_item_pair_der_LAM(IP,N_LAM);
    Rcpp::NumericMatrix cov_item_pair_der_PHI(IP,N_PHI);
    Rcpp::NumericMatrix cov_item_pair_der_PSI(IP,N_PSI);
    Rcpp::NumericMatrix var_item_der_LAM(I,N_LAM);
    Rcpp::NumericMatrix var_item_der_PHI(I,N_PHI);
    Rcpp::NumericMatrix tau_der_GAM(N_tau,N_GAM);
    Rcpp::NumericMatrix tau_der_LAM(N_tau,N_LAM);
    Rcpp::NumericMatrix tau_der_PHI(N_tau,N_PHI);
    Rcpp::NumericMatrix rho_der_LAM(N_rho,N_LAM);
    Rcpp::NumericMatrix rho_der_PHI(N_rho,N_PHI);
    Rcpp::NumericMatrix rho_der_PSI(N_rho,N_PSI);
    Rcpp::NumericVector x(D);
    Rcpp::NumericVector y(D);
    int item1, item2;
    double val, val1;

    // item-wise variances and standard deviations
    for (int ii=0; ii<I; ii++){
        x = LAM(ii,_);
        var_item[ii] = immer_bilinear_form( x, PHI, x) + 1;
        sd_item[ii] = std::sqrt( var_item[ii] );
    }

    // item-wise thresholds
    for (int ii=0; ii<I; ii++){
        tau(ii,0) = -999;
        for (int hh=1; hh<K-1; hh++){
            tau(ii,hh) = GAM(ii,hh) / sd_item[ii];
        }
        tau(ii,K-1)=999;
    }

    // item-pair-wise covariances and correlations
    for (int pp=0; pp<IP; pp++){
        item1 = item_pairs(pp,0);
        item2 = item_pairs(pp,1);
        x = LAM(item1,_);
        y = LAM(item2,_);
        cov_item_pair[pp] = immer_bilinear_form( x, PHI, y) + PSI(item1,item2);
        cor_item_pair[pp] = cov_item_pair[pp] / sd_item[item1] / sd_item[item2];
        rho(item1, item2) = cor_item_pair[pp];
        rho(item2, item1) = cor_item_pair[pp];
    }

    //--- derivatives covariances LAM
    for (int pp=0; pp<IP; pp++){
        item1 = item_pairs(pp,0);
        item2 = item_pairs(pp,1);
        // lambda derivative first item
        for (int dd=0; dd<D; dd++){
            val=0;
            for (int ee=0;ee<D;ee++){
                val += LAM(item2, ee)*PHI(dd,ee);
            }
            cov_item_pair_der_LAM(pp, LAM_index(item1,dd) ) = val;
        }
        // lambda derivative second item
        for (int ee=0; ee<D; ee++){
            val=0;
            for (int dd=0;dd<D;dd++){
                val += LAM(item1, dd)*PHI(ee,dd);
            }
            cov_item_pair_der_LAM(pp, LAM_index(item2,ee) ) = val;
        }
    }

    //--- derivatives covariances PHI
    for (int pp=0; pp<IP; pp++){
        item1 = item_pairs(pp,0);
        item2 = item_pairs(pp,1);
        for (int dd=0; dd<D; dd++){
            for(int ee=dd; ee<D; ee++){
                val = LAM(item1,dd)*LAM(item2,ee);
                if ( dd!= ee){
                    val += LAM(item1,ee)*LAM(item2,dd);
                }
                cov_item_pair_der_PHI(pp, PHI_index(dd,ee) ) = val;
            }
        }
    }

    //--- derivatives covariances PSI
    for (int pp=0; pp<IP; pp++){
        item1 = item_pairs(pp,0);
        item2 = item_pairs(pp,1);
        cov_item_pair_der_PSI(pp, PSI_index(item1,item2) ) = 1;
    }

    //--- derivatives variances LAM
    for (int ii=0; ii<I; ii++){
        for (int dd=0; dd<D; dd++){
            val=0;
            for (int ee=0; ee<D; ee++){
                val += 2 * LAM(ii,ee) * PHI(dd,ee);
            }
            var_item_der_LAM(ii, LAM_index(ii, dd) ) = val;
        }
    }

    //--- derivatives variances PHI
    for (int ii=0; ii<I; ii++){
        for (int dd=0; dd<D; dd++){
            val=0;
            for ( int ee=0; ee<D; ee++){
                val = LAM(ii,dd) * LAM(ii,ee);
                if (dd != ee){
                    val = 2*val;
                }
                var_item_der_PHI(ii, PHI_index(dd, ee) ) = val;
            }
        }
    }

    //--- derivatives tau with respect to GAM
    for (int ii=0; ii<I; ii++){
        val = 1 / sd_item[ii];
        for (int kk=1; kk<K-1; kk++){
            tau_der_GAM( tau_index(ii,kk), GAM_index(ii,kk) ) = val;
        }
    }

    //--- derivatives tau with respect to LAM and PHI
    for (int ii=0; ii<I; ii++){
        for (int kk=1; kk<K-1;kk++){
            val1 = - 0.5 * GAM(ii,kk) * std::pow( var_item[ii], -1.5 );
            //* derivatives LAM
            for (int dd=0; dd<D; dd++){
                val = val1 * var_item_der_LAM(ii, LAM_index(ii,dd) );
                tau_der_LAM( tau_index(ii,kk), LAM_index(ii,dd) ) = val;
            }
            //* derivatives PHI
            for (int dd=0; dd<D; dd++){
                for (int ee=dd; ee<D; ee++){
                    val = val1 * var_item_der_PHI(ii, PHI_index(dd,ee) );
                    tau_der_PHI( tau_index(ii,kk), PHI_index(dd,ee) ) = val;
                }
            }
        }  // end kk
    }  // end ii

    //--- derivatives rho with respect to LAM
    int LAM_index_dd;
    for (int pp=0; pp<IP; pp++){
        item1 = item_pairs(pp,0);
        item2 = item_pairs(pp,1);
        for (int dd=0; dd<D; dd++){
            for (int hh=0; hh<2;hh++){
                LAM_index_dd = LAM_index(item_pairs(pp,hh), dd);
                // apply differentiation chain rule for cov_ij * var_i^-.5 * var_j^-.5
                val = immer_derivative_correlation( cov_item_pair[pp], var_item[ item1 ], var_item[ item2 ],
                        cov_item_pair_der_LAM(pp, LAM_index_dd), var_item_der_LAM(item1, LAM_index_dd),
                        var_item_der_LAM(item2, LAM_index_dd) );
                rho_der_LAM( rho_index(item1, item2), LAM_index_dd ) = val;
            }
        }
    }

    //--- derivatives rho with respect to PSI
    int PSI_index_dd;
    for (int pp=0; pp<IP; pp++){
        item1 = item_pairs(pp,0);
        item2 = item_pairs(pp,1);
        PSI_index_dd = PSI_index(item1, item2);
        val = cov_item_pair_der_PSI(pp, PSI_index_dd) / sd_item[ item1 ] / sd_item[ item2 ];
        rho_der_PSI( rho_index(item1, item2), PSI_index_dd ) = val;
    }

    //--- derivatives rho with respect to PHI
    int PHI_index_dd;
    for (int pp=0; pp<IP; pp++){
        item1 = item_pairs(pp,0);
        item2 = item_pairs(pp,1);
        for (int dd=0; dd<D; dd++){
            for (int ee=dd; ee<D; ee++){
                PHI_index_dd = PHI_index(dd, ee);
                val = immer_derivative_correlation( cov_item_pair[pp], var_item[ item1 ], var_item[ item2 ],
                        cov_item_pair_der_PHI(pp, PHI_index_dd), var_item_der_PHI(item1, PHI_index_dd),
                        var_item_der_PHI(item2, PHI_index_dd) );
                rho_der_PHI( rho_index(item1, item2), PHI_index_dd ) = val;
            }
        }
    }

    //--- output
    return Rcpp::List::create(
            Rcpp::Named("var_item") = var_item,   // 1
            Rcpp::Named("sd_item") = sd_item,     // 2
            Rcpp::Named("cov_item_pair") = cov_item_pair,  // 3
            Rcpp::Named("cor_item_pair") = cor_item_pair,   // 4
            Rcpp::Named("rho") = rho,    // 5
            Rcpp::Named("tau") = tau,    // 6
            Rcpp::Named("cov_item_pair_der_LAM") = cov_item_pair_der_LAM,  // 7
            Rcpp::Named("cov_item_pair_der_PHI") = cov_item_pair_der_PHI,  // 8
            Rcpp::Named("cov_item_pair_der_PSI") = cov_item_pair_der_PSI,  // 9
            Rcpp::Named("var_item_der_LAM") = var_item_der_LAM,    // 10
            Rcpp::Named("var_item_der_PHI") = var_item_der_PHI,    // 11
            Rcpp::Named("tau_der_GAM") = tau_der_GAM,    // 12
            Rcpp::Named("tau_der_LAM") = tau_der_LAM,    // 13
            Rcpp::Named("tau_der_PHI") = tau_der_PHI,    // 14
            Rcpp::Named("rho_der_LAM") = rho_der_LAM,    // 15
            Rcpp::Named("rho_der_PHI") = rho_der_PHI,    // 16
            Rcpp::Named("rho_der_PSI") = rho_der_PSI    // 17
        );

}
///********************************************************************

///********************************************************************
///** immer_sum_product
// computes sum_v ( x[v] * y[v])
double immer_sum_product(Rcpp::NumericVector x, Rcpp::NumericVector y)
{
    int n=x.size();
    double val=0;
    for (int nn=0; nn<n; nn++){
        if ( ( x[nn] !=0 ) & ( y[nn] != 0) ){
            val += x[nn]*y[nn];
        }
    }
    return val;
}
///********************************************************************

///********************************************************************
///** immer_sum_matrix_product
Rcpp::NumericMatrix immer_sum_matrix_product(Rcpp::NumericMatrix x1, Rcpp::NumericMatrix x2,
        Rcpp::NumericMatrix y1, Rcpp::NumericMatrix y2)
{
    int IP=x1.nrow();
    int N_LAM=y1.ncol();
    Rcpp::NumericMatrix z(IP, N_LAM);
    double val=0;
    for (int pp=0; pp<IP; pp++){
        for (int vv=0; vv<N_LAM; vv++){
            val=0;
            val += immer_sum_product( x1(pp,_), y1(_,vv) );
            val += immer_sum_product( x2(pp,_), y2(_,vv) );
            z(pp, vv) = val;
        }
    }
    return z;
}
///********************************************************************



///********************************************************************
///** immer_cmml_calc_prob_pars
// [[Rcpp::export]]
Rcpp::List immer_cmml_calc_prob_pars(Rcpp::NumericMatrix LAM,
        Rcpp::NumericMatrix PHI, Rcpp::NumericMatrix PSI, Rcpp::NumericMatrix GAM,
        Rcpp::IntegerMatrix item_pairs, Rcpp::IntegerMatrix LAM_index,
        int N_LAM, Rcpp::IntegerMatrix PHI_index, int N_PHI, Rcpp::IntegerMatrix PSI_index,
        int N_PSI, Rcpp::IntegerMatrix tau_index, int N_tau, Rcpp::IntegerMatrix GAM_index,
        int N_GAM, Rcpp::IntegerMatrix rho_index, int N_rho, Rcpp::IntegerMatrix item_table )
{

    // compute parameters rho and tau
    Rcpp::List res = immer_cmml_trafo_irt_parameters(LAM, PHI, PSI, GAM, item_pairs, LAM_index,
            N_LAM, PHI_index, N_PHI, PSI_index, N_PSI, tau_index, N_tau, GAM_index,
            N_GAM, rho_index, N_rho );
    Rcpp::NumericMatrix rho = res["rho"];
    Rcpp::NumericMatrix tau = res["tau"];

    // compute probabilities
    Rcpp::NumericVector probs = immer_cmml_calc_probs(tau, rho, item_table);

    // output
    return Rcpp::List::create(
            Rcpp::Named("rho") = rho,
            Rcpp::Named("tau") = tau,
            Rcpp::Named("probs") = probs
        );
}
///********************************************************************



// derivatives and chain rules

// dF_ijhk / d(gam) = dF_ijhk / d(tau) * d(tau) / d(gam)
// dF_ijhk / d(lam) = dF_ijhk / d(tau) * d(tau) / d(lam) + dF_ijhk / d(rho) * d(rho) / d(lam)
// dF_ijhk / d(phi) = dF_ijhk / d(tau) * d(tau) / d(phi) + dF_ijhk / d(rho) * d(rho) / d(phi)
// dF_ijhk / d(psi) = dF_ijhk / d(rho) * d(rho) / d(phi)


///********************************************************************
///** immer_cmml_derivative_fa_parameters
// [[Rcpp::export]]
Rcpp::List immer_cmml_derivative_fa_parameters(Rcpp::NumericMatrix LAM,
        Rcpp::NumericMatrix PHI, Rcpp::NumericMatrix PSI, Rcpp::NumericMatrix GAM,
        Rcpp::IntegerMatrix item_pairs, Rcpp::IntegerMatrix LAM_index,
        int N_LAM, Rcpp::IntegerMatrix PHI_index, int N_PHI, Rcpp::IntegerMatrix PSI_index,
        int N_PSI, Rcpp::IntegerMatrix tau_index, int N_tau, Rcpp::IntegerMatrix GAM_index,
        int N_GAM, Rcpp::IntegerMatrix rho_index, int N_rho, Rcpp::IntegerMatrix item_table )
{
    // compute derivatives of tau, rho with respect to parameters
    Rcpp::List res = immer_cmml_trafo_irt_parameters_derivative(LAM, PHI, PSI, GAM, item_pairs, LAM_index,
                        N_LAM, PHI_index, N_PHI, PSI_index, N_PSI, tau_index, N_tau,
                        GAM_index, N_GAM, rho_index, N_rho);
    Rcpp::NumericMatrix rho = res["rho"];
    Rcpp::NumericMatrix tau = res["tau"];
    Rcpp::NumericMatrix tau_der_GAM = res["tau_der_GAM"];
    Rcpp::NumericMatrix tau_der_LAM = res["tau_der_LAM"];
    Rcpp::NumericMatrix tau_der_PHI = res["tau_der_PHI"];
    Rcpp::NumericMatrix rho_der_LAM = res["rho_der_LAM"];
    Rcpp::NumericMatrix rho_der_PHI = res["rho_der_PHI"];
    Rcpp::NumericMatrix rho_der_PSI = res["rho_der_PSI"];

    // calculate probabilities
    Rcpp::NumericVector probs = immer_cmml_calc_probs(tau, rho, item_table);

    // derivative of probabilities
    Rcpp::List res1 = immer_cmml_probs_derivatives_tau_rho(tau, rho, item_table, rho_index, N_rho, tau_index, N_tau );
    Rcpp::NumericMatrix probs_der_rho = res1["probs_der_rho"];
    Rcpp::NumericMatrix probs_der_tau = res1["probs_der_tau"];

    //************************************
    //****** derivatives
    int IP = item_table.nrow();
    int K = tau_index.ncol();
    Rcpp::NumericMatrix probs_der_GAM(IP, N_GAM);
    Rcpp::NumericMatrix probs_der_PSI(IP, N_PSI);
    int item, item1, item2, tau_index_temp, GAM_index_temp, PSI_index_temp, rho_index_temp;
    double val;

    //--- derivative with respect to GAM
    // d(probs) / d(gam) = d(probs) / d(tau) * d(tau) / d(gam)
    for (int pp=0; pp<IP; pp++){
        for (int ii=0; ii<2; ii++){
            item = item_table(pp,ii);
            for (int kk=1; kk < K-1; kk++){
                tau_index_temp = tau_index( item, kk);
                GAM_index_temp = GAM_index( item, kk);
                val = probs_der_tau(pp, tau_index_temp) * tau_der_GAM( tau_index_temp, GAM_index_temp);
                probs_der_GAM(pp, GAM_index_temp) = val;
            }
        }
    }

    //--- derivative with respect to LAM
    // d(probs) / d(lam) = d(probs) / d(tau) * d(tau) / d(lam)  +
    //                        d(probs) / d(rho) * d(rho) / d(lam)
    Rcpp::NumericMatrix probs_der_LAM = immer_sum_matrix_product( probs_der_tau,
                probs_der_rho, tau_der_LAM, rho_der_LAM);

    //--- derivative with respect to PHI
    // d(probs) / d(phi) = d(probs) / d(tau) * d(tau) / d(phi)  +
    //                        d(probs) / d(rho) * d(rho) / d(phi)
    Rcpp::NumericMatrix probs_der_PHI = immer_sum_matrix_product( probs_der_tau,
            probs_der_rho, tau_der_PHI, rho_der_PHI);

    //--- derivative with respect to PSI
    // d(probs) / d(psi) =  d(probs) / d(rho) * d(rho) / d(psi)
    for (int pp=0; pp<IP; pp++){
        item1 = item_table(pp,0);
        item2 = item_table(pp,1);
        PSI_index_temp = PSI_index(item1,item2);
        rho_index_temp = rho_index(item1,item2);
        val = probs_der_rho(pp, rho_index_temp)*rho_der_PSI(rho_index_temp, PSI_index_temp);
        probs_der_PSI(pp, rho_index_temp) = val;
    }

    // output
    return Rcpp::List::create(
            Rcpp::Named("rho") = rho,
            Rcpp::Named("tau") = tau,
            Rcpp::Named("probs") = probs,
            Rcpp::Named("probs_der_rho") = probs_der_rho,
            Rcpp::Named("probs_der_tau") = probs_der_tau,
            Rcpp::Named("probs_der_GAM") = probs_der_GAM,
            Rcpp::Named("probs_der_LAM") = probs_der_LAM,
            Rcpp::Named("probs_der_PHI") = probs_der_PHI,
            Rcpp::Named("probs_der_PSI") = probs_der_PSI,
            Rcpp::Named("tau_der_GAM") = tau_der_GAM,
            Rcpp::Named("tau_der_LAM") = tau_der_LAM,
            Rcpp::Named("tau_der_PHI") = tau_der_PHI,
            Rcpp::Named("rho_der_LAM") = rho_der_LAM,
            Rcpp::Named("rho_der_PHI") = rho_der_PHI,
            Rcpp::Named("rho_der_PSI") = rho_der_PSI
        );

}
///********************************************************************

///********************************************************************
///** immer_sparse_matrix_create_index
Rcpp::NumericMatrix immer_sparse_matrix_create_index(Rcpp::NumericMatrix x)
{
    int N=x.nrow();
    int P=x.ncol();
    int V=N*P;
    int index=0;
    Rcpp::NumericMatrix x_sparse(V,3);
    for (int nn=0; nn<N; nn++){
        for (int cc=0; cc<P; cc++){
            if ( x(nn,cc) != 0 ){
                x_sparse(index,0) = nn;
                x_sparse(index,1) = cc;
                x_sparse(index,2) = x(nn,cc);
                index ++;
            }
        }
    }
    x_sparse = x_sparse( Rcpp::Range(0, index), Rcpp::Range(0, 2) );
    return x_sparse;
}
///********************************************************************

///********************************************************************
///** immer_sparse_matrix_mat_mult_vec
Rcpp::NumericVector immer_sparse_matrix_mat_mult_vec(Rcpp::NumericMatrix x_sparse,
        Rcpp::NumericVector vec, int NX)
{
    int NS=x_sparse.nrow();
    Rcpp::NumericVector y(NX);
    y.fill(0);
    for (int ii=0; ii<NS-1; ii++){
        y[ x_sparse(ii,0) ] += x_sparse(ii,2) * vec[ x_sparse(ii,1) ];
    }
    // include argument for exponentiation
    // do y = exp(y);
    return y;
}
///********************************************************************


///********************************************************************
///** immer_matrix_mat_to_vec
Rcpp::NumericVector immer_matrix_mat_to_vec(Rcpp::NumericMatrix x, Rcpp::IntegerMatrix x_index )
{
    int NR = x.nrow();
    int NC = x.ncol();
    int Ny = NR*NC;
    Rcpp::NumericVector y(Ny);
    int maxval = 0;
    int x_temp;
    for (int rr=0; rr<NR; rr++){
        for (int cc=0; cc<NC; cc++){
            x_temp = x_index(rr,cc);
            if (x_temp >=0){
                y[ x_temp ] = x(rr,cc);
                if (x_temp > maxval){
                    maxval = x_temp;
                }
            }
        }
    }
    y = y[ Rcpp::Range(0, maxval) ];
    return y;
}
///********************************************************************


///********************************************************************
///** immer_matrix_vec_to_mat
Rcpp::NumericMatrix immer_matrix_vec_to_mat(Rcpp::NumericVector par,
        Rcpp::IntegerMatrix x_index, Rcpp::NumericMatrix x_init )
{
    int NR = x_index.nrow();
    int NC = x_index.ncol();
    int x_temp;
    Rcpp::NumericMatrix x(NR,NC);
    for (int rr=0; rr<NR; rr++){
        for (int cc=0; cc<NC; cc++){
            x_temp = x_index(rr,cc);
            x(rr,cc) = x_init(rr,cc);
            if (x_temp >=0){
                x(rr,cc) = par[ x_temp ] + x(rr,cc);
            }
        }
    }
    return x;
}
///********************************************************************


//  LAM = \sum_v w_v l_v
//  dF / d(l_v) = dF / d(lam) * d(lam) / d(l_v) = dF / d(lam) * w_v

//  LAM = exp( \sum_v w_v l_v )
//  dF / d(l_v) = dF / d(lam) * d(lam) / d(l_v) = dF / d(lam) * lam * w_v



///********************************************************************
///** immer_cmml_derivative_designmatrix
// [[Rcpp::export]]
Rcpp::NumericVector immer_cmml_derivative_designmatrix(
        Rcpp::NumericMatrix probs_der_par, Rcpp::NumericMatrix W_par,
        Rcpp::NumericVector par, Rcpp::NumericVector do_log )
{
    int IP = probs_der_par.nrow();
    int NP = W_par.ncol();
    int NW = W_par.nrow();
    Rcpp::NumericMatrix probs_der_basispar(IP,NP);

    Rcpp::NumericVector do_log1(NW);
    int NL = do_log.size();
    for (int ww=0; ww<NW; ww++){
        if (NL < NW ){
            do_log1[ww] = do_log[0];
        } else {
            do_log1[ww] = do_log[ww];
        }
    }

    probs_der_basispar.fill(0);
    double val=0;
    for (int pp=0; pp<IP; pp++){
        for (int hh=0; hh<NP; hh++){
            for (int vv=0; vv<NW; vv++){
                if (W_par(vv,hh) !=0 ){
                    if ( ! do_log1[vv] ){
                        val = probs_der_par(pp,vv) * W_par(vv,hh);
                    } else {
                        val = probs_der_par(pp,vv) * W_par(vv,hh) * par[ vv ];
                    }
                    probs_der_basispar(pp,hh) += val;
                }
            }
        }
    }
    return probs_der_basispar;
}
///********************************************************************



///********************************************************************
///** immer_cmml_basispar_types_to_full_pars
Rcpp::NumericMatrix immer_cmml_basispar_types_to_full_pars( Rcpp::NumericMatrix W_par,
        Rcpp::NumericVector par_basispar, Rcpp::NumericMatrix par_init, Rcpp::IntegerMatrix par_index,
        Rcpp::NumericVector do_log_par )
{
    int NR = W_par.nrow();
    // compute vector
    Rcpp::NumericVector par_fullvec(NR);
    for (int ww=0; ww<NR; ww++){
        par_fullvec[ww] = immer_sum_product( W_par(ww,_), par_basispar );
        if ( do_log_par[0] == 1){
            par_fullvec[ww] = std::exp( par_fullvec[ww] );
        }
    }
    // compute matrix
    Rcpp::NumericMatrix par_full = immer_matrix_vec_to_mat( par_fullvec, par_index, par_init );
    //-- output
    return par_full;

}
///********************************************************************


///********************************************************************
///** immer_cmml_basispar_to_pars_types_helper
Rcpp::NumericVector immer_cmml_basispar_to_pars_types_helper( Rcpp::NumericVector basispar,
    Rcpp::IntegerVector design_temp, int Nvv )
{
    if (Nvv==0){ Nvv=1;}
    Rcpp::NumericVector par_basispar(Nvv);
    int ND=design_temp.size();
    int gg=0;
    for (int dd=0; dd<ND; dd++){
        if (design_temp[dd] >= 0 ){
            par_basispar[gg] = basispar[ design_temp[dd] ];
            gg ++;
        }
    }
    return par_basispar;
}
///********************************************************************



///********************************************************************
///** immer_cmml_basispar_to_pars_types
Rcpp::List immer_cmml_basispar_to_pars_types( Rcpp::NumericVector basispar,
    Rcpp::IntegerMatrix basispar_design, Rcpp::IntegerVector basispar_length )
{
    // # LAM_basispar
    // # GAM_basispar
    // # PHI_basispar
    // # PSI_basispar
    int vv=0;
    // LAM
    vv = 0;
    Rcpp::NumericVector LAM_basispar = immer_cmml_basispar_to_pars_types_helper( basispar,
                        basispar_design(_,vv), basispar_length[vv] );
    // GAM
    vv = 1;
    Rcpp::NumericVector GAM_basispar = immer_cmml_basispar_to_pars_types_helper( basispar,
                        basispar_design(_,vv), basispar_length[vv] );
    // PHI
    vv = 2;
    Rcpp::NumericVector PHI_basispar = immer_cmml_basispar_to_pars_types_helper( basispar,
                        basispar_design(_,vv), basispar_length[vv] );
    // PSI
    vv = 3;
    Rcpp::NumericVector PSI_basispar = immer_cmml_basispar_to_pars_types_helper( basispar,
                        basispar_design(_,vv), basispar_length[vv] );

    //--- output
    return Rcpp::List::create(
            Rcpp::Named("basispar") = basispar,
            Rcpp::Named("LAM_basispar") = LAM_basispar,
            Rcpp::Named("GAM_basispar") = GAM_basispar,
            Rcpp::Named("PHI_basispar") = PHI_basispar,
            Rcpp::Named("PSI_basispar") = PSI_basispar
        );
}
///********************************************************************


///********************************************************************
///** immer_cmml_basispar_to_parameters
Rcpp::List immer_cmml_basispar_to_parameters( Rcpp::NumericVector basispar,
    Rcpp::IntegerMatrix basispar_design, Rcpp::IntegerVector basispar_length,
    Rcpp::NumericMatrix W_LAM, Rcpp::NumericMatrix LAM_init, Rcpp::IntegerMatrix LAM_index,
    Rcpp::NumericMatrix W_GAM, Rcpp::NumericMatrix GAM_init, Rcpp::IntegerMatrix GAM_index,
    Rcpp::NumericMatrix W_PHI, Rcpp::NumericMatrix PHI_init, Rcpp::IntegerMatrix PHI_index,
    Rcpp::NumericMatrix W_PSI, Rcpp::NumericMatrix PSI_init, Rcpp::IntegerMatrix PSI_index,
    Rcpp::NumericVector do_log_LAM, Rcpp::NumericVector do_log_GAM, Rcpp::NumericVector do_log_PHI,
    Rcpp::NumericVector do_log_PSI
    )
{
    // select basis parameters
    Rcpp::List res = immer_cmml_basispar_to_pars_types( basispar, basispar_design, basispar_length );
    Rcpp::NumericVector LAM_basispar = res["LAM_basispar"];
    Rcpp::NumericVector GAM_basispar = res["GAM_basispar"];
    Rcpp::NumericVector PHI_basispar = res["PHI_basispar"];
    Rcpp::NumericVector PSI_basispar = res["PSI_basispar"];

    // LAM
    Rcpp::NumericMatrix LAM = immer_cmml_basispar_types_to_full_pars( W_LAM, LAM_basispar, LAM_init,
                    LAM_index, do_log_LAM );
    // GAM
    Rcpp::NumericMatrix GAM = immer_cmml_basispar_types_to_full_pars( W_GAM, GAM_basispar, GAM_init,
                    GAM_index, do_log_GAM  );
    // PHI
    Rcpp::NumericMatrix PHI = immer_cmml_basispar_types_to_full_pars( W_PHI, PHI_basispar, PHI_init,
                    PHI_index, do_log_PHI );
    // PSI
    Rcpp::NumericMatrix PSI = immer_cmml_basispar_types_to_full_pars( W_PSI, PSI_basispar, PSI_init,
                    PSI_index, do_log_PSI );

    //--- output
    return Rcpp::List::create(
            Rcpp::Named("basispar") = basispar,
            Rcpp::Named("LAM_basispar") = LAM_basispar,
            Rcpp::Named("LAM") = LAM,
            Rcpp::Named("GAM_basispar") = GAM_basispar,
            Rcpp::Named("GAM") = GAM,
            Rcpp::Named("PHI_basispar") = PHI_basispar,
            Rcpp::Named("PHI") = PHI,
            Rcpp::Named("PSI_basispar") = PSI_basispar,
            Rcpp::Named("PSI") = PSI
        );
}
///********************************************************************

///********************************************************************
///** immer_cmml_basispar_to_probs
// [[Rcpp::export]]
Rcpp::List immer_cmml_basispar_to_probs( Rcpp::NumericVector basispar,
    Rcpp::IntegerMatrix basispar_design, Rcpp::IntegerVector basispar_length,
    Rcpp::NumericMatrix W_LAM, Rcpp::NumericMatrix LAM_init, Rcpp::IntegerMatrix LAM_index,
    Rcpp::NumericMatrix W_GAM, Rcpp::NumericMatrix GAM_init, Rcpp::IntegerMatrix GAM_index,
    Rcpp::NumericMatrix W_PHI, Rcpp::NumericMatrix PHI_init, Rcpp::IntegerMatrix PHI_index,
    Rcpp::NumericMatrix W_PSI, Rcpp::NumericMatrix PSI_init, Rcpp::IntegerMatrix PSI_index,
    int N_LAM, int N_GAM, int N_PHI, int N_PSI, Rcpp::IntegerMatrix item_pairs,
    Rcpp::IntegerMatrix item_table, Rcpp::IntegerMatrix tau_index, Rcpp::IntegerMatrix rho_index,
    int N_rho, int N_tau, Rcpp::NumericVector do_log_LAM, Rcpp::NumericVector do_log_GAM,
    Rcpp::NumericVector do_log_PHI,    Rcpp::NumericVector do_log_PSI )
{
    // convert basis parameters to factor-analytic model parameters
    Rcpp::List res = immer_cmml_basispar_to_parameters( basispar, basispar_design, basispar_length, W_LAM,
                        LAM_init, LAM_index, W_GAM, GAM_init, GAM_index, W_PHI, PHI_init, PHI_index,
                        W_PSI, PSI_init, PSI_index, do_log_LAM, do_log_GAM, do_log_PHI, do_log_PSI );
    Rcpp::NumericMatrix LAM = res["LAM"];
    Rcpp::NumericMatrix GAM = res["GAM"];
    Rcpp::NumericMatrix PHI = res["PHI"];
    Rcpp::NumericMatrix PSI = res["PSI"];

    // compute probabilities
    Rcpp::List res11 = immer_cmml_calc_prob_pars(LAM, PHI, PSI, GAM, item_pairs, LAM_index, N_LAM, PHI_index,
                    N_PHI, PSI_index, N_PSI, tau_index, N_tau, GAM_index, N_GAM, rho_index, N_rho,
                    item_table);
    Rcpp::NumericVector probs = res11["probs"];
    Rcpp::NumericMatrix tau = res11["tau"];
    Rcpp::NumericMatrix rho = res11["rho"];

    //--- output
    return Rcpp::List::create(
            Rcpp::Named("basispar") = basispar,
            Rcpp::Named("LAM") = LAM,
            Rcpp::Named("GAM") = GAM,
            Rcpp::Named("PHI") = PHI,
            Rcpp::Named("PSI") = PSI,
            Rcpp::Named("probs") = probs,
            Rcpp::Named("tau") = tau,
            Rcpp::Named("rho") = rho
        );
}
///********************************************************************

///********************************************************************
///** immer_cmml_basispar_to_derivatives
// [[Rcpp::export]]
Rcpp::List immer_cmml_basispar_to_derivatives( Rcpp::NumericVector basispar,
    Rcpp::IntegerMatrix basispar_design, Rcpp::IntegerVector basispar_length,
    Rcpp::NumericMatrix W_LAM, Rcpp::NumericMatrix LAM_init, Rcpp::IntegerMatrix LAM_index,
    Rcpp::NumericMatrix W_GAM, Rcpp::NumericMatrix GAM_init, Rcpp::IntegerMatrix GAM_index,
    Rcpp::NumericMatrix W_PHI, Rcpp::NumericMatrix PHI_init, Rcpp::IntegerMatrix PHI_index,
    Rcpp::NumericMatrix W_PSI, Rcpp::NumericMatrix PSI_init, Rcpp::IntegerMatrix PSI_index,
    int N_LAM, int N_GAM, int N_PHI, int N_PSI, Rcpp::IntegerMatrix item_pairs,
    Rcpp::IntegerMatrix item_table, Rcpp::IntegerMatrix tau_index, Rcpp::IntegerMatrix rho_index,
    int N_rho, int N_tau, Rcpp::NumericVector do_log_LAM, Rcpp::NumericVector do_log_GAM,
    Rcpp::NumericVector do_log_PHI,    Rcpp::NumericVector do_log_PSI )
{
    // convert basis parameters to factor-analytic model parameters
    Rcpp::List res = immer_cmml_basispar_to_parameters( basispar, basispar_design, basispar_length, W_LAM,
                        LAM_init, LAM_index, W_GAM, GAM_init, GAM_index, W_PHI, PHI_init, PHI_index,
                        W_PSI, PSI_init, PSI_index, do_log_LAM, do_log_GAM, do_log_PHI, do_log_PSI);
    Rcpp::NumericMatrix LAM = res["LAM"];
    Rcpp::NumericMatrix GAM = res["GAM"];
    Rcpp::NumericMatrix PHI = res["PHI"];
    Rcpp::NumericMatrix PSI = res["PSI"];

    // compute derivatives with respect to FA parameters
    Rcpp::List res11 = immer_cmml_derivative_fa_parameters(LAM, PHI, PSI, GAM, item_pairs, LAM_index, N_LAM, PHI_index,
                N_PHI, PSI_index, N_PSI, tau_index, N_tau, GAM_index, N_GAM, rho_index, N_rho,
                item_table);
    Rcpp::NumericMatrix probs_der_LAM = res11["probs_der_LAM"];
    Rcpp::NumericMatrix probs_der_GAM = res11["probs_der_GAM"];
    Rcpp::NumericMatrix probs_der_PHI = res11["probs_der_PHI"];
    Rcpp::NumericMatrix probs_der_PSI = res11["probs_der_PSI"];
    Rcpp::NumericVector probs = res11["probs"];

    // compute derivatives with respect to basis parameters
    // 1 # LAM_basispar
    // 2 # GAM_basispar
    // 3 # PHI_basispar
    // 4 # PSI_basispar
    int vv=0;
    int NP = basispar.size();
    int IP = item_table.nrow();
    Rcpp::NumericMatrix probs_der_basispar(IP, NP);
    int col_temp = 0;
    int NBL=0;

    //** LAM
    vv=0;
    Rcpp::NumericVector LAM_parvec = immer_matrix_mat_to_vec( LAM, LAM_index);
    Rcpp::NumericVector LAM_probs_der_basispar = immer_cmml_derivative_designmatrix( probs_der_LAM,
                                W_LAM, LAM_parvec, do_log_LAM );
    NBL = basispar_length[vv];
    for (int hh=0; hh<NBL; hh++){
        for (int ii=0; ii<IP; ii++){
            probs_der_basispar(ii, col_temp) = LAM_probs_der_basispar(ii, hh);
        }
        col_temp ++;
    }
    //** GAM
    vv=1;
    Rcpp::NumericVector GAM_parvec = immer_matrix_mat_to_vec( GAM, GAM_index);
    Rcpp::NumericVector GAM_probs_der_basispar = immer_cmml_derivative_designmatrix( probs_der_GAM,
                                W_GAM, GAM_parvec, do_log_GAM );
    NBL = basispar_length[vv];
    for (int hh=0; hh<NBL; hh++){
        for (int ii=0; ii<IP; ii++){
            probs_der_basispar(ii, col_temp) = GAM_probs_der_basispar(ii, hh);
        }
        col_temp ++;
    }
    //** PHI
    vv=2;
    Rcpp::NumericVector PHI_parvec = immer_matrix_mat_to_vec( PHI, PHI_index);
    Rcpp::NumericVector PHI_probs_der_basispar = immer_cmml_derivative_designmatrix( probs_der_PHI,
                                W_PHI, PHI_parvec, do_log_PHI );
    NBL = basispar_length[vv];
    for (int hh=0; hh<NBL; hh++){
        for (int ii=0; ii<IP; ii++){
            probs_der_basispar(ii, col_temp) = PHI_probs_der_basispar(ii, hh);
        }
        col_temp ++;
    }
    //** PSI
    vv=3;
    Rcpp::NumericVector PSI_parvec = immer_matrix_mat_to_vec( PSI, PSI_index);
    Rcpp::NumericVector PSI_probs_der_basispar = immer_cmml_derivative_designmatrix( probs_der_PSI,
                                W_PSI, PSI_parvec, do_log_PSI );
    NBL = basispar_length[vv];
    for (int hh=0; hh<NBL; hh++){
        for (int ii=0; ii<IP; ii++){
            probs_der_basispar(ii, col_temp) = PSI_probs_der_basispar(ii, hh);
        }
        col_temp ++;
    }

    //--- output
    return Rcpp::List::create(
            Rcpp::Named("probs") = probs,
            Rcpp::Named("LAM_probs_der_basispar") = LAM_probs_der_basispar,
            Rcpp::Named("GAM_probs_der_basispar") = GAM_probs_der_basispar,
            Rcpp::Named("PHI_probs_der_basispar") = PHI_probs_der_basispar,
            Rcpp::Named("PSI_probs_der_basispar") = PSI_probs_der_basispar,
            Rcpp::Named("probs_der_basispar") = probs_der_basispar
        );
}
///********************************************************************


// setup PML function


// GAM_par, W_GAM, c_GAM
// LAM_par, W_LAM, c_LAM
// PHI_par, W_PHI, c_PHI
// PSI_par, W_PSI, c_PSI

// par ... parameter vector with all free elements
// GAM_index to par_index
// LAM_index to par_index
// PHI_index to par_index
// PSI_index to par_index

// use corresponding design matrices ...



///********************************************************************
///** immer_cmml_proc_freq
// [[Rcpp::export]]
Rcpp::NumericMatrix immer_cmml_proc_freq( Rcpp::IntegerMatrix dat,
    Rcpp::IntegerMatrix dat_resp, int K, Rcpp::NumericVector weights)
{

    // item1: 0
    // item2: 1
    // cat1: 2
    // cat2: 3
    // n: 4
    // ll_index: 5

    int I = dat.ncol();
    int K1 = K+1;
    int N_pair = I * (I-1) / 2;
    int NR = N_pair * K1 * K1;
    int N = dat.nrow();
    Rcpp::NumericMatrix dfr(NR,6);
    Rcpp::NumericMatrix freq(K+1,K+1);
    int zz=0;
    int ll_temp=1;

    for (int ii1=0; ii1< (I-1); ii1++){
        for( int jj1=ii1+1; jj1 < I; jj1++){
            for (int nn=0; nn<N; nn++){
                if ( ( dat_resp(nn,ii1) == 1 ) & ( dat_resp(nn,jj1) == 1 ) ){
                    freq( dat(nn,ii1), dat(nn,jj1) ) += weights[nn];
                }
            }
            for (int rr=0; rr<K1; rr++){
                for (int cc=0; cc<K1; cc++){
                    dfr(zz,0) = ii1;
                    dfr(zz,1) = jj1;
                    dfr(zz,2) = rr;
                    dfr(zz,3) = cc;
                    dfr(zz,4) = freq(rr,cc);
                    freq(rr,cc) = 0;
                    dfr(zz,5) = ll_temp;
                    zz ++;
                }
            }
            ll_temp ++;
        }  // end jj1
    }    // end ii1

    //-- output
    return dfr;
}
///********************************************************************


///********************************************************************
///** immer_cmml_proc_generate_rho_index
// [[Rcpp::export]]
Rcpp::IntegerMatrix immer_cmml_proc_generate_rho_index( int I)
{
    Rcpp::IntegerMatrix rho_index(I,I);
    int ind = 1;
    for (int ii=0; ii<I-1; ii++){
        for (int jj=ii+1; jj<I; jj++){
            rho_index(ii,jj) = ind;
            ind ++;
        }
    }
    //-- output
    return rho_index;
}
///********************************************************************


///********************************************************************
///** immer_cmml_proc_generate_tau
// [[Rcpp::export]]
Rcpp::List immer_cmml_proc_generate_tau(Rcpp::IntegerVector maxK, int K)
{
    int I = maxK.size();
    Rcpp::IntegerMatrix tau(I,K+2);
    Rcpp::IntegerMatrix tau_index(I,K+2);
    tau.fill(999);
    tau_index.fill(1);
    int K_ii=0;
    double K2=0;
    int temp_count=2;

    for (int ii=0; ii<I; ii++){
        tau(ii,0) = -999;
        K_ii = maxK[ii];
        K2 = ( 1 + K) / 2;
        for (int kk=1; kk < ( K_ii+1); kk++){
            tau(ii,kk) = kk - K2;
            tau_index(ii,kk) = temp_count;
            temp_count ++;
        }
    }

    //--- output
    return Rcpp::List::create(
            Rcpp::Named("tau") = tau,
            Rcpp::Named("tau_index") = tau_index
        );
}
///********************************************************************

///********************************************************************
///** immer_cmml_proc_generate_LAM
// [[Rcpp::export]]
Rcpp::List immer_cmml_proc_generate_LAM(Rcpp::IntegerMatrix Q)
{
    int I = Q.nrow();
    int D = Q.ncol();
    Rcpp::IntegerMatrix LAM(I,D);
    Rcpp::IntegerMatrix LAM_index(I,D);
    int temp=1;
    int NW=0;

    for (int ii=0; ii<I; ii++){
        for (int dd=0; dd<D; dd++){
            LAM(ii,dd) = Q(ii,dd);
            if (Q(ii,dd) > 0){
                NW++;
            }
            LAM_index(ii,dd) = temp;
            temp ++;
        }
    }
    int ND1 = I*D;
    Rcpp::NumericMatrix W_LAM(ND1,NW);
    Rcpp::NumericVector LAM_basispar(NW);
    int ww=0;

    for (int ii=0; ii<I; ii++){
        for (int dd=0; dd<D; dd++){
            if (Q(ii,dd) > 0){
                W_LAM( LAM_index(ii,dd)-1, ww ) = 1;
                ww ++;
            }
        }
    }
    LAM_basispar.fill(.75);

    //--- output
    return Rcpp::List::create(
            Rcpp::Named("LAM") = LAM,
            Rcpp::Named("LAM_index") = LAM_index,
            Rcpp::Named("W_LAM") = W_LAM,
            Rcpp::Named("LAM_basispar") = LAM_basispar
        );
}
///********************************************************************

///********************************************************************
///** immer_cmml_proc_generate_PHI
// [[Rcpp::export]]
Rcpp::List immer_cmml_proc_generate_PHI(int D, bool use_diag)
{
    Rcpp::NumericMatrix PHI(D,D);
    Rcpp::IntegerMatrix PHI_index(D,D);
    Rcpp::NumericMatrix PHI_init(D,D);
    int H=D*(D+1)/2;
    if (! use_diag){ H = H - D; }
    Rcpp::IntegerMatrix W_PHI(H,H);
    Rcpp::NumericVector PHI_basispar(H);
    int hh=0;
    double val=1;
    int nd=0;
    if ( ! use_diag){ nd = 1; }

    int temp=1;
    for (int dd=0; dd<D; dd++){
        for (int ee=0; ee<=dd-nd; ee++){
            PHI_index(dd,ee) = temp;
            PHI_index(ee,dd) = temp;
            if (dd == ee){
                val=1;
            } else {
                val=0.05;
            }
            PHI(dd,ee)=val;
            PHI(ee,dd)=val;
            PHI_basispar[hh] = val;
            temp++;
            hh++;
        }
    }
    for (int hh=0; hh<H; hh++){
        W_PHI(hh,hh)=1;
    }

    //--- output
    return Rcpp::List::create(
            Rcpp::Named("PHI") = PHI,
            Rcpp::Named("PHI_index") = PHI_index,
            Rcpp::Named("PHI_init") = PHI_init,
            Rcpp::Named("W_PHI") = W_PHI,
            Rcpp::Named("PHI_basispar") = PHI_basispar
        );
}
///********************************************************************


//     return Rcpp::List::create(
//         Rcpp::Named("rho") = rho,
//         Rcpp::Named("thresh")=thresh
//         );

// Rcout << "rho1 " << rho1 << std::endl;
// Rcout << "rho2 " << rho2 << std::endl;

