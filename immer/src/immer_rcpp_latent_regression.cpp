//// File Name: immer_rcpp_latent_regression.cpp
//// File Version: 1.02


// [[Rcpp::depends(RcppArmadillo)]]

// #include <RcppArmadillo.h>
#include <Rcpp.h>
#include <Rmath.h>

using namespace Rcpp;
// using namespace arma;



///********************************************************************
///** immer_latent_regression_prior_normal
// [[Rcpp::export]]
Rcpp::NumericMatrix immer_latent_regression_prior_normal(Rcpp::NumericVector mu,
    Rcpp::NumericVector sigma, Rcpp::NumericVector theta)
{
    int TP = theta.size();
    int N = mu.size();
    Rcpp::NumericMatrix prior(N,TP);
    for (int nn=0; nn<N; nn++){
        for (int tt=0; tt<TP; tt++){
            prior(nn,tt) = R::dnorm( theta[tt], mu[nn], sigma[nn], FALSE);
        }
    }
    //-- output
    return prior;
}
///********************************************************************

///********************************************************************
///** immer_latent_regression_posterior
// [[Rcpp::export]]
Rcpp::List immer_latent_regression_posterior(Rcpp::NumericMatrix like,
    Rcpp::NumericMatrix prior, Rcpp::NumericVector weights )
{
    int TP = like.ncol();
    int N = like.nrow();
    Rcpp::NumericMatrix post(N,TP);
    Rcpp::NumericMatrix post_unnorm(N,TP);
    Rcpp::NumericVector indloglike(N);
    double temp=0;
    double loglike=0;
    double eps=1e-300;
    for (int nn=0; nn<N; nn++){
        temp=0;
        for (int tt=0; tt<TP; tt++){
            post(nn,tt) = like(nn,tt) * prior(nn,tt);
            post_unnorm(nn,tt) = post(nn,tt);
            temp += post(nn,tt);
        }
        indloglike[nn] = std::log( temp + eps);
        loglike += weights[nn] * indloglike[nn];
        for (int tt=0; tt<TP; tt++){
            post(nn,tt) = post(nn,tt) / temp;
        }
    }
    //-- output
    return Rcpp::List::create(
                Rcpp::Named("post") = post,
                Rcpp::Named("post_unnorm") = post_unnorm,
                Rcpp::Named("loglike") = loglike,
                Rcpp::Named("indloglike") = indloglike
            );
}
///********************************************************************

///********************************************************************
///** immer_latent_regression_calc_mu_sigma
// [[Rcpp::export]]
Rcpp::List immer_latent_regression_calc_mu_sigma(Rcpp::NumericMatrix X,
    Rcpp::IntegerVector group, int G, Rcpp::NumericVector beta,
    Rcpp::NumericVector gamma )
{
    int N = X.nrow();
    int NB = X.ncol();
    Rcpp::NumericVector mu(N);
    Rcpp::NumericVector sigma(N);
    double temp=0;
    for (int nn=0; nn<N; nn++){
        sigma[nn] = gamma[ group[nn] -1 ];
        temp=0;
        for (int bb=0;bb<NB; bb++){
            temp += X(nn,bb) * beta[bb];
        }
        mu[nn] = temp;
    }
    //-- output
    return Rcpp::List::create(
                Rcpp::Named("mu") = mu,
                Rcpp::Named("sigma") = sigma
            );
}
///********************************************************************

///********************************************************************
///** immer_add_elements
// [[Rcpp::export]]
Rcpp::NumericVector immer_add_elements( Rcpp::NumericVector x, int pos, double h)
{
    int V = x.size();
    Rcpp::NumericVector y(V);
    for (int vv=0; vv<V; vv++){
        y[vv] = x[vv];
    }
    if (pos >=0){
        y[pos] = x[pos] + h;
    }
    //-- output
    return y;
}
///********************************************************************

///********************************************************************
///** immer_latent_regression_calc_individual_likelihood
// [[Rcpp::export]]
Rcpp::NumericVector immer_latent_regression_calc_individual_likelihood(Rcpp::NumericMatrix X,
    Rcpp::IntegerVector group, int G, Rcpp::NumericVector pars,
    Rcpp::NumericVector theta, Rcpp::NumericVector weights,    Rcpp::NumericMatrix like)
{
    int NB = X.ncol();
    Rcpp::NumericVector beta(NB);
    Rcpp::NumericVector gamma(G);
    for (int bb=0; bb<NB; bb++){
        beta[bb] = pars[bb];
    }
    for (int gg=0; gg<G; gg++){
        gamma[gg] = pars[gg+NB];
    }

    // calculate predicted means and standard deviations
    Rcpp::List res1 = immer_latent_regression_calc_mu_sigma(X, group, G, beta, gamma);
    Rcpp::NumericVector mu = res1["mu"];
    Rcpp::NumericVector sigma = res1["sigma"];

    // compute prior
    Rcpp::NumericMatrix prior = immer_latent_regression_prior_normal( mu, sigma, theta);

    // calculate individual log-likelihood
    Rcpp::List res2 = immer_latent_regression_posterior(like, prior, weights);
    Rcpp::NumericVector indloglike = res2["indloglike"];

    //-- output
    return indloglike;
}
///********************************************************************

///********************************************************************
///** immer_latent_regression_calc_individual_likelihood_increment
// [[Rcpp::export]]
Rcpp::NumericVector immer_latent_regression_calc_individual_likelihood_increment(
    Rcpp::NumericMatrix X, Rcpp::IntegerVector group, int G, Rcpp::NumericVector pars,
    Rcpp::NumericVector theta, Rcpp::NumericVector weights,
    Rcpp::NumericMatrix like, int pos1, double h1, int pos2, double h2)
{
    Rcpp::NumericVector pars1 = immer_add_elements(pars, pos1, h1);
    pars1 = immer_add_elements(pars1, pos2, h2);
    Rcpp::NumericVector indloglike = immer_latent_regression_calc_individual_likelihood(
                X, group, G, pars1, theta, weights, like );
    //-- output
    return indloglike;
}
///********************************************************************


///********************************************************************
///** immer_latent_regression_vcov_xpd
// [[Rcpp::export]]
Rcpp::NumericMatrix immer_latent_regression_vcov_xpd(
    Rcpp::NumericMatrix X, Rcpp::IntegerVector group, int G,
    Rcpp::NumericVector pars, Rcpp::NumericVector theta, Rcpp::NumericVector weights,
    Rcpp::NumericMatrix like, double h)
{
    int NP = pars.size();
    int N = X.nrow();
    Rcpp::NumericMatrix indgrad(N,NP);
    Rcpp::NumericMatrix infomat(NP,NP);

    //-- compute individual gradients
    Rcpp::NumericVector indloglike1a(N);
    Rcpp::NumericVector indloglike1b(N);
    for (int pp=0; pp<NP; pp++){
        indloglike1a = immer_latent_regression_calc_individual_likelihood_increment(
                    X, group, G, pars, theta, weights, like, pp, h, -1,0 );
        indloglike1b = immer_latent_regression_calc_individual_likelihood_increment(
                    X, group, G, pars, theta, weights, like, pp, -h, -1,0 );
        for (int nn=0; nn<N; nn++){
            indgrad(nn,pp) = ( indloglike1a[nn] - indloglike1b[nn] ) / (2*h);
        }
    }
    //-- compute information matrix
    for (int pp1=0; pp1<NP; pp1++){
        for (int pp2=pp1; pp2<NP; pp2++){
            for (int nn=0; nn<N; nn++){
                infomat(pp1,pp2) += weights[nn] * indgrad(nn,pp1) * indgrad(nn,pp2);
            }
            infomat(pp2,pp1) = infomat(pp1,pp2);
        }
    }
    //-- output
    return infomat;
}
///********************************************************************


//    return Rcpp::List::create(
//                Rcpp::Named("incr") = incr,
//                Rcpp::Named("xsi") = xsi_new,
//            );

