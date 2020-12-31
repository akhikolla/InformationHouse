#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
NumericMatrix matVecProd(NumericMatrix m, NumericVector v){
    mat M = as<mat>(m);
    vec V = as<vec>(v);
    M.each_col() %= V;
    return wrap(M);
}

// [[Rcpp::export]]
NumericMatrix matVecProdSum(NumericMatrix m, NumericVector ext, NumericVector v, NumericVector group){
    mat M = as<mat>(m);
    vec Ext = as<vec>(ext);
    vec V = as<vec>(v);
    ivec G = as<ivec>(group);
    for(uword i=0; i<Ext.n_elem; i++) M.col(M.n_cols-Ext.n_elem+i).fill(Ext[i]);
    M.each_col() %= V;
    if(G.is_empty()) return wrap(M);
    uword nGroup = G.n_elem - 1;
    mat z = zeros<mat>(nGroup, M.n_cols);
    for(uword i=0; i< nGroup; i++) z.row(i) = sum(M.rows(G[i], G[i+1]-1), 0);
    return wrap(z);
}

// [[Rcpp::export]]
NumericMatrix matVecProdSumExt(NumericMatrix m, NumericVector ext, NumericVector ext2, NumericVector v, NumericVector group){
    mat M = as<mat>(m);
    vec Ext = as<vec>(ext);
    vec Ext2 = as<vec>(ext2);
    vec V = as<vec>(v);
    ivec G = as<ivec>(group);
    for(uword i=0; i<Ext.n_elem; i++) M.col(M.n_cols-1-Ext.n_elem+i).fill(Ext[i]);
    M.col(M.n_cols-1) = Ext2;
    M.each_col() %= V;
    if(G.is_empty()) return wrap(M);
    uword nGroup = G.n_elem - 1;
    mat z = zeros<mat>(nGroup, M.n_cols);
    for(uword i=0; i< nGroup; i++) z.row(i) = sum(M.rows(G[i], G[i+1]-1), 0);
    return wrap(z);
}

// [[Rcpp::export]]
NumericMatrix groupProd(NumericVector v, NumericVector group){
    vec V = as<vec>(v);
    ivec G = as<ivec>(group);
    uword nGroup = G.n_elem - 1;
    vec z = zeros<vec> (nGroup);
    for(uword i=0; i< nGroup; i++) z[i] = prod(V.subvec(G[i], G[i+1]-1));
    return wrap(z);
}


