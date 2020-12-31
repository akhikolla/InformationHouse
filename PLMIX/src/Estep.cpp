#include <Rcpp.h>
using namespace Rcpp;
///' E-step in the EM algorithm for MAP estimation of a Bayesian mixture of Plackett-Luce models
///'
///' The function \code{Estep} updates the posterior component membership probabilities in the EM algorithm for MAP estimation of a Bayesian mixture of Plackett-Luce models.
///'
///' @param p Numeric \eqn{G}\eqn{\times}{x}\eqn{K} matrix of component-specific support parameters.
///' @param ref_order Numeric \eqn{G}\eqn{\times}{x}\eqn{K} matrix of component-specific reference orders.
///' @param weights Numeric vector of the \eqn{G} mixture weights.
///' @param pi_inv Numeric \eqn{N}\eqn{\times}{x}\eqn{K} data matrix of partial orderings.
///' @return Numeric \eqn{N}\eqn{\times}{x}\eqn{G} matrix of estimated posterior component membership probabilities.
// [[Rcpp::export]]
NumericMatrix Estep(NumericMatrix p, NumericMatrix ref_order, NumericVector weights, NumericMatrix pi_inv) {
int N = pi_inv.nrow();
int K = pi_inv.ncol();
int G = p.nrow();

NumericMatrix z(N,G);
 
int    s ;
int    j ;

int    slot ;
int    slot3 ;

int    tt ;

double f ;
double g ;

double temp_den ;

double h ;

for(s=0; s<N; s++){

h =0.0 ; 

for(j=0; j<G; j++){

f = 0.0 ;

slot  = 0 ;
tt    = 1 ; 
while(tt>-1 && slot<K){       
tt = pi_inv(s,slot)-1 ;
f = f + log(p(j,tt)) ; 
slot = slot+1 ; 
          if(slot<K){
              tt = pi_inv(s,slot)-1 ;
            }
}

g = 0.0 ;

slot  = 0 ;
tt    = 1 ;

/* FIRST compute the INITIAL denominator g */
/* i.e. sum of all support parameters  */

temp_den = 0.0;

for( slot3=0; slot3<K; slot3++){
temp_den = temp_den + p(j,slot3);
}

while(tt>-1 && slot<K){       
tt = pi_inv(s,slot)-1 ;

g = g + log(temp_den) ;

/* UPDATE the denominator temp_den removing the support for the current item */

temp_den = temp_den - p(j,tt) ;

if(temp_den<0.0){
}

slot = slot+1 ; 
          if(slot<K){
              tt = pi_inv(s,slot)-1 ;
            }

}

f = f - g ;

z(s,j) = weights[j]*exp(f) ; // numerator 
if(ISNAN(z(s,j))){
z(s,j) = 0.000000000001;
}
if(z(s,j)<=0.0){
z(s,j) = 0.000000000001;
}
h = h + z(s,j) ;
}

for(j=0; j<G; j++){
z(s,j) = z(s,j)/h ; // numerator divided by sum of numerators
if(ISNAN(z(s,j))){
z(s,j) = 0.000000000001;
}
if(z(s,j)<0.000000000000001){
z(s,j) = 0.0000000001;
}
}

/* recompute z(s,j) normalization (sum up to 1 for each s) */
h =0.0 ; 
for(j=0; j<G; j++){
h = h + z(s,j) ;
}

for(j=0; j<G; j++){
z(s,j)  = z(s,j)/h ;
}

}


return z;

}
