#include <Rcpp.h>
using namespace Rcpp;

//' @rdname loglikelihood
//' @export
// [[Rcpp::export]]
double loglikPLMIX(NumericMatrix p, NumericMatrix ref_order, NumericVector weights, NumericMatrix pi_inv) {

int N = pi_inv.nrow();
int K = pi_inv.ncol();
int G = p.nrow();

int    s ;
int    j ;

int    slot ;
int    slot3 ;

int    tt ;

double f  ;
double g  ;
double h  ;
double ll ;

ll = 0.0 ;


for(s=0; s<N; s++){

h = 0.0 ;

for(j=0; j<G; j++){

f = 0.0 ;

slot  = 0 ;
tt    = 1 ; 


/* FIRST compute the INITIAL denominator g */
/* i.e. sum of all support parameters  */

g = 0.0;

for( slot3=0; slot3<K; slot3++){
g = g + p(j,slot3);
}

while(tt>-1 && slot<K){       

tt = pi_inv(s,slot)-1 ;
f = f + (log(p(j,tt)) - log(g)) ;

/* UPDATE the denominator g removing the support for the current item */

g = g - p(j,tt);

if(g<0){
/* printf(" SOMETHING IS WRONG WITH THE DENOMINATOR !!! \\n"); */
}

slot = slot+1 ; 
          if(slot<K){
              tt = pi_inv(s,slot)-1 ;
            }
}

h = h + weights[j]*exp(f) ;

}

ll = ll + log(h) ;

}

return ll;

}
