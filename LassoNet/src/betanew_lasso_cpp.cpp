#include <Rcpp.h>


using namespace Rcpp;



// [[Rcpp::export]]
List betanew_lasso_cpp(NumericMatrix xx, NumericVector xy, NumericVector beta, NumericMatrix M,NumericVector y, double Lambda1, double Lambda2, double iter, double tol){




int p          = beta.size();
int N          = y.size();

double stop = 0;
double convsteps = 0;
NumericVector betanew(p);
NumericVector betaprev(p);
double abssum;



while  (stop < iter) {
    for (int j = 0; j <= p-1; j++) {
        betaprev[j] =  betanew[j]  ;
    }
    
    
    for (int i = 0; i < p; i++) {
       
        double xxprodb    = 0;
        double num        = 0;
        double den        = 0;
        double Mjjprodbjj = 0;
        
       
        NumericVector xtemp     = xx(_, i ) ;
        
        for (int j = 0; j <= p-1; j++) {
            xxprodb +=  (xtemp[j])*(betanew[j])   ;
        }
        
        NumericVector Mjj = M(i,_);
        
        
        for (int j = 0; j <= p-1; j++) {
            if(j==i)
            {
                
            }
            else
            {
                Mjjprodbjj += (Mjj[j])*(betanew[j]) ;
            }
        }
        
        num = 2*(N*betanew[i] + xy[i] - xxprodb) - 2*Lambda2*Mjjprodbjj  ;
        den = 2*N + 2*Lambda2*M(i, i) ;
        
        
        if (num > Lambda1)
        {
            betanew[i] = (num - Lambda1) / den ;
            
        }
        
        else if(-num > Lambda1)
        {
            betanew[i] = (num + Lambda1) / den ;
        }
        else
        {
            betanew[i] = 0 ;
        }
    }
    // check convergence
    abssum = 0;
    
    // L1 convergence
    for (int j = 0; j <= p-1; j++) {
        abssum += std::abs(betanew[j]-betaprev[j]);
    }
    
    abssum  = abssum / N;
    

    
    stop = stop + 1;
    
    if (abssum  < tol)
    {
        stop = iter;
    }
    convsteps = convsteps + 1;
    
}
    
return  Rcpp::List::create(Rcpp::Named("beta") = betanew ,Rcpp::Named("steps") = convsteps);
}
