#include <Rcpp.h>
using namespace Rcpp;

// Helper functions

// [[Rcpp::export]]
NumericMatrix ansol(IntegerMatrix awm, IntegerVector maxsc) {
// awm = response matrix
// maxsc = max score possible

int npers = awm.nrow();
int nitem = awm.ncol();

NumericMatrix pperg(npers,2);

for(int pe = 0; pe < npers; pe++)
  {
  
  int sumresp = 0; // sum of responses
  int summax  = 0; // sum of maximal scores
  
  for(int it = 0; it < nitem; it++)
  {
    
  if(IntegerVector::is_na(awm(pe,it)))
    {
    continue;
    } else 
        {
        sumresp += awm(pe,it);
        summax += maxsc(it);
        }
    
  }
  
  // wenn max oder 0
  
  if(sumresp == summax)
    {
    pperg(pe,0)  = R_PosInf;
    pperg(pe,1)  = NA_REAL;
      
    } else if(sumresp == 0)
      {
      pperg(pe,0)  =  R_NegInf;
      pperg(pe,1)  = NA_REAL;  
        
      } else 
        {
        continue;
        }

  }

return pperg;
}




