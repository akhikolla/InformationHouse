#include <math.h>
#include <iostream>
#include <Rcpp.h>

using namespace Rcpp;

unsigned long zunif(){
  double r = floor(runif(1, 0, 4294967296.0))[0];
  return (unsigned long)r;
}

// [[Rcpp::export(".breed")]]
IntegerVector breed(IntegerVector Parents, int ns, int Ns, int nLoci){
  int s;
  int subpopStart;
  int i;
  int nNumAlleles = ns * Ns * 2 * nLoci;
  
  IntegerVector Children(nNumAlleles);
  
  for(s = 0; s < ns; s++){
    subpopStart = s * Ns * 2 * nLoci;
    
    //i2 = subpopStart + *Ns - 1;
    
    int p1, p2;
    IntegerVector::const_iterator iParent1, iParent2;
    IntegerVector::iterator iChild;
    
    int ctr = 0;
    unsigned int g, u = zunif();
    int loc;
    NumericVector z1 = floor(runif(Ns, 0, Ns));
    NumericVector z2 = floor(runif(Ns, 0, Ns));
    
    for(i = 0; i < Ns; i++){
      p1 = subpopStart + 2 * nLoci * (int)z1[i];
      p2 = subpopStart + 2 * nLoci * (int)z2[i];
      
      iChild = Children.begin() + subpopStart + (2 *  nLoci * i);
      iParent1 = Parents.begin() + p1;
      iParent2 = Parents.begin() + p2;
      
      for(loc = 0; loc < nLoci; loc++){
        g = u & 1; // this pops a bit off u
        u >>= 1; // and this shifts to the right so the bit is removed
        
        
        // g will be 0 or 1
        int a1  = iParent1[2 * loc + g];
        
        g = u & 1; // this pops a bit off u
        u >>= 1; // and this shifts to the right so the bit is removed
        ctr += 2;
        
        int a2  = iParent2[2 * loc + g];
        
        if(a1 > a2){
          int nSwap = a1;
          a1 = a2;
          a2 = nSwap;
        }
        
        iChild[2 * loc] = a1;
        iChild[2 * loc + 1] = a2;
        
        if(ctr == 14){ // g and u should be 32 bits long but I don't trust them
          u = zunif();
          ctr = 0;
        }
      }
    }
  }
  
  return Children;
}
