#include <Rcpp.h>
#include <fstream>
#include <iostream>
using namespace Rcpp;

// [[Rcpp::export]]
IntegerMatrix read_bed_cpp(std::string file, int nind, int nsnp){
  // first three bytes of plink bed file are fixed (SNP-major format)
  IntegerMatrix first3bytes(3, 8);
  first3bytes.row(0) = IntegerVector::create(0, 1, 1, 0, 1, 1, 0, 0);
  first3bytes.row(1) = IntegerVector::create(0, 0 ,0, 1, 1, 0, 1, 1);
  first3bytes.row(2) = IntegerVector::create(0, 0 ,0, 0, 0, 0, 0, 1);
  
  // initialize output matrix of genotypes  
  IntegerMatrix geno_mat(nind, nsnp);
  
  std::ifstream f(file.c_str(), std::ios::binary | std::ios::in);

  char c;
  int bytecount = 0;
  int snpcount = 0;
  int indcount = 0;
  
  while(f.get(c)){
    if(bytecount < 3){
      for(int i = 7; i >= 0; i--){
        if(first3bytes(bytecount, 7-i) != ((c >> i) & 1)){
          f.close();
          Rcerr << "Check input data.\n";
          return 0;
        }
      }
      bytecount++;
      continue;
    }
    
    for(int i = 0; i < 4; i++){
      int geno;
      if(((c >> 2*i) & 1) == 0){
        if(((c >> (2*i+1)) & 1) == 0){
          geno = 2;
        }else{
          geno = 1;
        }
      }else{
        if(((c >> (2*i+1)) & 1) == 0){
          geno = -9;
        }else{
          geno = 0;
        }
      }
      
      geno_mat(indcount, snpcount) = geno;
      indcount++;
      if(indcount == nind){
        indcount = 0;
        snpcount++;
        break;
      }
    }
  }
  f.close();
  
  return geno_mat;
}
