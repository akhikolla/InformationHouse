#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List process_snptext_cpp(DataFrame df, int snpmajor = 0, int maref = 0, int inhaplo = 0, int outallele = 1, int outhaplo = 0) {
  int nind;
  int nsnp;
	
  if(snpmajor){
    nind = df.length()/2;
    nsnp = df.nrows();
  }else{
    if(inhaplo){
      nind = df.nrows()/2;
      nsnp = df.length();
    }else{
      nind = df.nrows();
      nsnp = df.length()/2;
    }
  }
  
  if(outallele){	// output both alleles
    if(outhaplo){
      IntegerMatrix haplo_mat(nind*2, nsnp); // init haplotype matrix
      CharacterVector a1(nsnp), a2(nsnp); // init ref allele vectors
      
      // gather SNP information
      if(snpmajor){
        IntegerVector allele_vec(nsnp);	// number of alleles identified
        for(int i = 0; i < nind; i++){
          // alleles of individual i for all SNPs
          CharacterVector vec1 = as<CharacterVector>(df[2*i]);
          CharacterVector vec2 = as<CharacterVector>(df[2*i+1]);
          
          for(int j = 0; j < nsnp; j++){
            // first allele of an individual
            if(vec1[j] == "NA"){ 
              haplo_mat(i*2, j) = -9;
            }else{
              if(allele_vec[j] == 2){
                if(vec1[j] == a1[j]){
                  haplo_mat(i*2, j) = 1;
                }else if(vec1[j] == a2[j]){
                  haplo_mat(i*2, j) = 2;
                }
              }else if(allele_vec[j] == 1){
                if(vec1[j] == a1[j]){
                  haplo_mat(i*2, j) = 1;
                }else{
                  haplo_mat(i*2, j) = 2;
                  a2[j] = vec1[j];
                  allele_vec[j]++;
                }
              }else{
                haplo_mat(i*2, j) = 1;
                a1[j] = vec1[j];
                allele_vec[j]++;
              }
            }	// end first allele of an individual
            
            // second allele of an individual
            if(vec2[j] == "NA"){ 
              haplo_mat(i*2+1, j) = -9;
            }else{
              if(allele_vec[j] == 2){
                if(vec2[j] == a1[j]){
                  haplo_mat(i*2+1, j) = 1;
                }else if(vec2[j] == a2[j]){
                  haplo_mat(i*2+1, j) = 2;
                }
              }else if(allele_vec[j] == 1){
                if(vec2[j] == a1[j]){
                  haplo_mat(i*2+1, j) = 1;
                }else{
                  haplo_mat(i*2+1, j) = 2;
                  a2[j] = vec2[j];
                  allele_vec[j]++;
                }
              }else{
                haplo_mat(i*2+1, j) = 1;
                a1[j] = vec2[j];
                allele_vec[j]++;
              }
            }	// end second allele of an individual
          }
        }
      }else{
        if(inhaplo){
          for(int j = 0; j < nsnp; j++){
            // alleles of all individuals for SNP j
            CharacterVector vec = as<CharacterVector>(df[j]);
            int set_allele = 0;
            
            for(int i = 0; i < nind*2; i++){
              // first allele from an individual
              if(vec[i] == "NA"){ // first allele from an individual
                haplo_mat(i, j) = -9;
              }else{
                if(set_allele == 2){
                  if(vec[i] == a1[j]){
                    haplo_mat(i, j) = 1;
                  }else if(vec[i] == a2[j]){
                    haplo_mat(i, j) = 2;
                  }
                }else if(set_allele == 1){
                  if(vec[i] == a1[j]){
                    haplo_mat(i, j) = 1;
                  }else{
                    haplo_mat(i, j) = 2;
                    a2[j] = vec[i];
                    set_allele++;
                  }
                }else{
                  haplo_mat(i, j) = 1;
                  a1[j] = vec[i];
                  set_allele++;
                }
              } // end first allele from an individual
            } // end loop through individuals
          } // end loop through snps
        }else{
          for(int j = 0; j < nsnp; j++){
            // alleles of all individuals for SNP j
            CharacterVector vec1 = as<CharacterVector>(df[2*j]);
            CharacterVector vec2 = as<CharacterVector>(df[2*j+1]);
            int set_allele = 0;
            
            for(int i = 0; i < nind; i++){
              // first allele from an individual
              if(vec1[i] == "NA"){ // first allele from an individual
                haplo_mat(i*2, j) = -9;
              }else{
                if(set_allele == 2){
                  if(vec1[i] == a1[j]){
                    haplo_mat(i*2, j) = 1;
                  }else if(vec1[i] == a2[j]){
                    haplo_mat(i*2, j) = 2;
                  }
                }else if(set_allele == 1){
                  if(vec1[i] == a1[j]){
                    haplo_mat(i*2, j) = 1;
                  }else{
                    haplo_mat(i*2, j) = 2;
                    a2[j] = vec1[i];
                    set_allele++;
                  }
                }else{
                  haplo_mat(i*2, j) = 1;
                  a1[j] = vec1[i];
                  set_allele++;
                }
              } // end first allele from an individual
              
              // second allele from an individual
              if(vec2[i] == "NA"){ 
                haplo_mat(i*2+1, j) = -9;
              }else{
                if(set_allele == 2){
                  if(vec2[i] == a1[j]){
                    haplo_mat(i*2+1, j) = 1;
                  }else if(vec2[i] == a2[j]){
                    haplo_mat(i*2+1, j) = 2;
                  }
                }else if(set_allele == 1){
                  if(vec2[i] == a1[j]){
                    haplo_mat(i*2+1, j) = 1;
                  }else{
                    haplo_mat(i*2+1, j) = 2;
                    a2[j] = vec2[i];
                    set_allele++;
                  }
                }else{
                  haplo_mat(i*2+1, j) = 1;
                  a1[j] = vec2[i];
                  set_allele++;
                }
              }	// end second allele from an individual
            }
          }
        }
      }
      
      // minor allele as the reference allele, to be consistent with plink binary format
      for(int j = 0; j < nsnp; j++){
        // number of ref alleles, and total number of alleles
        int ra_count = 0, nm_count = 0;
        for(int i = 0; i < nind; i++){
          if(haplo_mat(i*2, j) >= 0){
            nm_count++;
            if(haplo_mat(i*2, j) == 1){
              ra_count++;
            }
          }
          
          if(haplo_mat(i*2+1, j) >= 0){
            nm_count++;
            if(haplo_mat(i*2+1, j) == 1){
              ra_count++;
            }
          }
        }
        
        if(ra_count == nm_count){
          a2[j] = "NA";
        }
        
        if(maref){	// make minor allele as the reference allele
          if(2*ra_count >= nm_count){
            String a3 = a1[j];
            a1[j] = a2[j];
            a2[j] = a3;
            for(int i = 0; i < nind; i++){
              if(haplo_mat(i*2, j) >= 0){
                haplo_mat(i*2, j) = 3 - haplo_mat(i*2, j);
              }
              if(haplo_mat(i*2+1, j) >= 0){
                haplo_mat(i*2+1, j) = 3 - haplo_mat(i*2+1, j);
              }
            }
          }
        }else{
          if(a1[j] == "2"){	
            if(a2[j] == "1" || a2[j] == "NA"){
              a1[j] = a2[j];
              a2[j] = "2";
              for(int i = 0; i < nind; i++){
                if(haplo_mat(i*2, j) >= 0){
                  haplo_mat(i*2, j) = 3 - haplo_mat(i*2, j);
                }
                if(haplo_mat(i*2+1, j) >= 0){
                  haplo_mat(i*2+1, j) = 3 - haplo_mat(i*2+1, j);
                }
              }  
            }   
          }
        }
      }
      
      List output = List::create(Named("data") = haplo_mat, Named("alleles") = DataFrame::create(Named("allele_1") = a1, Named("allele_2") = a2));
      return output;
    }else{
      IntegerMatrix haplo_mat(nind, nsnp*2); // init genotype matrix
      CharacterVector a1(nsnp), a2(nsnp); // init ref allele vectors
      
      // gather SNP information
      if(snpmajor){
        IntegerVector allele_vec(nsnp);	// number of alleles identified
        for(int i = 0; i < nind; i++){
          // alleles of individual i for all SNPs
          CharacterVector vec1 = as<CharacterVector>(df[2*i]);
          CharacterVector vec2 = as<CharacterVector>(df[2*i+1]);
          
          for(int j = 0; j < nsnp; j++){
            // first allele of an individual
            if(vec1[j] == "NA"){ 
              haplo_mat(i, j*2) = -9;
            }else{
              if(allele_vec[j] == 2){
                if(vec1[j] == a1[j]){
                  haplo_mat(i, j*2) = 1;
                }else if(vec1[j] == a2[j]){
                  haplo_mat(i, j*2) = 2;
                }
              }else if(allele_vec[j] == 1){
                if(vec1[j] == a1[j]){
                  haplo_mat(i, j*2) = 1;
                }else{
                  haplo_mat(i, j*2) = 2;
                  a2[j] = vec1[j];
                  allele_vec[j]++;
                }
              }else{
                haplo_mat(i, j*2) = 1;
                a1[j] = vec1[j];
                allele_vec[j]++;
              }
            }	// end first allele of an individual
            
            // second allele of an individual
            if(vec2[j] == "NA"){ 
              haplo_mat(i, j*2+1) = -9;
            }else{
              if(allele_vec[j] == 2){
                if(vec2[j] == a1[j]){
                  haplo_mat(i, j*2+1) = 1;
                }else if(vec2[j] == a2[j]){
                  haplo_mat(i, j*2+1) = 2;
                }
              }else if(allele_vec[j] == 1){
                if(vec2[j] == a1[j]){
                  haplo_mat(i, j*2+1) = 1;
                }else{
                  haplo_mat(i, j*2+1) = 2;
                  a2[j] = vec2[j];
                  allele_vec[j]++;
                }
              }else{
                haplo_mat(i, j*2+1) = 1;
                a1[j] = vec2[j];
                allele_vec[j]++;
              }
            }	// end second allele of an individual
          }
        }
      }else{
        if(inhaplo){
          for(int j = 0; j < nsnp; j++){
            // alleles of all individuals for SNP j
            CharacterVector vec = as<CharacterVector>(df[j]);
            int set_allele = 0;
            
            for(int i = 0; i < nind; i++){
              // first allele from an individual
              if(vec[i*2] == "NA"){ // first allele from an individual
                haplo_mat(i, j*2) = -9;
              }else{
                if(set_allele == 2){
                  if(vec[i*2] == a1[j]){
                    haplo_mat(i, j*2) = 1;
                  }else if(vec[i*2] == a2[j]){
                    haplo_mat(i, j*2) = 2;
                  }
                }else if(set_allele == 1){
                  if(vec[i*2] == a1[j]){
                    haplo_mat(i, j*2) = 1;
                  }else{
                    haplo_mat(i, j*2) = 2;
                    a2[j] = vec[i*2];
                    set_allele++;
                  }
                }else{
                  haplo_mat(i, j*2) = 1;
                  a1[j] = vec[i*2];
                  set_allele++;
                }
              } // end first allele from an individual
              
              // second allele from an individual
              if(vec[i*2+1] == "NA"){ // first allele from an individual
                haplo_mat(i, j*2+1) = -9;
              }else{
                if(set_allele == 2){
                  if(vec[i*2+1] == a1[j]){
                    haplo_mat(i, j*2+1) = 1;
                  }else if(vec[i*2+1] == a2[j]){
                    haplo_mat(i, j*2+1) = 2;
                  }
                }else if(set_allele == 1){
                  if(vec[i*2+1] == a1[j]){
                    haplo_mat(i, j*2+1) = 1;
                  }else{
                    haplo_mat(i, j*2+1) = 2;
                    a2[j] = vec[i*2+1];
                    set_allele++;
                  }
                }else{
                  haplo_mat(i, j*2+1) = 1;
                  a1[j] = vec[i*2+1];
                  set_allele++;
                }
              } // end second allele from an individual
            } // end loop through individuals
          } // end loop through snps
        }else{
          for(int j = 0; j < nsnp; j++){
            // alleles of all individuals for SNP j
            CharacterVector vec1 = as<CharacterVector>(df[2*j]);
            CharacterVector vec2 = as<CharacterVector>(df[2*j+1]);
            int set_allele = 0;
            
            for(int i = 0; i < nind; i++){
              // first allele from an individual
              if(vec1[i] == "NA"){ // first allele from an individual
                haplo_mat(i, j*2) = -9;
              }else{
                if(set_allele == 2){
                  if(vec1[i] == a1[j]){
                    haplo_mat(i, j*2) = 1;
                  }else if(vec1[i] == a2[j]){
                    haplo_mat(i, j*2) = 2;
                  }
                }else if(set_allele == 1){
                  if(vec1[i] == a1[j]){
                    haplo_mat(i, j*2) = 1;
                  }else{
                    haplo_mat(i, j*2) = 2;
                    a2[j] = vec1[i];
                    set_allele++;
                  }
                }else{
                  haplo_mat(i, j*2) = 1;
                  a1[j] = vec1[i];
                  set_allele++;
                }
              } // end first allele from an individual
              
              // second allele from an individual
              if(vec2[i] == "NA"){ 
                haplo_mat(i, j*2+1) = -9;
              }else{
                if(set_allele == 2){
                  if(vec2[i] == a1[j]){
                    haplo_mat(i, j*2+1) = 1;
                  }else if(vec2[i] == a2[j]){
                    haplo_mat(i, j*2+1) = 2;
                  }
                }else if(set_allele == 1){
                  if(vec2[i] == a1[j]){
                    haplo_mat(i, j*2+1) = 1;
                  }else{
                    haplo_mat(i, j*2+1) = 2;
                    a2[j] = vec2[i];
                    set_allele++;
                  }
                }else{
                  haplo_mat(i, j*2+1) = 1;
                  a1[j] = vec2[i];
                  set_allele++;
                }
              }	// end second allele from an individual
            }
          }
        }
      }
      
      // minor allele as the reference allele, to be consistent with plink binary format
      for(int j = 0; j < nsnp; j++){
        // number of ref alleles, and total number of alleles
        int ra_count = 0, nm_count = 0;
        for(int i = 0; i < nind; i++){
          if(haplo_mat(i, j*2) >= 0){
            nm_count++;
            if(haplo_mat(i, j*2) == 1){
              ra_count++;
            }
          }
          
          if(haplo_mat(i, j*2+1) >= 0){
            nm_count++;
            if(haplo_mat(i, j*2+1) == 1){
              ra_count++;
            }
          }
        }
        
        if(ra_count == nm_count){
          a2[j] = "NA";
        }
        
        if(maref){	// make minor allele as the reference allele
          if(2*ra_count >= nm_count){
            String a3 = a1[j];
            a1[j] = a2[j];
            a2[j] = a3;
            for(int i = 0; i < nind; i++){
              if(haplo_mat(i, j*2) >= 0){
                haplo_mat(i, j*2) = 3 - haplo_mat(i, j*2);
              }
              if(haplo_mat(i, j*2+1) >= 0){
                haplo_mat(i, j*2+1) = 3 - haplo_mat(i, j*2+1);
              }
            }
          }
        }else{
          if(a1[j] == "2"){	
            if(a2[j] == "1" || a2[j] == "NA"){
              a1[j] = a2[j];
              a2[j] = "2";
              for(int i = 0; i < nind; i++){
                if(haplo_mat(i, j*2) >= 0){
                  haplo_mat(i, j*2) = 3 - haplo_mat(i, j*2);
                }
                if(haplo_mat(i, j*2+1) >= 0){
                  haplo_mat(i, j*2+1) = 3 - haplo_mat(i, j*2+1);
                }
              }  
            }   
          }
        }
      }
      
      List output = List::create(Named("data") = haplo_mat, Named("alleles") = DataFrame::create(Named("allele_1") = a1, Named("allele_2") = a2));
      return output;
    }
  }else{	// output genotype coded as 0, 1 or 2
    IntegerMatrix geno_mat(nind, nsnp); // init genotype matrix
    CharacterVector a1(nsnp), a2(nsnp); // init ref allele vectors
    
    // gather SNP information
    if(snpmajor){
      IntegerVector allele_vec(nsnp);
      for(int i = 0; i < nind; i++){
        CharacterVector vec1 = as<CharacterVector>(df[2*i]);
        CharacterVector vec2 = as<CharacterVector>(df[2*i+1]);
        
        for(int j = 0; j < nsnp; j++){
          if(vec1[j] == "NA" || vec2[j] == "NA"){
            geno_mat(i, j) = -9;
          }else{
            if(allele_vec[j] == 2){
              if(vec1[j] != vec2[j]){
                geno_mat(i, j) = 1;
              }else{
                if(vec1[j] == a1[j]){
                  geno_mat(i, j) = 2;
                }else if(vec1[j] == a2[j]){
                  geno_mat(i, j) = 0;
                }else{
									geno_mat(i, j) = 9;
                }
              }
            }else if(allele_vec[j] == 1){
              if(vec1[j] != vec2[j]){
                geno_mat(i, j) = 1;
                if(vec1[j] == a1[j]){
                  a2[j] = vec2[j];
                }else{
                  a2[j] = vec1[j];
                }
                allele_vec[j]++;
              }else{
                if(vec1[j] == a1[j]){
                  geno_mat(i, j) = 2;
                }else{
                  geno_mat(i, j) = 0;
                  a2[j] = vec1[j];
                  allele_vec[j]++;
                }
              }
            }else{
              if(vec1[j] != vec2[j]){
                geno_mat(i, j) = 1;
                a1[j] = vec1[j];
                a2[j] = vec2[j];
                allele_vec[j]+=2;
              }else{
                geno_mat(i, j) = 2;
                a1[j] = vec1[j];
                allele_vec[j]++;
              }
            }
          }
        }
      }
    }else{
      if(inhaplo){
        for(int j = 0; j < nsnp; j++){
          CharacterVector vec = as<CharacterVector>(df[j]);
          int set_allele = 0;
          
          for(int i = 0; i < nind; i++){
            if(vec[i*2] == "NA" || vec[i*2+1] == "NA"){
              geno_mat(i, j) = -9;
            }else{
              if(set_allele == 2){
                if(vec[i*2] != vec[i*2+1]){
                  geno_mat(i, j) = 1;
                }else{
                  if(vec[i*2] == a1[j]){
                    geno_mat(i, j) = 2;
                  }else if(vec[i*2] == a2[j]){
                    geno_mat(i, j) = 0;
                  }else{
                    geno_mat(i, j) = 9;
                  }
                }
              }else if(set_allele == 1){
                if(vec[i*2] != vec[i*2+1]){
                  geno_mat(i, j) = 1;
                  if(vec[i*2] == a1[j]){
                    a2[j] = vec[i*2+1];
                  }else{
                    a2[j] = vec[i*2];
                  }
                  set_allele++;
                }else{
                  if(vec[i*2] == a1[j]){
                    geno_mat(i, j) = 2;
                  }else{
                    geno_mat(i, j) = 0;
                    a2[j] = vec[i*2];
                    set_allele++;
                  }
                }
              }else{
                if(vec[i*2] != vec[i*2+1]){
                  geno_mat(i, j) = 1;
                  a1[j] = vec[i*2];
                  a2[j] = vec[i*2+1];
                  set_allele+=2;
                }else{
                  geno_mat(i, j) = 2;
                  a1[j] = vec[i*2];
                  set_allele++;
                }
              }
            }
          }
        }
      }else{
        for(int j = 0; j < nsnp; j++){
          CharacterVector vec1 = as<CharacterVector>(df[2*j]);
          CharacterVector vec2 = as<CharacterVector>(df[2*j+1]);
          int set_allele = 0;
          
          for(int i = 0; i < nind; i++){
            if(vec1[i] == "NA" || vec2[i] == "NA"){
              geno_mat(i, j) = -9;
            }else{
              if(set_allele == 2){
                if(vec1[i] != vec2[i]){
                  geno_mat(i, j) = 1;
                }else{
                  if(vec1[i] == a1[j]){
                    geno_mat(i, j) = 2;
                  }else if(vec1[i] == a2[j]){
                    geno_mat(i, j) = 0;
                  }else{
                    geno_mat(i, j) = 9;
                  }
                }
              }else if(set_allele == 1){
                if(vec1[i] != vec2[i]){
                  geno_mat(i, j) = 1;
                  if(vec1[i] == a1[j]){
                    a2[j] = vec2[i];
                  }else{
                    a2[j] = vec1[i];
                  }
                  set_allele++;
                }else{
                  if(vec1[i] == a1[j]){
                    geno_mat(i, j) = 2;
                  }else{
                    geno_mat(i, j) = 0;
                    a2[j] = vec1[i];
                    set_allele++;
                  }
                }
              }else{
                if(vec1[i] != vec2[i]){
                  geno_mat(i, j) = 1;
                  a1[j] = vec1[i];
                  a2[j] = vec2[i];
                  set_allele+=2;
                }else{
                  geno_mat(i, j) = 2;
                  a1[j] = vec1[i];
                  set_allele++;
                }
              }
            }
          }
        }
      }
    }
    
    // minor allele as the reference allele, to be consistent with plink binary format
    for(int j = 0; j < nsnp; j++){
      int ra_count = 0, nm_count = 0;
      for(int i = 0; i < nind; i++){
        if(geno_mat(i, j) >= 0){
          nm_count++;
          ra_count+= geno_mat(i, j);
        }
      }
      
      if(ra_count == (2*nm_count)){
        a2[j] = "NA";
      }
      
      if(maref){
        if(ra_count >= nm_count){
          String a3 = a1[j];
          a1[j] = a2[j];
          a2[j] = a3;
          for(int i = 0; i < nind; i++){
            if(geno_mat(i, j) >= 0){
              geno_mat(i, j) = 2 - geno_mat(i, j);
            }
          }
        }
      }else{
        if(a1[j] == "2"){
          if(a2[j] == "1" || a2[j] == "NA"){
            a1[j] = a2[j];
            a2[j] = "2";
            for(int i = 0; i < nind; i++){
              if(geno_mat(i, j) >= 0){
                geno_mat(i, j) = 2 - geno_mat(i, j);
              }
            }
          }
        }
      }
    }
    
    List output = List::create(Named("data") = geno_mat, Named("alleles") = DataFrame::create(Named("allele_1") = a1, Named("allele_2") = a2));
    return output;
  }
}


