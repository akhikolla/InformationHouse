# rres
Realized Relatedness Estimation and Simulation

This package contains useful functions for studying realized genetic relatedness between individuals who can trace back to some recent common ancestors. Functions in this package fall into four main categories:

1. Simulation of inheritance pattern and SNP marker data    
-sim.recomb()   
-sim.haplotype()   
-populate.snp()

2. Scoring identity by descent information from simulated inheritance patterns   
-ibd.length()   
-ibd.proportion()  
-ibd.segment()  
-ibd.marker()

3. Estimation of realized relatedness from SNP marker data using the GRM estimators   
-grm.pair()  
-grm.matrix()  
-ld.weights()

4. Various utility functions   
-check.pedinfo()  
-get.pedindex()  
-fgl2ibd()  
-fgl2relatedness()  
-recode.ibd()  
-recode.snpdata()  
-read.plink.text()  
-read.plink.binary()  
-write.ibdhaplo()

