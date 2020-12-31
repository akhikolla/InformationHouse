#ifndef BICCLUST_H 
#define BICCLUST_H 

#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <Rmath.h>
#include <Rdefines.h>

#include <sys/time.h>
#include <sys/types.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <iostream>
#include <fstream>
#include <cstring>
#include <sstream>
#include <vector>
#include <math.h>
#include <algorithm>

using namespace std;
using namespace Rcpp;
using namespace arma;



class CritClust{
  
  //Attributs
  string crit;
  int k; 
  S4 m;                  //classe et modele
  NumericMatrix data;        //matrice des données  
  IntegerVector knownlabels; //les labels pour la AD
  bool DA;  
 public:
  //Constructeurs
  CritClust();
  CritClust(int k, S4 m, NumericMatrix data, string crit, IntegerVector knownlabels, bool DA);
  List ClustBestModel(vector<int> numExp);
  
};

#endif 





