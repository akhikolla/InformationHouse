#ifndef VECT_H
#define VECT_H

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


class Vect {
 
private:
  NumericMatrix Data;       //path to the data file 
  
  //bool is_readable(const string& file);

public:
  vector<int> experiments;   //variable numbers
 
  /**** Constructors ****/
  Vect();
  Vect(NumericMatrix Data,vector<int> experiments);
  Vect(NumericMatrix Data);

  //initialisation of the variables
  void initExperiments();                
  mat const_matrix(vector<int> vecteur);
  
  //calculation of the Bicreg 
  List bicReggen(vector<int> vectH, vector<int> vectY, int numr);
  //remove a part of a vector 
  vector<int> enlever_var(vector<int>& vecteur, vector<int>& varenv);
  //concatenate and order two variable vectors 
  vector<int> ajouter_var(vector<int>& vecteur, vector<int>& varajout);
  
};//class Vect
#endif
