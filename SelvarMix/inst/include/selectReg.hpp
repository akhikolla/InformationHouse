#ifndef SELECTREG_H 
#define SELECTREG_H 

#include "vect.hpp"

class SelectReg{

  /**** Attribute ****/
  Vect v;    


public:
  /**** Constructors ****/
  SelectReg();  
  SelectReg(Vect v); 

  /**** Methods ****/
  
  //Exclusion step
  void exclusion_reg(vector<int>& varSelectReg, vector<int>& varNonSig,vector<int>& jE, vector<int>& jI, int& stop, int& InitialProjectsNb);

  //Inclusion step
  void inclusion_reg(vector<int> varSelect, vector<int>& varSelectReg, vector<int>& varNonSig,vector<int>& jE, vector<int>& jI, int& arret, int& InitialProjectsNb);

  //backward stepwise algorithm for the regression
  vector<int> selectReg(vector<int> varSelect,vector<int>& varNonSig, int& InitialProjectsNb);

};

#endif
