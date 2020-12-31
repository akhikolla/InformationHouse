#ifndef SELECT_H 
#define SELECT_H 
#include "vect.hpp" 
#include "critClust.hpp"
#include "selectReg.hpp"

class Select{
  
  Vect v;
  CritClust b;
  SelectReg sReg;
  int packSize;
public:

  /**** Constructors ****/
  Select();
  Select(Vect v, CritClust b, SelectReg sReg, int packSize);
  Select(Vect v, SelectReg sReg, int packSize);

  /**** Methods ****/
  List selectS(vector<int> Order);
  vector<int> selectW(vector<int> Order, vector<int> OtherVar);
  
};

#endif
