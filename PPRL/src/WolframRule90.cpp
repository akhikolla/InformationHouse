// Enable C++11 via this plugin
// [[Rcpp::plugins("cpp11")]]

#include <Rcpp.h>
#include "WolframRule90.h"

using namespace Rcpp;
using namespace std;

//[[Rcpp::export]]
CharacterVector WolframRule90(CharacterVector ID, CharacterVector data, int lenBloom, int t){
  vector<string> CLKin =Rcpp::as<std::vector<string> > (data);
  vector<string> IDc =Rcpp::as<std::vector<string> > (ID);

  CLK* clkin = new CLK(lenBloom);
  CLK *clkout= new CLK(lenBloom);
  char *id = new char[lenBloom + 1];;
  char *str = new char[lenBloom + 1];;
 //
    CharacterVector CLKout;
 // //for all lines
  for (int i =0; i < data.size() ; i++){
    strcpy(id, IDc[i].c_str());
    strcpy(str, CLKin[i].c_str());
    clkin->copyFromString(id, str);
    for (int j = 0 ; j < t ; j++){
      WolframRule90c(clkin, clkout, t);
      clkin->copy(clkout);
    }
    clkout->copyToString(str, lenBloom);
    CLKout.push_back(string(str));
   }
   delete[] id;
   delete[] str;
   delete clkin;
   delete clkout;
   return CLKout;
}

void WolframRule90c(CLK* clkin, CLK* clkout, int t) {
  int left, right, length, l, r, m;
  length = clkin->getLength();
  left = length-1;
  right = 1;
  clkout->clear();
   //check lenght
  if (length != clkout->getLength()){
    Rcpp::Rcerr << "length problem "<< length << " " <<clkout->getLength()<< " \n";
  }

    // if (j>0)
    // clkin = clkout;
   for(int i =0; i < length ; i++){
      l=clkin->getBit(left);
      r=clkin->getBit(right);
      m=clkin->getBit(i);
     if (((l == 1 )&& (m==1)&&(r == 0))|| ((l == 1 )&&(m==0)&& (r == 0))||((l == 0 )&& (m==1)&&(r == 1))||((l == 0)&& (m==0)&&(r == 1))){
        clkout->setBit(i);
     }
      else{
        clkout->unsetBit(i);
     }
     left=i;
      right=(i+2)%length;

   }
 // }
}


