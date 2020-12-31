//============================================================================
// Name        : ElegantPairingFunctions.cpp
// Author      : Rukasz
// Version     :
// Copyright   : Your copyright notice
// Description : Ansi-style
//============================================================================

#include "ElegantPairingFunctions.h"

unsigned unorderedPairing(int x, int y){
  if (x<0||y <0){
    Rcpp::Rcerr << "Input has to be non-negative."<< endl;
    return -1;
  }
	 return x*y+floor(pow ( (abs(x-y)-1),2)/4);
}


using namespace std;
