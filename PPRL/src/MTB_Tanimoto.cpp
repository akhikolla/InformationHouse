/*
 * MTB_Tanimoto.cpp
 *
 *  Created on: 16.02.2018
 *      Author: Dorothea Rukasz
 */

#include "MTB_Tanimoto.h"


string MTB_Tanimoto::getName() {
  return "Taminoto";
}

double MTB_Tanimoto::getRelativeValue(string var1, string var2) {
  // if (var1.length() != var2.length()) {
  //   Rcpp::Rcout <<
  //     "Bloomfilter have different length. The tanimoto distance will be calculated up to the minimum length of both bloomfilters."<< endl;
  // }
  int minStringLength = var1.length() < var2.length() ? var1.length() : var2.length();
  int identicalCharacters = 0;
  int differentCharacters = abs((int)(var1.length() - var2.length()));
  for (int stringPos = 0; stringPos < minStringLength; stringPos++) {
    if ((var1[stringPos] == '1') && (var2[stringPos] == '1')) {
      identicalCharacters++;
    } else if (((var1[stringPos] == '0') && (var2[stringPos] == '1'))
                 || (((var1[stringPos] == '1') && (var2[stringPos] == '0')))) {
      differentCharacters++;
    }
  }
  return (double) (identicalCharacters / ( identicalCharacters + differentCharacters));

}

double MTB_Tanimoto::getAbsoluteValue(string var1, string var2) {

  return getRelativeValue(var1, var2);
}




