/*
 * MTB_Tanimoto.h
 *
 *  Created on: 16.02.2018
 *      Author: Dorothea Rukasz
 */

#ifndef MTB_TANIMOTO_H_
#define MTB_TANIMOTO_H_

#include <vector>
#include <string>
#include <string.h>
#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <assert.h>
#include <iostream>
#include <Rcpp.h>
#include <cmath>
#include "MTB_Similarity.h"

using namespace std;


class MTB_Tanimoto: public MTB_Similarity {

public:
  MTB_Tanimoto() {
  }
  ;
  ~MTB_Tanimoto() {
  }
  ;

  string getName();
  double getRelativeValue(string o1, string o2);
  double getAbsoluteValue(string o1, string o2);

};

#endif /* MTB_TANIMOTO_H_ */
