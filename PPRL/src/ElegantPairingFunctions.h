/*
 * ElegantPairingFunctions.h
 *
 *  Created on: 17.06.2016
 *      Author: schnell-42
 */

#ifndef ELEGANTPAIRINGFUNCTIONS_H_
#define ELEGANTPAIRINGFUNCTIONS_H_
#include <Rcpp.h>
#include <stdlib.h>
#include <math.h>
using namespace std;


/**
* Unordered Pairing Function:
*
* elegantPairing(x,y)=x*y+floor( ( (abs(x-y)-1)**2)/4)
*
* elegantPairing(x,y) = elegantPairing(y,x)
*
*
 When x and y are non−negative integers, elegantPairing(x,y) outputs a single
 non−negative integer that is uniquely associated with that unorderd pair.
 */
unsigned unorderedPairing(int x, int y);

#endif /* FUNCTIONS_ELEGANTPAIRINGFUNCTIONS_H_ */
