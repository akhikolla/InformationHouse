/*
 * MTB_JaroWinklerMcLaughlinWinklerLynch.h
 *
 *  Created on: 10.04.2017
 *      Author: schnell-42
 */

#ifndef MTB_JAROWINKLERMCLAUGHLINWINKLERLYNCH_H_
#define MTB_JAROWINKLERMCLAUGHLINWINKLERLYNCH_H_

#include <string.h>
#include <ctype.h>
#include <iostream>
#include <cstring>
#include <vector>
#include "MTB_Similarity.h"

using namespace std;

double strcmp95(char *ying, char *yang, long y_length, bool *ind_c);
double JWMcLWL(string word1, string word2, bool *ind_c);




class MTB_JaroWinklerMcLaughlinWinklerLynchAlgorithm: public MTB_Similarity {

public:
	MTB_JaroWinklerMcLaughlinWinklerLynchAlgorithm() {
	}
	;
	~MTB_JaroWinklerMcLaughlinWinklerLynchAlgorithm() {
	}
	;

	//
	string getName();
	double getRelativeValue(string o1, string o2);
	double getAbsoluteValue(string o1, string o2);

};

#endif /* MTB_JAROWINKLERMCLAUGHLINWINKLERLYNCH_H_ */
