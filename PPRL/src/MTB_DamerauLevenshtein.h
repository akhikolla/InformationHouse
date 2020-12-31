/*
 * MTB_DamerauLevenshtein.h
 *
 *  Created on: 21.02.2017
 *      Author: schnell-42
 */

#ifndef MTB_DAMERAULEVENSHTEIN_H_
#define MTB_DAMERAULEVENSHTEIN_H_

#include <iostream>
#include <cstring>
#include <algorithm>
#include "MTB_Similarity.h"
#include "MTB_Levenshtein.h"

using namespace std;



class MTB_DamerauLevenshteinAlgorithm: public MTB_Similarity {

public:
	MTB_DamerauLevenshteinAlgorithm() {
	}
	;
	~MTB_DamerauLevenshteinAlgorithm() {
	}
	;


	 string getName();
	 double getRelativeValue(string o1, string o2);
	 double getAbsoluteValue(string o1, string o2);
//	 double getThreshold();
//	 bool isAbsolute();
//	 bool isSimilar(string o1, string o2);
//	 void setAbsolute(bool absolute);
//	 void setThreshold(double threshold);
//	 void Free(){ delete this; }
};




int damlevdist( string a,  string b);


#endif /* MTB_DAMERAULEVENSHTEIN_H_ */
