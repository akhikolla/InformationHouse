/*
 * MTB_JaccardAlgorithm.h
 *
 *  Created on: 05.07.2017
 *      Author: schnell-42
 */

#ifndef MTB_JACCARDALGORITHM_H_
#define MTB_JACCARDALGORITHM_H_

#include <algorithm>
#include <string>
#include <ctype.h>
#include <iostream>
#include <cstring>
#include <vector>
#include "MTB_Similarity.h"

class MTB_Jaccard: public MTB_Similarity {

public:
	MTB_Jaccard() {
	}
	;
	~MTB_Jaccard() {
	}
	;

	//
	string getName();
	double getRelativeValue(string o1, string o2);
	double getAbsoluteValue(string o1, string o2);
//	 double getThreshold();
//	 bool isAbsolute();
//	 bool isSimilar(string o1, string o2);
//	 void setAbsolute(bool absolute);
//	 void setThreshold(double threshold);


};



#endif /* MTB_JACCARDALGORITHM_H_ */
