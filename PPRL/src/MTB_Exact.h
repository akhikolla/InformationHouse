/*
 * MTB_Exact.h
 *
 *  Created on: 04.04.2017
 *      Author: schnell-42
 */

#ifndef MTB_EXACT_H_
#define MTB_EXACT_H_

#include <string>
#include <algorithm>
#include "MTB_Similarity.h"

using namespace std;


class MTB_ExactAlgorithm: public MTB_Similarity {

public:
	MTB_ExactAlgorithm() {
	}
	;
	~MTB_ExactAlgorithm() {
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

};

bool exact(string var1, string var2);

class MTB_ExactCapitalLettersAlgorithm: public MTB_Similarity {

public:
	MTB_ExactCapitalLettersAlgorithm() {
	}
	;
	~MTB_ExactCapitalLettersAlgorithm() {
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

};
bool exactCapitalLetters(string var1, string var2);
string StringToUpper(string var);

#endif /* MTB_EXACT_H_ */
