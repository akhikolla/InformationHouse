/*
 * MTB_DoubleMetaphone.h
 *
 *  Created on: 06.07.2017
 *      Author: schnell-42
 */

#ifndef MTB_DOUBLEMETAPHONE_H_
#define MTB_DOUBLEMETAPHONE_H_

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
#include "MTB_Similarity.h"

using namespace std;


class MTB_DoubleMetaphone: public MTB_Similarity {

public:
	MTB_DoubleMetaphone() {
	}
	;
	~MTB_DoubleMetaphone() {
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

void doubleMetaphone(const string &str,
vector<string> *codes);


#endif /* MTB_DOUBLEMETAPHONE_H_ */
