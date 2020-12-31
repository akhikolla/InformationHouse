/*
 * soundex.cpp
 *
 *  Created on: 25.01.2017
 *      Author: schnell-42
 */
#include "MTB_Soundex.h"

double MTB_Soundex::getRelativeValue(string o1, string o2) {
	soundexC(o1);
	soundexC(o2);

	if (o1.compare(o2) == 0) {
		return 1.0;
	}
	return 0.0;

}

double MTB_Soundex::getAbsoluteValue(string o1, string o2) {
	return getRelativeValue(o1, o2);

}

string MTB_Soundex::getName() {
	return "Soundex";
}


