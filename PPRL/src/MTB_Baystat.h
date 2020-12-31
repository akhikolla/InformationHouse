/*
 * MTB_Baystat.h
 *
 *  Created on: 05.07.2017
 *      Author: schnell-42
 */

#ifndef MTB_BAYSTAT_H_
#define MTB_BAYSTAT_H_


#include <algorithm>
#include <string>
#include <ctype.h>
#include <iostream>
#include <cstring>
#include <vector>
#include "MTB_Similarity.h"
#include "Util.h"
#include "Standardisation.h"


/**
 * 		The "Phonetikmodul" of the Landesamt fuer Statistik und Datenverarbeitung as described in `Zusammenfuehrung von Datenbestaenden ohne
 *         numerische Identifikatoren', Bayern in Zahlen, Juli 2002, Heft 7 Release 1.0 match = baystat($string1, $string2 [, $minSubstringLength, $rangeLeft,
 *         $rangeRight] ) Input: string1, string2 Optional: minSubstringLength (if left out or set to -1, DEFAULT_MINSUBSTRINGLENGTH is chosen) Optional:
 *         rangeLeft, rangeRight (if left out or set to -1 DEFAULT_RANGELEFT/DEFAULT_RANGERIGHT are chosen) Returns match = 0..100 (0 means no common
 *         Substrings, 100 means identical...)
 *
 *         Baystat (Fürnrohr, Rimmelspacher, & Roncador, 2002)
 *
 */

class MTB_Baystat: public MTB_Similarity {

public:
	MTB_Baystat() {
	}
	;
	~MTB_Baystat() {
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
private:
	int	minSubstringLengthLongerEqual7	= 2;
	int	minSubstringLengthSmaller7		= 1;
	int	rangeLeftLongerEqual7			= 2;
	int	rangeRightLongerEqual7			= 2;
	int	rangeLeftSmaller7				= 1;
	int	rangeRightSmaller7				= 1;

};

int maxLength(string str1, string str2);

#endif /* MTB_BAYSTAT_H_ */
