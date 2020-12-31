/*
 * Levenshtein.h
 *
 *  Created on: 30.01.2017
 *      Author: schnell-42
 */

#ifndef MTB_LEVENSHTEIN_H_
#define MTB_LEVENSHTEIN_H_

#include <string.h>
#include <vector>
#include <numeric>
#include <iostream>
#include <algorithm>
#include "MTB_Similarity.h"

using namespace std;

class MTB_LevenshteinAlgorithm: public MTB_Similarity {

public:
	MTB_LevenshteinAlgorithm() {
	}
	;
	~MTB_LevenshteinAlgorithm() {
	}
	;

	 string getName();
	 double getRelativeValue(string o1, string o2);
	 double getAbsoluteValue(string o1, string o2);
//in case of non ascii characters errors may occure, so a flag is set for the user

//	 double getThreshold();
//	 bool isAbsolute();
//	 bool isSimilar(string o1, string o2);
//	 void setAbsolute(bool absolute);
//	 void setThreshold(double threshold);
//	 void Free(){ delete this; }

};

int levenshteinDistance(const std::string &s1, const std::string &s2, bool &nonAscii);

/**
 * Returns average length of the Strings in the given String Array
 *
 * @param strings
 *            Vector of strings
 * @return Average length
 */
double averageLengthVec(vector<string> strings);

/**
* Returns average length of 2 given Strings
*
* @param str1
*            String 1
* @param str2
*            String 2
* @return Average length
*/
double averageLength(string str1, string str2);
/**
	 * used for Levenshtein and Damerau-Levensthein
	 *
	 *
	 * @param str1
	 * @param str2
	 * @param result
	 * @return
	 */
double calcFromDistanceToInterval(string str1, string str2, double result);

#endif /* MTB_LEVENSHTEIN_H_ */
