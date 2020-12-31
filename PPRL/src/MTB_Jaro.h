/*
 * MTB_Jaro.h
 *
 *  Created on: 10.04.2017
 *      Author: schnell-42
 */

#ifndef MTB_JARO_H_
#define MTB_JARO_H_

#include <string>
#include <iostream>
#include <typeinfo>
#include <algorithm>
#include "MTB_Similarity.h"

using namespace std;

const double JARO_WEIGHT_STRING_A(1.0/3.0);
const double JARO_WEIGHT_STRING_B(1.0/3.0);
const double JARO_WEIGHT_TRANSPOSITIONS(1.0/3.0);

const unsigned long int JARO_WINKLER_PREFIX_SIZE(4);
const double JARO_WINKLER_SCALING_FACTOR(0.1);
const double JARO_WINKLER_BOOST_THRESHOLD(0.7);


double jaroWinklerDistance(const std::string&a, const std::string& b);
double jaroDistance(const std::string& a, const std::string& b);


class MTB_JaroAlgorithm: public MTB_Similarity {

public:
	MTB_JaroAlgorithm() {
	}
	;
	~MTB_JaroAlgorithm() {
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

class MTB_JaroWinklerAlgorithm: public MTB_Similarity {

public:
	MTB_JaroWinklerAlgorithm() {
	}
	;
	~MTB_JaroWinklerAlgorithm() {
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

#endif /* MTB_JARO_H_ */
