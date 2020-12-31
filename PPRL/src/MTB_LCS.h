/*
 * LCS.h
 *
 *  Created on: Feb 17, 2017
 *      Author: rukasz
 */

#ifndef MTB_LCS_H_
#define MTB_LCS_H_

/** Dynamic Programming C/C++ implementation of LCS problem
 * (Hirschberg, 1977)
 * */

#include<vector>
#include <string>
#include "MTB_Similarity.h"


using namespace std;

int max(int a, int b);
int lcsHelper(string X, string Y);
/** normalize result
 *
 * @param str1
 * @param str2
 * @param result LCS
 * @return
 */
double calculateCoeff(string str1, string str2, double result);

class MTB_LongestCommonSubsequenceAlgorithm: public MTB_Similarity {

public:
	MTB_LongestCommonSubsequenceAlgorithm() {
	}
	;
	~MTB_LongestCommonSubsequenceAlgorithm() {
	}
	;

	//
	string getName();
	double getRelativeValue(string o1, string o2);
	double getAbsoluteValue(string o1, string o2);
};

#endif /* MTB_LCS_H_ */
