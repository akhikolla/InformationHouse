/*
 * Util.h
 *
 *  Created on: 05.07.2017
 *      Author: schnell-42
 *
 *      Help routines
 */

#ifndef UTIL_H_
#define UTIL_H_

#include <algorithm>
#include <string>
#include <ctype.h>
#include <iostream>
#include <cstring>
#include <vector>
#include "Standardisation.h"
#include <math.h>

using namespace std;


//void replaceAll(string& str, const string& from,const string& to);
//void replaceSharpS(string& str);
//void replaceUmlauts(string& str);
//void toUpperCase(string& str);
string replaceEnding(string fullString, string ending, char repl);
vector<int> convertToSingleFrequencies(vector<int> patternFrequencies);
int sum(vector<int> values);

#endif /* UTIL_H_ */
