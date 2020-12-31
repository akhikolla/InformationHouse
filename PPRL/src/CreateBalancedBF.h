/*
 * CreateBalancedBF.h
 *
 *  Created on: 30.08.2016
 *      Author: schnell-42
 */

#ifndef CREATEBALANCEDBF_H_
#define CREATEBALANCEDBF_H_

#include <iostream>
#include <algorithm>
#include <bitset>
#include <random>
#include <string>
#include <functional>
//#include <regex>
#include "ElegantPairingFunctions.h"
#include "CreateCLKBigramSeed.h"
using namespace std;

string CreateBalancedBloomfilterHelper(string clk, string password);
bool dateSalt(unsigned saltSize, unsigned pairingSalt1Size, unsigned pairingSalt2Size, unsigned dateSaltSize, unsigned varsSize);
bool saltVectors(unsigned saltSize, unsigned pairingSalt1Size, unsigned pairingSalt2Size, unsigned dateSaltSize, unsigned varsSize);
void checkVectors(unsigned varsSize,  vector<int>* padding,  vector<int> *qgram);
#endif /* CREATEBALANCEDBF_H_ */
