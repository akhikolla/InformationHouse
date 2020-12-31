//============================================================================
// Name        : BloomHelper.h
// Author      : Rukasz
// Version     :
// Copyright   : Your copyright notice
// Description : Ansi-style
//============================================================================

#ifndef BLOOMHELPER_H
#define BLOOMHELPER_H

#include <iostream>
#include <stdio.h>
#include <cstring>
#include <memory>
#include <vector>
#include <random>
#include <algorithm>
#include <iterator>

#include "MTB_Ngram.h"
#include "Standardisation.h"

using namespace std;


//Definition of the used alphabet
#define ALPHABET " ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789"
/**
*
* Creats a vector of bloomfilters for according to the vector of input-strings.
*
* @param vars is a vector of variables/attributes
* @param lenBloom is the length of the bloomfilters
* @param k is the number of 1's for each bloomfilter
* @param lenQgram is vector of integers, indicating how long the qgrams should be for each attribute
* @param padding is a vector of integers, indicating how many spaces should be padded for each attribute
* @param lookupTable is vector of lookup tables, one for attribute
* @param uniBiGrams is a vector of uni- and bigrams
* @return a vector of bloomfilters, one for each attribute
*/
vector<string> CreateBloom(vector<string> vars, unsigned lenBloom, int k,
                         vector<int> lenQgram, vector<int> padding, vector<vector<int> > lookupTable,
                         vector<string> uniBiGrams);

/**
* Creates a vector of all possible uni- and bigrams.
* Basis is the above defined ALPHABET.
*
* @param uniBiGrams
* @return a vector of all possible uni- and bigrams
*/
vector<string> CreateUniBiGrams(vector<string> uniBiGrams);

/**
* Creates a vector of lookup tables for random hashing.
* The lookup table is a matrix of size lenDictionary times (k*number of variables) is created,
* where lenDictionary is the number of all possible unigrams and bigrams
* and k is the number of entries in the bloomfilter.
* For every entry in the Dictionary k random values are created.
* The matrix contains numbers from 0 to lenBloom (length of the Bloomfilter).
* A key/password/seed is necessary to build the matrix. For the same key the same matrix is build.
* The first entry for each variable is at k*position of the variable in the vector of variable.
*
* @param password is used to set a seed for each attribute, from those seeds further seeds
* 		are generated for each lookup table
* @param k is the number of 1's for each bloomfilter
* @param lenBloom is the length of the bloomfilters
* @param uniBiGrams
* @return a vector of lookup tables
*/
vector<vector<int> >  LookupTableBloom(int varsSize, string password, int k, unsigned lenBloom,
                                       vector<string> uniBiGrams);

/**
* Padding puts n spaces before and after the input string.
*
* @param var is the variable/attibute to be padded
* @param n number of spaces to be add
* @return padded string var
*/
void Padding(string &var, int n);



/**
 * Creates a single bloomfilter from the q-grams of size lenBloom
 *
 * @param qgrams is a vector of qgrams of the attribute
 * @param lenBloom is the length of the bloomfilters
 * @param k is the number of 1's for each bloomfilter
 * @param lookupTable, see LookupTableBloom
 * @param uniBiGrams is a vector of all possible uni- and bigrams
 * @return a vector of qgrams of the attribute
 */
string CreateBloomfilter(vector<string> qgrams, unsigned lenBloom, int k,
                       vector<vector<int> > lookupTable, vector<string> uniBiGrams, int i);

/**
 * Creates an empty bloomfilter of length lenBloom
 *
 * @param lenBloom
 * @return empty bloomfilter
 */
string CreateEmptyBloomfilter(int lenBloom);


#endif
