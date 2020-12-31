#ifndef CREATEALC_H
#define CREATEALC_H

#include <iostream>
#include <stdio.h>
#include <vector>
#include <cstring>
#include <memory>
#include <random>
#include <iterator>
#include <Rcpp.h>
#include "Standardisation.h"
#include "Hashing.h"
using namespace std;

/**
 * Creating soundex/alc from am input vector,
 * concatenating the results of soundex (if desired) of each variable of the input vector.
 *
 * @param vars input vector with attributes
 * @param soundex boolean vector indicating for each attribute if soundex should be used
 * @param password
 * @return alc
 */
string createALC(vector<string> vars, vector<bool> soundex, string password);

/**
 *  Applies soundex if desired on the input string.
 *
 * @param var
 * @param soundex
 * @return
 */
void createALCHelper(string &var, bool soundex);

/**
 * Deletes vowels.
 *
 * @param var
 * @return
 */
void deleteVowels(string& var);

/**
 * Deletes Y, W and H if not at the beginning of the string.
 *
 * @param var
 * @return
 *
 */
void deleteYWH(string &var);

/**
 * Codes consonants except the first letter.
 * Deletes duplicates.
 *
 * @param var
 * @return
 */
void codeConsonants(string &var);

/**
 * Replaces dublicates in coded numbers
 * @param var
 */
void replaceDuplicates(string &var);

/**
 * Fills with zeros if code is shorter than three.
 *
 * @param var
 * @return
 */
void fillZero(string& var);

/**
 *  Cuts to code to three if it is longer.
 * @param var
 *
 */
void cutToThree(string &var);

/**
 *  Soundex is used on var.
 * @param var
 *
 */
void soundexC(string &var);
#endif

