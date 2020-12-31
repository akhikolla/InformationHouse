#ifndef CREATEESL_H
#define CREATEESL_H

#include <iostream>
#include <stdio.h>
#include <vector>
#include <cstring>
#include <memory>
#include <random>
#include <iterator>
#include <string.h>
#include "Standardisation.h"
#include "Hashing.h"

using namespace std;

/**
 * Creating esl from am input vector,
 * concatenating the results of esl (if desired) of each variable of the input vector.
 *
 * @param vars
 * @param code
 * @param password
 * @return
 */
string createESL(vector<string> vars, vector<vector<int>> code, string password);

/**
 * Applies esl if desired on the input string according to the given code.
 * If the code is {0} nothing will be done and the output is the same as the input, for letter in upper case.
 *
 * @param var
 * @param code
 * @return
 */
string createESLHelper(string var, vector<int> code);


#endif
