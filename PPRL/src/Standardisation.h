//============================================================================
// Name        : Standardisation.h
// Author      : Rukasz
// Version     :
// Copyright   : Your copyright notice
// Description : Ansi-style
//============================================================================


#ifndef STANDARDISATION_H
#define STANDARDISATION_H

#include <iostream>
#include <vector>
#include <algorithm>
#include <string>

using namespace std;

//void replaceAll(string &s, char x, char y);
//static inline void ReplaceAllSubstr(string str, const string& from, const string& to);
void ReplaceAllSubstr( string& source, const string& from, const string& to );
//void removeAll(string &s, char x);
void toUpper(string &s);
void ReplaceSpecials(string &CaseString);
bool invalidChar(char c) ;
bool invalidCharPunct(char c);

/**
 * Makes ALL CAPS.
 *
 * @param var
 * @return standardized var
 */
 void preprocess(string& var);

/**
 * Makes ALL CAPS,
 * replaces Titel, Fuellsel and all non
 * resolve umlauts, foreign and Special Chars,
 * only ASCII is kept.
 *
 * @param var
 * @return standardized var
 */
void standardisation(string& var);

/**
 * Makes ALL CAPS,
 * resolve umlauts, foreign and Special Chars,
 * only ASCII is kept.
 *
 * @param var
 * @return standardized var
 */
void replaceNonAscii(string& var);

/**
 * Deletes titles.
 *
 * @param var
 * @return  var without titles
 */
void delTitel(string& var) ;

/**
 * Deletes all unwanted elements of names ('fuellsel').
 *
 * @param var
 * @return var without fuellsel
 */
void delFuell(string& var);

/**
 * Deletes punct, blanc and control signs.
 * @param var
 * @return
 */
void deletePunct(string& var);

#endif
