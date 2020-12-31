//============================================================================
// Name        : prepareData.h
// Author      : Rukasz
// Version     :
// Copyright   : Your copyright notice
// Description : Ansi-style
//============================================================================




#ifndef PREPAREDATA_H
#define PREPAREDATA_H

#include <algorithm> // for transform
#include <functional> // for plus
#include <Rcpp.h>

using namespace Rcpp;
using namespace std;


/**
 * gets data from different types
 */
vector<string> prepareData(SEXP data_, string dataName, bool em) ;

bool isAscii(string str, string dataName);

#endif
