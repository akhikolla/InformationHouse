//============================================================================
// Name        : createCLKWrapper.h
// Author      : Rukasz
// Version     :
// Copyright   : Your copyright notice
// Description : Ansi-style
//============================================================================
// Enable C++11 via this plugin
// [[Rcpp::plugins("cpp11")]]



#ifndef CREATECLKWRAPPER_H
#define CREATECLKWRAPPER_H

#include <Rcpp.h>
#include "CreateRecordLevelBF.h"
#include "CreateCLKBigramSeed.h"
#include <numeric>
using namespace std;
using namespace Rcpp;
double delPunktTest(std::string test);

DataFrame CreateCLKc(SEXP ID, DataFrame data, SEXP password, int k, IntegerVector padding, IntegerVector qgram, unsigned lenBloom );

DataFrame CreateBFc(SEXP ID, SEXP data, SEXP password, int k, int padding, int qgram, unsigned lenBloom);

void freeCLKArray(CLK** clk, int n);

#endif
