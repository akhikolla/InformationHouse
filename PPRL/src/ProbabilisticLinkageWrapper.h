//============================================================================
// Name        : ProbabilisticLinkageWrapper.h
// Author      : Rukasz
// Version     :
// Copyright   : Your copyright notice
// Description : Ansi-style
//============================================================================




#ifndef PROBABILISTICLINKAGEWRAPPER_H
#define PROBABILISTICLINKAGEWRAPPER_H
#include <Rcpp.h>
#include <string>
#include "MTB_LCS.h"
#include "MTB_Jaro.h"
#include "MTB_JaroWinklerMcLaughlinWinklerLynch.h"
#include "MTB_Levenshtein.h"
#include "MTB_DamerauLevenshtein.h"
#include "MTB_Exact.h"
#include "MTB_EMAlgorithm.h"
#include "MTB_ProbabilityCalculation.h"
#include "MTB_InternalProbabilisticComparisonsObservation.h"
#include "MTB_Blocking.h"
#include "MTB_Result.h"
#include "MTB_MergeData.h"
#include "prepareData.h"
#include <algorithm> // for transform
#include <functional> // for plus

using namespace Rcpp;
using namespace std;


/**
 * Calculates qualities with em
 */MTB_Result emWrapper(MTB_StringVectorData d1, MTB_StringVectorData d2, vector<MergingConfiguration> mc, int counterSim);


#endif //PROBABILISTICLINKAGEWRAPPER_H
