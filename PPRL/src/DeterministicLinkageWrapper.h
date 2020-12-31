//============================================================================
// Name        : DeterministicLinkageWrapper.h
// Author      : Rukasz
// Version     :
// Copyright   : Your copyright notice
// Description : Ansi-style
//============================================================================




#ifndef DETERMINISTICLINKAGEWRAPPER_H
#define DETERMINISTICLINKAGEWRAPPER_H
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
 * Calculates distance withou em
 */
MTB_Result distanceWrapper(MTB_StringVectorData d1, MTB_StringVectorData d2, vector<MergingConfiguration> mc, int counterSim);
/**
 * Calculates qualities withou em, when more than one distance is calculated
 */
MTB_Result multiDistanceWrapper(MTB_StringVectorData d1, MTB_StringVectorData d2,
                           vector<MergingConfiguration> mc);


#endif //
