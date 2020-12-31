/*
 * InternalProbabilisticComparisonsObservation.cpp
 *
 *  Created on: 20.06.2017
 *      Author: schnell-42
 */

/**
* Calculating the merge
*/
#include <fstream>      // std::ofstream

#include "MTB_InternalProbabilisticComparisonsObservation.h"

MTB_InternalProbabilisticComparisonsObservation::MTB_InternalProbabilisticComparisonsObservation() {
}

MTB_InternalProbabilisticComparisonsObservation::~MTB_InternalProbabilisticComparisonsObservation() {
}

double MTB_InternalProbabilisticComparisonsObservation::calculateQuality(
    MergingConfiguration mergingConfiguration, string val1, string ID1, string val2, string ID2,
    double jaroWeightFactor, int position, MTB_Result *res1) {

  MTB_Similarity *algorithm = mergingConfiguration.getAlgorithm();
  if (algorithm == NULL) {
    Rcpp::Rcerr << "no algorithm set!" << algorithm << endl;
    delete algorithm;
    return 0.0;
  }

  double m = mergingConfiguration.getM();
  double u = mergingConfiguration.getU();
  //double maxM = mergingConfiguration.getMaxM();
  //double minU = mergingConfiguration.getMinU();
  double missingValue = calculateMissingValue(mergingConfiguration);

  if ((val1.length() < 1) || (val2.length() < 1)) {
    // missing value
    res1->addResult(ID1, ID2, missingValue);
    delete algorithm;
    return missingValue;
  }

  // calculate similarity
  double agreement = algorithm->getRelativeValue(val1, val2);
  double wa = log2(m / u);
  double wd = log2((1 - m) / (1 - u));

  // fellegi-sunter, if similarity is 0 or 1
  if (agreement == 1) {
    res1->addResult(ID1, ID2, wa);
    delete algorithm;
    return wa;
  } else if (agreement == 0) {
    res1->addResult(ID1, ID2, wd);
    delete algorithm;
    return wd;
  }

  double res = max(wa - ((wa - wd) * (1 - agreement) * jaroWeightFactor), wd);
  // if (res > 0.0) {
  // } else if (res <= 0.0) {
  // } else {
  // }
  res1->addResult(ID1, ID2, res);
  delete algorithm;
  return res;
}

double MTB_InternalProbabilisticComparisonsObservation::calculateMissingValue(
    MergingConfiguration mergingConfiguration) {
  if (mergingConfiguration.getMissingWeight() != MISSING_VALUE_AVERAGE) {
    return mergingConfiguration.getMissingWeight();
  }
  double m = mergingConfiguration.getM();
  double u = mergingConfiguration.getU();
  double wa = log2(m / u);
  double wd = log2((1 - m) / (1 - u));
  return (wa + wd) / 2;
}

double MTB_InternalProbabilisticComparisonsObservation::calculateDistance(
    MergingConfiguration mergingConfiguration, string val1, string val2) {


  MTB_Similarity *algorithm = mergingConfiguration.getAlgorithm();
  if (algorithm == NULL) {
    Rcpp::Rcerr << "no algorithm set!" << algorithm << endl;
    delete algorithm;
    return 0.0;
  }

  // calculate similarity
  double distance = algorithm->getAbsoluteValue(val1, val2);
  delete algorithm;
  return distance;
}

