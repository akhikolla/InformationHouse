/*
 * MergingConfiguration.cpp
 *
 *  Created on: 30.05.2017
 *      Author: schnell-42
 */

#include "MergingConfiguration.h"

MergingConfiguration::MergingConfiguration(double m, double u,
                                                   string algorithm) {
  this->m = m;
  this->u = u;
  this->algorithm = algorithm;
}

MergingConfiguration::MergingConfiguration() {
}

MergingConfiguration::~MergingConfiguration() {
}

void MergingConfiguration::setAlgorithm(string algorithm_){
  this->algorithm = algorithm_;
}

MTB_Similarity* MergingConfiguration::getAlgorithm() {

  if (this->algorithm == "lv") {
    this->similarity = new MTB_LevenshteinAlgorithm();
  }
  if (this->algorithm == "dl") {
    this->similarity = new MTB_DamerauLevenshteinAlgorithm();
  }
  if (this->algorithm == "exact") {
    this->similarity = new MTB_ExactAlgorithm();
  }
  if (this->algorithm == "exactCL") {
    this->similarity = new MTB_ExactCapitalLettersAlgorithm();
  }
  if (this->algorithm == "jaro") {
    this->similarity = new MTB_JaroAlgorithm();
  }
  if (this->algorithm == "jw") {
    this->similarity = new MTB_JaroWinklerAlgorithm();
  }
  if (this->algorithm == "ngram") {
    this->similarity = new MTB_NgramDistanceAlgorithm(this->lenNgram, this->lenNgram);
  }
  if (this->algorithm == "jw2") {
    this->similarity = new MTB_JaroWinklerMcLaughlinWinklerLynchAlgorithm();
  }
  if (this->algorithm == "LCS") {
    this->similarity = new MTB_LongestCommonSubsequenceAlgorithm();
  }
  if (this->algorithm == "Gcp") {
    this->similarity = new MTB_Baystat();
  }
  if (this->algorithm == "Reth") {
    this->similarity = new MTB_Reth();
  }
  if (this->algorithm == "Soundex") {
    this->similarity = new MTB_Soundex();
  }
  if (this->algorithm == "DoubleMetaphone") {
    this->similarity = new MTB_DoubleMetaphone();
  }
  if (this->algorithm == "Metaphone") {
    this->similarity = new MTB_Metaphone();
  }
  if (this->algorithm == "tanimoto") {
    this->similarity = new MTB_Tanimoto();
  }
  if (this->similarity == NULL) {
    Rcpp::Rcerr << " no valid algorithm set!" << endl;
    return NULL;
  }
  return this->similarity;
}

