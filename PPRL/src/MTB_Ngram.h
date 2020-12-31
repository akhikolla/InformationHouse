/*
 * Ngram.h
 *
 *  Created on: Feb 20, 2017
 *      Author: rukasz
 *      see (Ukkonen, 1992)
 */

#ifndef MTB_NGRAM_H_
#define MTB_NGRAM_H_
#include <iostream>
#include <stdio.h>
#include <vector>
#include <algorithm>
#include <Rcpp.h>
#include "MTB_Similarity.h"
#include "BloomHelper.h"

using namespace std;
/**
* Creates q-Grams of size lenQgram from the input string.
* Only uni-, bi-, tri- and fourgram are allowed (for now)
*
* @param var is the input string
* @param lenQgram is the size of the qgram
* @return the input string transformed into a vector of qgrams
*/
vector<string> CreateQgrams(string var, unsigned lenQgram);

/**
 * Computes the distance of the q-Grams of two strings according to Ukkonen, 1992.
 * Only uni-, bi-, tri- and fourgram are allowed (for now).
 *
 *
 * @param var is the input string
 * @param lenQgram is the size of the qgram
 * @return the distance of the q-Grams of two strings
 */
double NgramDistance(vector<string> qgram1s, vector<string> qgrams2s);


class MTB_NgramDistanceAlgorithm: public MTB_Similarity {

public:
  MTB_NgramDistanceAlgorithm(int l1_, int l2_) {
    this->l1 = l1_;
    this->l2 = l2_;
  }
  ;
  ~MTB_NgramDistanceAlgorithm() {
  }
  ;

  //
  string getName();
  double getRelativeValue(string o1, string o2);
  double getAbsoluteValue(string o1, string o2);

  int getL1(){return l1;};
  int getL2(){return l2;};
  void setL1(int l1_){this->l1 = l1_;};
  void setL2(int l2_){this->l2 = l2_;};
  //	 double getThreshold();
  //	 bool isAbsolute();
  //	 bool isSimilar(string o1, string o2);
  //	 void setAbsolute(bool absolute);
  //	 void setThreshold(double threshold);
private:
  int l1 = 2, l2 = 2;

};

#endif /* MTB_NGRAM_H_ */
