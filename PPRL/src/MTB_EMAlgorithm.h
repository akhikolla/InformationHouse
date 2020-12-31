/*
 * MTB_EMAlgorithm.h
 * see de.joergreiher.mergeToolBox.core;
 *  Created on: 16.05.2017
 *      Author: schnell-42
 */

#ifndef MTB_EMALGORITHM_H_
#define MTB_EMALGORITHM_H_

#include <iostream>
#include <vector>
#include <cmath>       /* pow */
#include <bitset> // for binary representation of an int
#include <stdlib.h>     /* abs */
#include <string>
#include <algorithm>
#include <climits>

using namespace std;

class MTB_EMAlgorithm {


  /**
  * Constructor for em algorithm
  *
  * @param Frequencies
  * @param m
  * @param u
  * @param p
  * @param epsilon
  */
public:
  MTB_EMAlgorithm(vector<int> frequencies, vector<double> m, vector<double> u, double p, double epsilon);

  //destuctor
  ~MTB_EMAlgorithm();
  /**
  *
  * @param vn
  * @return
  */
  vector<vector<int>> getPatternMatrix(int vn);
  bool calculate();
  string getException();
  int getIterations();
  vector<double> getMArray();
  bool checkValidity();
  double delta();
  // --------- E-Step --------- //
  void expectation(int i);
  // --------- M-Step --------- //
  void maximization(int i);

private:
  int	MAX_ITERATIONS	= 500;
  //float tol = 0.0;
  //double maxM = 0.9999;
  //double minU = 1.0E-4;

  string	EX_CONVERGENCY	= "EM doesn't converge";
  string	EX_M_NAN		= "Some m value is not a number";
  string	EX_M_0			= "Some m value is 0";
  // if difference of m and u of previous iteration smaller epsilon, stop
  // algorithm
  double	epsilon = 1.0e-4;
  // number of patterns
  int		patternCount;
  // number of record pairs
  int		recordPairs;
  // number of variables
  int		variableCount;
  vector<vector<double>>	m;
  vector<vector<double>>	u;
  vector<double>	p;
  // pattern frequencies
  vector<int>	patternFrequencies;
  // helper g_m, g_u
  vector<vector<double>>	gm;
  vector<vector<double>>	gu;
  vector<vector<int>>		patternMatrix;
  // number of iterationCount
  int		iterationCount = 0;
  string	exception = "";

  //from MTBSetThread
  //int	value = 0;
  //public:

};

#endif /* MTB_EMALGORITHM_H_ */
