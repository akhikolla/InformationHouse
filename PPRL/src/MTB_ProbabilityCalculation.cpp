/*
 * MTB_ProbabilityCalculation.cpp
 *
 *  Created on: 16.05.2017
 *      Author: schnell-42
 */

#include "MTB_ProbabilityCalculation.h"
#include "MergingConfiguration.h"

MTB_ProbabilityCalculation::MTB_ProbabilityCalculation() {
}

MTB_ProbabilityCalculation::~MTB_ProbabilityCalculation() {
}

vector<int> MTB_ProbabilityCalculation::countFrequencies(
		vector<string> data1, vector<string> data2) {
  // get the mergeset
	vector<int> frequencyCount(round(pow(2, 1)));
	int totalComparisons =  data1.size()*data2.size();// each entry of DataA is compared to each entry in DataB

	this->value = 0;

	int doneComparisons = 0;
	int onePercent = totalComparisons / 100;
	if (onePercent == 0) {
		onePercent = 1;
	}


		for (unsigned j = 0; j < data1.size(); j++) {

//

			for (unsigned k = 0; k < data2.size(); k++) {

				int pattern = 0;
//								// only the first
//								// variable of each array is taken
					string val1 = data1[j];
					string val2 = data2[k];

					// check for equality
					if (!val1.empty() && !val2.empty() && (val1 == val2)) {
											pattern += round(pow(2, 0));
										}

				doneComparisons++;
				frequencyCount[pattern]++;

				if ((doneComparisons % 1000) == 0) {
					this->value = doneComparisons / onePercent;

				}
			}
		}

	return frequencyCount;
}

vector<MergingConfiguration> MTB_ProbabilityCalculation::run(vector<MergingConfiguration> mc,
		vector<string> data1, vector<string> id1, vector<string> data2, vector<string> id2, MTB_EMAlgorithm em) {
//for blocks, get initial m array, get initial u array
	// only calculate m, if "use global m" is set
	if (mc.size() > 0 && mc[0].getUseGlobalM()) {
		//calculates the Expectation-Maximisation Algorithm
		if (!em.calculate()) {
			return mc;
		}
		mc = setMArray(em.getMArray(), mc);
	}
//		// only calculate u, if "use global u" is set
	if (mc.size() > 0 && mc[0].getUseGlobalU()) {
		if (data1.size() <= 0 && data2.size() <= 0) {
			return mc;
		}
		vector<int> allFrequencies = countFrequencies(data1, data2);
		vector<double> uArray = getUArray(convertToSingleFrequencies(allFrequencies), sum(allFrequencies));
		for (double u : uArray) {

			if (u == 1) {
				Rcpp::Rcerr << "Some u value is 1" << endl;
				return mc;
			}
		}
		mc = setUArray(uArray, mc);
	}
	return mc;
}

vector<MergingConfiguration> MTB_ProbabilityCalculation::setMArray(vector<double> mArray,
		vector<MergingConfiguration> mergingConfigurationList) {
	if ((mergingConfigurationList.size() != 0) && (mArray.size() != 0)) {
		for (unsigned i = 0; (i < mergingConfigurationList.size()) && (i < mArray.size()); i++) {
			mergingConfigurationList[i].setM(mArray[i]);
		}
	}
	return mergingConfigurationList;
}

vector<MergingConfiguration> MTB_ProbabilityCalculation::setUArray(vector<double> uArray,
		vector<MergingConfiguration> mergingConfigurationList) {
	if ((mergingConfigurationList.size() > 0) && (uArray.size() > 0)) {
		for (unsigned i = 0; (i < mergingConfigurationList.size()) && (i < uArray.size()); i++) {
			mergingConfigurationList[i].setU(uArray[i]);
		}
	}
	return mergingConfigurationList;
}

vector<double> MTB_ProbabilityCalculation::getUArray(vector<int> frequencies, int compFreq) {
	vector<double> uArray(frequencies.size());
	if ((frequencies.size() == 0) || (compFreq == 0)) {
		return uArray;
	}

	for (unsigned i = 0; i < uArray.size(); i++) {
		uArray[i] = (double) frequencies[i] / (double) compFreq;
	}
	return uArray;
}
