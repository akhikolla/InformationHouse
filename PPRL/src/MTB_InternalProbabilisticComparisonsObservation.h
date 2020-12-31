/*
 * MTB_InternalProbabilisticComparisonsObservation.h
 *
 *  Created on: 20.06.2017
 *      Author: schnell-42
 */

#ifndef MTB_INTERNALPROBABILISTICCOMPARISONSOBSERVATION_H_
#define MTB_INTERNALPROBABILISTICCOMPARISONSOBSERVATION_H_

#include <iostream>
#include <stdio.h>
#include <vector>
#include "MergingConfiguration.h"
#include "MTB_MergeData.h"
#include <math.h>       /* log2 */
#include "MTB_Similarity.h"
#include <Rcpp.h>
#include "MTB_Result.h"
using namespace std;


/**
 * Calculating the merge
 */
class MTB_InternalProbabilisticComparisonsObservation {
public:
	/** Constructor
	 *
	 */
	MTB_InternalProbabilisticComparisonsObservation();
	/** Destuctor
	 *
	 */
	~MTB_InternalProbabilisticComparisonsObservation();

	string CLASS_MATCH = "M";
	string CLASS_POSSIBLE_MATCH = "PM";
	string CLASS_NON_MATCH = "NM";


	/**
	 * Calculates the quality value
	 * @param mergingConfiguration
	 * @param var1
	 * @param var2
	 * @param jaroWeightFactor
	 * @param position
	 * @return
	 */
	double calculateQuality( MergingConfiguration mergingConfiguration,
			string var1, string ID1, string var2, string ID2, double jaroWeightFactor, int position, MTB_Result *res);
	/**
		 * Calculates the absolute distance according to the algorithm
		 * @param mergingConfiguration
		 * @param var1
		 * @param var2
		 * @return
		 */
	double calculateDistance( MergingConfiguration mergingConfiguration,
				string var1, string var2);

private:
	double MISSING_VALUE_AVERAGE = -1;
	vector<double> values;
	vector<bool> maximumAgreement;
	double quality = 0.0;
	string classification;


	double calculateArrayQuality(MergingConfiguration mergingConfiguration,
			int position);

	double calculateMissingValue(MergingConfiguration mergingConfiguration);

	double matchArrayVariables(MergingConfiguration mergingConfiguration,

			double jaroWeightFactor, bool checkMissing);
};

#endif /* MTB_INTERNALPROBABILISTICCOMPARISONSOBSERVATION_H_ */
