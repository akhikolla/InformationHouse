/*
 * MTB_MergeData.h
 *
 *  Created on: 24.05.2017
 *      Author: schnell-42
 */

#ifndef MTB_MERGEDATA_H_
#define MTB_MERGEDATA_H_

#include <string>
#include <vector>
#include <map>
#include <Rcpp.h>
#include "MTB_Data.h"
#include "MergingConfiguration.h"
#include "MTB_Blocking.h"

using namespace std;

/**
 * Holds merge data
 */

class MTB_MergeData {
public:
	/**
	 * constructor
	 */
	MTB_MergeData();
	~MTB_MergeData();

	/**
	 * Returns data1
	 *
	 * @param
	 * @param
	 */
	MTB_StringVectorData getData1();

	/**
		 * Returns data2
		 *
		 * @param
		 * @param
		 */
	MTB_StringVectorData getData2();

	/**
	 * Sets data
	 *
	 * @param
	 * @param
	 */
	void setData(MTB_StringVectorData data1, MTB_StringVectorData data2);






private:
	MTB_StringVectorData data1;
	MTB_StringVectorData data2;
};

/**
 * Initializes Merge Data
 * Applies blocking
 *
 * @param d1
 * @param d2
 * @param mc
 * @return
 */
vector<MTB_MergeData> initMergeData(MTB_StringVectorData d1, MTB_StringVectorData d2,
		vector<MergingConfiguration> mc);

#endif /* MTB_MERGEDATA_H_ */
