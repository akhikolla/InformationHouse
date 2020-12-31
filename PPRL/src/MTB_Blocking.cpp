/*
 * MTB_Blocking.cpp
 *
 *  Created on: 28.07.2017
 *      Author: rukasz
 */

#include "MTB_Blocking.h"

/**
 * sorts blocking data
 * removes duplicates
 */
vector<MTB_StringVectorData> exactBlocking(MTB_StringVectorData data) {
	vector<MTB_StringVectorData> res;
	vector<size_t> idx = sort_blockingData(data.getBlockingData());
	MTB_StringVectorData tmp;
	vector<string> tmpData;
	vector<string> tmpID;
	vector<string> tmpBlockingData;
	int counter = 0;

	for (unsigned i = 0; i < data.getBlockingData().size() - 1; i++) {

		if (data.getBlockingData()[idx[i]] == data.getBlockingData()[idx[i + 1]]) {
			tmpData.push_back(data.getData()[idx[i]]);
			tmpID.push_back(data.getID()[idx[i]]);
			tmpBlockingData.push_back(data.getBlockingData()[idx[i]]);
			counter++;
			if (i == data.getBlockingData().size() - 2) { //last value
				tmpData.push_back(data.getData()[idx[i + 1]]);
				tmpID.push_back(data.getID()[idx[i + 1]]);
				tmpBlockingData.push_back(data.getBlockingData()[idx[i + 1]]);
			}
		} else if (data.getBlockingData()[idx[i]] != data.getBlockingData()[idx[i + 1]]) {
			counter = 0;
			tmpData.push_back(data.getData()[idx[i]]);
			tmpID.push_back(data.getID()[idx[i]]);
			tmpBlockingData.push_back(data.getBlockingData()[idx[i]]);
			tmp.setData(tmpData);
			tmp.setID(tmpID);
			tmp.setBlockingData(tmpBlockingData);
			res.push_back(tmp);
			tmpData.clear();
			tmpID.clear();
			tmpBlockingData.clear();
			if (i == data.getBlockingData().size() - 2) { //last value
				tmpData.push_back(data.getData()[idx[i + 1]]);
				tmpID.push_back(data.getID()[idx[i + 1]]);
				tmpBlockingData.push_back(data.getBlockingData()[idx[i + 1]]);
				tmp.setData(tmpData);
				tmp.setID(tmpID);
				tmp.setBlockingData(tmpBlockingData);
				res.push_back(tmp);
			}
		}

	}
	//last values
	if (counter > 0) {
		tmp.setData(tmpData);
		tmp.setID(tmpID);
		tmp.setBlockingData(tmpBlockingData);
		res.push_back(tmp);
		tmpData.clear();
		tmpID.clear();
		tmpBlockingData.clear();
	}
	return res;
}


/**
 *
 */
vector<MTB_StringVectorData> exactCLBlocking(MTB_StringVectorData data) {
	vector<MTB_StringVectorData> res;
	vector<size_t> idx = sort_blockingData(data.getBlockingData());
		MTB_StringVectorData tmp;
		vector<string> tmpData;
		vector<string> tmpID;
		vector<string> tmpBlockingData;
		int counter = 0;
		for (unsigned i = 0; i < data.getBlockingData().size() - 1; i++) {

			if ( StringToUpper(data.getBlockingData()[idx[i]]).compare(StringToUpper( data.getBlockingData()[idx[i + 1]]))==0) {
				tmpData.push_back(data.getData()[idx[i]]);
				tmpID.push_back(data.getID()[idx[i]]);
				tmpBlockingData.push_back(data.getBlockingData()[idx[i]]);
				counter++;
				if (i == data.getBlockingData().size() - 2) { //last value
					tmpData.push_back(data.getData()[idx[i + 1]]);
					tmpID.push_back(data.getID()[idx[i + 1]]);
					tmpBlockingData.push_back(data.getBlockingData()[idx[i + 1]]);
				}
			} else if (StringToUpper(data.getBlockingData()[idx[i]]).compare(StringToUpper( data.getBlockingData()[idx[i + 1]]))!=0) {
				counter = 0;
				tmpData.push_back(data.getData()[idx[i]]);
				tmpID.push_back(data.getID()[idx[i]]);
				tmpBlockingData.push_back(data.getBlockingData()[idx[i]]);
				tmp.setData(tmpData);
				tmp.setID(tmpID);
				tmp.setBlockingData(tmpBlockingData);
				res.push_back(tmp);
				tmpData.clear();
				tmpID.clear();
				tmpBlockingData.clear();
				if (i == data.getBlockingData().size() - 2) { //last value
					tmpData.push_back(data.getData()[idx[i + 1]]);
					tmpID.push_back(data.getID()[idx[i + 1]]);
					tmpBlockingData.push_back(data.getBlockingData()[idx[i + 1]]);
					tmp.setData(tmpData);
					tmp.setID(tmpID);
					tmp.setBlockingData(tmpBlockingData);
					res.push_back(tmp);
				}
			}

		}
		//last values
		if (counter > 0) {
			tmp.setData(tmpData);
			tmp.setID(tmpID);
			tmp.setBlockingData(tmpBlockingData);
			res.push_back(tmp);
			tmpData.clear();
			tmpID.clear();
			tmpBlockingData.clear();
		}
		return res;
}

/**
 * sorts bocking vectors, remembers the original indexes
 */
vector<size_t> sort_blockingData(const vector<string> &v) {

	// initialize original index locations
	vector<size_t> idx(v.size());
	iota(idx.begin(), idx.end(), 0);

	// sort indexes based on comparing values in v
	sort(idx.begin(), idx.end(),
			[&v](size_t i1, size_t i2) {return v[i1] < v[i2];});

	return idx;
}
