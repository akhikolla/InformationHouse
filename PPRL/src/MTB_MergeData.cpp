/*
 * MTB_MergeData.cpp
 *
 *  Created on: 29.05.2017
 *      Author: rukasz
 */

#include <iostream>
#include "MTB_MergeData.h"


MTB_MergeData::MTB_MergeData(){}
MTB_MergeData::~MTB_MergeData(){}

MTB_StringVectorData MTB_MergeData::getData1(){return data1;}
MTB_StringVectorData MTB_MergeData::getData2(){return data2;}

void MTB_MergeData::setData(MTB_StringVectorData data1, MTB_StringVectorData data2){
	this->data1 = data1;
	this->data2 = data2;
}

vector<MTB_MergeData> initMergeData(MTB_StringVectorData d1, MTB_StringVectorData d2,
		vector<MergingConfiguration> mc){
	//set MergeData

		vector<MTB_MergeData> mergeData;
    mergeData.clear();
		MTB_MergeData tmp;
		 if (mc[0].getBlocking()!= "0"&&mc[0].getBlocking()!= "mbt") {
		 	vector<MTB_StringVectorData> blocked1;
		 	vector<MTB_StringVectorData> blocked2;
		 	if (mc[0].getBlocking()== "exact") {
		 	  blocked1 = exactBlocking(d1);
		 		blocked2 = exactBlocking(d2);
			}
			if (mc[0].getBlocking()== "exactCL") {
						blocked1 = exactCLBlocking(d1);
						blocked2 = exactCLBlocking(d2);
					}
			for (unsigned i = 0; i < blocked1.size(); i++) {
				for (unsigned j = 0; j < blocked2.size(); j++) {
					if (blocked1[i].getData().size() > 0
							&& blocked2[j].getData().size() > 0) {

						if (blocked1[i].getBlockingData()[0]
								== blocked2[j].getBlockingData()[0]) {
							tmp.setData(blocked1[i], blocked2[j]);
							mergeData.push_back(tmp);
						}
					}
				}
		}
		 } else { // without blocking
		 	tmp.setData(d1, d2);
		 	mergeData.push_back(tmp);
		 }
		 if (mergeData.size() == 0){
		   Rcpp::Rcerr << "Nothing to link since blocking filtered everything!" << endl;
		 }
		return mergeData;
}
