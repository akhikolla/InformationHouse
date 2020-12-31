/*
 * MTB_Similarity.h
 *
 *  Created on: 21.06.2017
 *      Author: rukasz
 */

#ifndef MTB_SIMILARITY_H_
#define MTB_SIMILARITY_H_

#include <stdio.h>
#include <string>
#include <vector>

using namespace std;

/**
 * similarities
 */
class MTB_Similarity {

public:
	MTB_Similarity() {
	}
	virtual ~MTB_Similarity() {
	}

	string getName();

	virtual double getRelativeValue(string o1, string o2)=0;
	virtual double getAbsoluteValue(string o1, string o2)=0;
//	virtual double getThreshold()=0;
//	virtual bool isAbsolute()=0;
//	virtual bool isSimilar(string o1, string o2)=0;
//	virtual void setAbsolute(bool absolute)=0;
//	virtual void setThreshold(double threshold)=0;


};

#endif /* MTB_SIMILARITY_H_ */
