/*
 * MTB_Reth.h
 *
 *  Created on: 05.07.2017
 *      Author: schnell-42
 */

#ifndef MTB_RETH_H_
#define MTB_RETH_H_


#include<vector>
#include <string>
#include "MTB_Similarity.h"
#include "Util.h"
#include "Standardisation.h"

class MTB_Reth: public MTB_Similarity {

public:
	MTB_Reth() {
	}
	;
	~MTB_Reth() {
	}
	;

	//
	string getName();
	double getRelativeValue(string o1, string o2);
	double getAbsoluteValue(string o1, string o2);

//	 double getThreshold();
//	 bool isAbsolute();
//	 bool isSimilar(string o1, string o2);
//	 void setAbsolute(bool absolute);
//	 void setThreshold(double threshold);

};
string encodeReth(string input);



#endif /* MTB_RETH_H_ */
