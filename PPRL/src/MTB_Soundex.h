/*
 * MTB_Soundex.h
 *
 *  Created on: 06.07.2017
 *      Author: schnell-42
 */

#ifndef MTB_SOUNDEX_H_
#define MTB_SOUNDEX_H_


#include <vector>
#include <string>
#include "MTB_Similarity.h"
#include "CreateALC.h"


class MTB_Soundex: public MTB_Similarity {

public:
	MTB_Soundex() {
	}
	;
	~MTB_Soundex() {
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



#endif /* MTB_SOUNDEX_H_ */
