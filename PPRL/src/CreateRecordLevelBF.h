//============================================================================
// Name        : CreateRecordLevelBF.h
// Author      : Rukasz
// Version     :
// Copyright   : Your copyright notice
// Description : Ansi-style
//============================================================================


#ifndef CREATERECORDLEVELBF_H
#define CREATERECORDLEVELBF_H

#include "BloomHelper.h"
#include "CreateCLKBigramSeed.h"
#include <math.h>
#include <algorithm>    // std::min_element, std::max_element
#include "CLK.h"

using namespace std;


void CreateRBF(CLK* rlbf,vector<string> vars, int k, vector<int> padding, vector<int> lenQgram, unsigned lenBloom,int lenRLBF, vector<string> password);


#endif
