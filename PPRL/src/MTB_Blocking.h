/*
 * MTB_Blocking.h
 *
 *  Created on: 28.07.2017
 *      Author: rukasz
 */

#ifndef MTB_BLOCKING_H_
#define MTB_BLOCKING_H_

#include <iostream>
#include <vector>
#include <string>
#include <numeric>      // std::iota
#include <algorithm>    // std::sort
#include "MTB_Data.h"
#include "MTB_Exact.h"
#include <time.h>

using namespace std;
/**
 * Blocking
 */

vector<MTB_StringVectorData> exactBlocking(MTB_StringVectorData);
vector<MTB_StringVectorData> exactCLBlocking(MTB_StringVectorData);

/**
 * sorts bocking vectors, remembers the original indexes
 */
vector<size_t> sort_blockingData(const vector<string> &v);

#endif /* MTB_BLOCKING_H_ */
