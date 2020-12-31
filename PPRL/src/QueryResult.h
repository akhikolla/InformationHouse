// QueryResult.h
//
// Copyright (c) 2016
// Prof. Dr. Rainer Schnell
// Universitaet Duisburg-Essen
// Campus Duisburg
// Institut fuer Soziologie
// Lotharstr. 65
// 47057 Duisburg
//
// This file is part of the command line application "mbtSearch".
//
// "mbtSearch" is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// "mbtSearch" is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with "mbtSearch". If not, see <http://www.gnu.org/licenses/>.

#ifndef QUERYRESULT_H
#define QUERYRESULT_H

#include <stdlib.h>
#include <stdio.h>
#include <pthread.h>
#include <string.h>
#include "CLK.h"

// Instances of QueryResultNode store one CLK and
// the corresponding similarity coefficient within a sorted
// tree. If sorting is disabled, QueryResultNodes are simply
// used as nodes of a linked list.

typedef struct QueryResultNodeStruct {
	char *mQueryId;				// ID of query CLK
  char *mMatchId;				// ID of match CLK
	CLK *mPrint;			// pointer to matching CLK
	float mSimilarity;			// corresponding similarity coefficient
	struct QueryResultNodeStruct *mLeft;	// left subtree (higher similarity coeffs)
	struct QueryResultNodeStruct *mRight;	// right subtree (lower similarity coeffs)
} QueryResultNode;

// Objects of class QueryResult store search results,
// which consist of a CLK and the computed
// similarity coefficient. The QueryResult is designed as
// a sorted tree. On construction the QueryResult can be
// flagged to sort the results while collecting them by
// the similarity coefficient or not.

class QueryResult {
	private:

	QueryResultNode *mRootNode;		// root node of sorted tree or linked list
	int64_t mSize;				// size of result set
	int mSort;				// flag if the query results have ro be sorted
	pthread_mutex_t mAddMutex;		// thread save insertion mutex
	pthread_mutex_t mCntXORMutex;		// thread save counter mutex
	pthread_mutex_t mCntSimMutex;		// thread save counter mutex
	const char *mSeperator;			// column seperator for csv output
	int64_t mCntXOR;			// statistic counter for XOR-Hash comparisons
	int64_t mCntSim;			// statistic counter for similarity comparisons
	int64_t mSizeTree;			// statistic tree size info
	int64_t mSizeSearch;			// statistic search size info

	void deleteNode(QueryResultNode *node);	// delete a node and all its subnodes

	public:

	// constructor
	QueryResult(int sort, int64_t sizeTree);

	// destructor
	~QueryResult();

  // add Id of the CLKs and  the corresponding similarity coefficient to the query result
  void add(char *queryId, char *matchId, double similarity);
	// add a CLK and the corresponding similarity coefficient to the query result
	void add(CLK *queryCLK, CLK *clk, float similarity);

	// return root node for reading
	inline QueryResultNode *getRootNode() {
		return mRootNode;
	}

	// return size of result set
	inline int64_t getSize() {
		return mSize;
	}

	// get query statistics
	void getStatistics(double *valuesPtr, double *percentsPtr);

	// increment XOR counter
	inline void countXOR() {
		if (mSizeTree > 0) {
			// lock mutex
			pthread_mutex_lock(&mCntXORMutex);

			// increment counter
			mCntXOR++;

			// unlock mutex
			pthread_mutex_unlock(&mCntXORMutex);
		}
	}

	// get XOR counter
	inline int64_t getCntXOR() {
		return mCntXOR;
	}

	// increment similarity counter
	inline void countSim() {
		if (mSizeTree > 0) {
			// lock mutex
			pthread_mutex_lock(&mCntSimMutex);

			// increment counter
			mCntSim++;

			// unlock mutex
			pthread_mutex_unlock(&mCntSimMutex);
		}
	}

	// get similarity counter
	inline int64_t getCntSim() {
		return mCntSim;
	}

        // set size of last search;
	inline void setSizeLastSearch(int64_t sls) {
		mSizeSearch = sls;
	}
};
#endif
