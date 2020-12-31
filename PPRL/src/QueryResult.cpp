// QueryResult.cpp
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

#include "QueryResult.h"

// constructor
QueryResult::QueryResult(int sort,int64_t sizeTree) {
	mSize = 0;
	mRootNode = NULL;
	mSort = sort;

	mCntXOR = 0;
	mCntSim = 0;
	mSizeTree = sizeTree;
	mSizeSearch = 0;

	pthread_mutex_init(&mAddMutex, NULL);
	pthread_mutex_init(&mCntXORMutex, NULL);
	pthread_mutex_init(&mCntSimMutex, NULL);
}

// destructor
QueryResult::~QueryResult() {
	QueryResultNode *node, *nextNode;
	if (mSort) {
		deleteNode(mRootNode);		// discard the whole result set recursivly
	} else {
		// discard the whole result iteratively for large scale results
		node = mRootNode;
		while (node != NULL) {
			if (node->mQueryId != NULL) {
				delete[] node->mQueryId;
			}
			nextNode = node->mRight;
			delete node;
			node = nextNode;
		}
	}
}

// delete a node and all its subnodes recursively
void QueryResult::deleteNode(QueryResultNode *node) {
  if (node != NULL) {
    if (node->mQueryId != NULL) {
      delete[] node->mQueryId;
    }
    if (node->mMatchId != NULL) {
      delete[] node->mMatchId;
    }
    deleteNode(node->mLeft);
    deleteNode(node->mRight);
    delete node;
  }
}

// add a CLK and the corresponding similarity coefficient to the query result
void QueryResult::add(CLK *queryClk, CLK *clk, float similarity) {
	//char resStr[clk->getLength()];
	// lock mutex
	pthread_mutex_lock(&mAddMutex);
		// if not collect results in memory
		QueryResultNode *newNode;
		QueryResultNode **nextNodePtr;
		// create new node
		newNode = new QueryResultNode;

		// copy query id, clk pointer and similarity
		if (queryClk->getId() != NULL) {
			newNode->mQueryId = new char[strlen(queryClk->getId()) + 1];
			strcpy(newNode->mQueryId, queryClk->getId());
		} else {
			newNode->mQueryId = NULL;
		}

		newNode->mPrint = clk;
		newNode->mSimilarity = similarity;
		newNode->mLeft = NULL;
		newNode->mRight = NULL;

		if (mSort) {
			// if query results have to be sorted, find the right location
			// to insert the new node
			nextNodePtr = &mRootNode;

			while (*nextNodePtr != NULL) {
				if ((*nextNodePtr)->mSimilarity >= similarity) {
					nextNodePtr = &((*nextNodePtr)->mRight);
				} else {
					nextNodePtr = &((*nextNodePtr)->mLeft);
				}
			}

			// insert new node
			*nextNodePtr = newNode;
		} else {
			// if query results need not to be sorted
			// insert at end of list
			if (mRootNode == NULL) {
				mRootNode = newNode;
			} else {
				mRootNode->mLeft->mRight = newNode;
			}

			mRootNode->mLeft = newNode; // use mLeft as pointer to end node
		}

	mSize++;

	// unlock mutex
	pthread_mutex_unlock(&mAddMutex);
}

// add a matching pair and the corresponding similarity coefficient to the query result
void QueryResult::add(char *queryId, char *matchId, double similarity) {
  QueryResultNode *newNode;
  QueryResultNode **nextNodePtr;

  // create new node
  newNode = new QueryResultNode;

  // copy query id, match id and similarity
  if (queryId != NULL) {
    newNode->mQueryId = new char[strlen(queryId) + 1];
    strcpy(newNode->mQueryId, queryId);
  } else {
    newNode->mQueryId = NULL;
  }

    newNode->mMatchId = new char[strlen(matchId) + 1];
    strcpy(newNode->mMatchId, matchId);

  newNode->mSimilarity = similarity;
  newNode->mLeft = NULL;
  newNode->mRight = NULL;

  if (mSort) {
    // if query results have to be sorted, find the right location
    // to insert the new node
    nextNodePtr = &mRootNode;

    while (*nextNodePtr != NULL) {
      if ((*nextNodePtr)->mSimilarity >= similarity) {
        nextNodePtr = &((*nextNodePtr)->mRight);
      } else {
        nextNodePtr = &((*nextNodePtr)->mLeft);
      }
    }

    // insert new node
    *nextNodePtr = newNode;
  } else {
    // if query results need not to be sorted
    // insert at end of list
    if (mRootNode == NULL) {
      mRootNode = newNode;
    } else {
      mRootNode->mLeft->mRight = newNode;
    }

    mRootNode->mLeft = newNode; // use mLeft as pointer to end node
  }

  mSize++;
}
