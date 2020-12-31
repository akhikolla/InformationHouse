/*
 * MTB_Result.h
 *
 *  Created on: 02.08.2017
 *      Author: rukasz
 */

#ifndef MTB_RESULT_H_
#define MTB_RESULT_H_

#include <string>
#include <vector>
#include <iostream>

using namespace std;

/**
* Holds result data
*/

class MTB_Result {
public:
  /**
  * constructor
  */
  MTB_Result() {

  }
  ;
  ~MTB_Result() {
  }
  ;

  void addResult(string data1, string ID1, string blockingData1, string data2,
                 string ID2, string blockingData2, float res) {
    this->data1V.push_back(data1);
    this->ID1V.push_back(ID1);
    this->blockingData1V.push_back(blockingData1);
    this->data2V.push_back(data2);
    this->ID2V.push_back(ID2);
    this->blockingData2V.push_back(blockingData2);
    this->res.push_back(res);
  }

  void addResult(string data1, string ID1, string data2,
                 string ID2, float res) {
    this->data1V.push_back(data1);
    this->ID1V.push_back(ID1);
    this->data2V.push_back(data2);
    this->ID2V.push_back(ID2);
    this->res.push_back(res);
  }

  void addResult(string data1, string ID1, string blockingData1, string data2,
                 string ID2, string blockingData2, float res, string match) {
    this->data1V.push_back(data1);
    this->ID1V.push_back(ID1);
    this->blockingData1V.push_back(blockingData1);
    this->data2V.push_back(data2);
    this->ID2V.push_back(ID2);
    this->blockingData2V.push_back(blockingData2);
    this->res.push_back(res);
    this->match.push_back(match);
  }

  void addResult(string ID1, string ID2,  float res) {
     this->ID1V.push_back(ID1);
     this->ID2V.push_back(ID2);
     this->res.push_back(res);
   }

  vector<string> getdata1() {
    return this->data1V;
  }

  vector<string> getID1() {
    return this->ID1V;
  }

  vector<string> getdata2() {
    return this->data2V;
  }

  vector<string> getID2() {
    return this->ID2V;
  }

  vector<float> getRes() {
    return this->res;
  }

  void setRes(vector<float> resIn) {
    this->res = resIn;
  }

  vector<string> getblockingData1() {
    return this->blockingData1V;
  }

  vector<string> getblockingData2() {
    return this->blockingData2V;
  }

  vector<string> getMatch() {
    return this->match;
  }

private:
  vector<string> data1V, ID1V, data2V, ID2V, blockingData1V, blockingData2V;
  vector<float> res;
  vector<string> match;
  vector<int> card1, card2;
};

#endif /* MTB_RESULT_H_ */
