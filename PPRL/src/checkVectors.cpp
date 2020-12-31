/*
 * checkVectors.cpp
 *
 *  Created on: 13.09.2016
 *      Author: schnell-42
 */

#include "checkVectors.h"

void checkVectors(unsigned varsSize,  vector<int>* padding,  vector<int> *qgram){
  //Sizes of Vectors are checked and adapted
  if (varsSize>padding->size()){
    Rcpp::Rcerr << "Vector padding must have the same size as the input data.frame. Padding will be fill with zeros.\n";
    for (unsigned i = padding->size() ; i< varsSize; i++){
      padding->push_back(0);
    }
    //return NULL; Alternative zu Auffüllen ist Programm abzubrechen.
  }
  if (varsSize<padding->size()){
    Rcpp::Rcerr << "Vector padding must have the same size as the input data.frame. Padding will be cut.\n";
    while((unsigned int)padding->size()> varsSize){
      padding->pop_back();
    }
    //return NULL; Alternative zu Auffüllen ist Programm abzubrechen.
  }
  if (varsSize>qgram->size()){
    Rcpp::Rcerr << "Vector qgrams must have the same size as the input data.frame. Qgrams will be fill with 2s.\n" ;
    for (unsigned i = qgram->size() ; i< varsSize; i++){
      qgram->push_back(2);
    }
    //return NULL; Alternative zu Auffüllen ist Programm abzubrechen.
  }
  if (varsSize<qgram->size()){
    Rcpp::Rcerr << "Vector qgram must have the same size as the input data.frame. Qgram will be cut.\n";
    while((unsigned int)qgram->size()> varsSize){
      qgram->pop_back();
    }
    //return NULL; Alternative zu Auffüllen ist Programm abzubrechen.
  }
}
