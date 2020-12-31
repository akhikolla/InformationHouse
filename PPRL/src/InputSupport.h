// InputSupport.h
//
// Copyright (c) 2017
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

#ifndef INPUTSUPPORT_H
#define INPUTSUPPORT_H

#define __STDC_FORMAT_MACROS
#include <inttypes.h>
#include <vector>
//#include "Grid1D.h"

// maximal line size = length of ascii representation of CLK
#define STRSIZE 4000
using namespace std;

// objects of template class FileSupport provide functionalities to
// read CLK files in csv like format

template<class T>
class InputSupport {
private:

  // check, if a character is considered to be a white-space or seperator
  inline int isWS(char c) {
    return ((c == '"') || (c == '\'') || (c == ',') || (c == ';') || (c == ' ') || (c == '\t'));
  }

  // check, if a character is considered to be end of line
  inline int isEOL(char c) {
    return ((c == 10) || (c == 13) || (c == 0));
  }

  int parseLine(FILE *in, char *str, int64_t *idx1, int64_t *end1, int64_t *idx2, int64_t *end2);

public:

  T* parseCLK(string IDc, string CLKin);
  T* parseCLK(string IDc, string CLKin, int origin);
  void parseAllCLKs(vector<string> IDc, vector<string> CLKin, T ***clks, int64_t *size, int *nBits);
  void parseAllCLKs(vector<string> IDc, vector<string> CLKin, T ***clks, int64_t *size, int *nBits, int origin);
};

// create a CLK
template<class T>
T* InputSupport<T>::parseCLK(string IDc, string CLKin) {
  char *idStr = new char[13];
  char CLKstr[STRSIZE];

  strcpy(CLKstr, CLKin.c_str());
  strcpy(idStr, IDc.c_str());
  return new T(idStr, CLKstr);
}

// create a CLK
template<class T>
T* InputSupport<T>::parseCLK(string IDc, string CLKin, int origin) {
  char *idStr = new char[13];
  char CLKstr[STRSIZE];

  strcpy(CLKstr, CLKin.c_str());
  strcpy(idStr, IDc.c_str());
  return new T(idStr, CLKstr, origin);
}

// read CLK array
template<class T>
void InputSupport<T>::parseAllCLKs(vector<string> IDc, vector<string> CLKin, T ***clks, int64_t *size, int *nBits) {
  int64_t sizeCLKs;
  int maxLen;
  int len;
  T **readCLKs;

  sizeCLKs = CLKin.size(); //number of CLKs
  maxLen = 0;

  if (sizeCLKs == 0) {
    *size = 0;
    *nBits = 0;
    return;
  }

  // create array of CLKs
  readCLKs = new T*[sizeCLKs];

  // read CLKs
  for (int64_t i = 0; i < sizeCLKs; i++) {
    // parse line
    readCLKs[i] = parseCLK(IDc[i], CLKin[i]);
   //Rcpp::Rcout << "parseAllCLKs " << readCLKs[i]->getId() << endl; hier ok
    // compute maximal length of CLKs
    len = readCLKs[i]->getLength();
    if (len > maxLen) {
      maxLen = len;
    }
  }

  *nBits = maxLen;
  *size = sizeCLKs;
  *clks = readCLKs;
}

// read CLK array
template<class T>
void InputSupport<T>::parseAllCLKs(vector<string> IDc, vector<string> CLKin, T ***clks, int64_t *size, int *nBits, int origin) {
  int64_t sizeCLKs;
  int maxLen;
  int len;
  int64_t sizeA = *size;
  T **readCLKs;

  sizeCLKs = CLKin.size(); //number of CLKs
  // second pass: count number of CLKs
  maxLen = 0;

  if (sizeCLKs == 0) {
    *size = 0;
    *nBits = 0;
    return;
  }

  // create array of CLKs
  readCLKs = new T*[sizeCLKs];

  // read CLKs
  for (int64_t i = 0; i < sizeA; i++) {
    // parse line
    readCLKs[i] = parseCLK(IDc[i], CLKin[i], 1);

    if (readCLKs[i] == NULL) {
      sizeCLKs = i;
      break;
    }

    // compute maximal length of CLKs
    len = readCLKs[i]->getLength();
    if (len > maxLen) {
      maxLen = len;
    }
  }

  // read CLKs
  for (int64_t i = sizeA; i < sizeCLKs; i++) {
    // parse line
    readCLKs[i] = parseCLK(IDc[i], CLKin[i], 2);

    if (readCLKs[i] == NULL) {
      sizeCLKs = i;
      break;
    }

    // compute maximal length of CLKs
    len = readCLKs[i]->getLength();
    if (len > maxLen) {
      maxLen = len;
    }
  }

  *nBits = maxLen;
  *size = sizeCLKs;
  *clks = readCLKs;
}
#endif
