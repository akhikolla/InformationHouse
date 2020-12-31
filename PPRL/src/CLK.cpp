// CLK.cpp
//
// Copyright (c) 2015
// Prof. Dr. Rainer Schnell
// Universitaet Duisburg-Essen
// Campus Duisburg
// Institut fuer Soziologie
// Lotharstr. 65
// 47057 Duisburg
//
// This file is part of the R-Package "multibitTree".
//
// "multibitTree" is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// "multibitTree" is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with "multibitTree". If not, see <http://www.gnu.org/licenses/>.



#include "CLK.h"

// static class member for fast cardinality calculation
// array contains pre-computed cardinalities for all possible 16bit-words

int CLK::sCardinalityMap[0x10000]; //0x10000 = 65536

// static class function to initialise cardinality-map

void CLK::init() {
  for (int i = 0; i < 0x10000; i++) {
    sCardinalityMap[i] = 0;
    for (int j = 0; j < 16; j++) {
      if ((i & (1 << j)) != 0) {
        sCardinalityMap[i]++;
      }
    }
  }
}

// convert fingerprint to ascii-string of size n

void CLK::copyToString(char *str, int n) {
  int l;

  l = MIN(mLength, n - 1);
  for (int i = 0; i < l; i++) {
    if (getBit(i) == 0) {
      str[i] ='0';
    } else {
      str[i] ='1';
    }
  }
  str[l] = 0;
}


// convert CLK to int array of size n

void CLK::copyToInt(int* out, int n){
  int l;

  l = MIN(mLength, n);
  for (int i = 0; i < l; i++) {
    if (getBit(i) == 0) {
      out[i]=0;
    } else {
      out[i]=1;
    }
  }
}

// get CLK from string

void CLK::copyFromString(char * id, const char *str)
{
  int l = 0;

  mId = id;

  while ((str[l] == '0') || (str[l] == '1')) {
    l++;
  }
  mLength = l;
  clear();

  for (int i = 0; i < l; i++) {
    if (str[i] != '0') {
      setBit(i);
    }
  }

  //init();
  //fold();
}

// copy

void CLK::copy(CLK *clk) {

  // copy id
  if (clk->mId != NULL) {
    mId = new char[strlen(clk->mId) + 1];
    strcpy(mId, clk->mId);
  } else {
    mId = NULL;
  }

  // copy word-array
  mLength = clk->mLength;

  for (int i = 0; i < arrayLength(); i++) {
    mArray[i] = clk->mArray[i];
  }

  // // copy hash-key
  // for (int i = 0; i < FOLDED_WORDS; i++) {
  //   mHashArray[i] = clk->mHashArray[i];
  // }
}

