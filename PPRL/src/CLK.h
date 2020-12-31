// CLK.h
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

#ifndef CLK_H
#define CLK_H

#include <iostream>
#include <limits>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "Misc.h"
#include <vector>
#include <random> // std::default_random_engine
#include <algorithm> // std::shuffle

using namespace std;

typedef unsigned int WORDTYPE;			// 32bit-words
#define WORD_LEN 32				// word-length in bits
#define BIT1 1u					// unsigned int literal 1

#define FOLDED_WORDS (128/WORD_LEN)		// array-length for hash-key

// Objects of class CLK hold a single bit-vector of arbitrary length
// The class also provides functions for modifying and checking single bits of
// the bit-vector and to compute its cardinality.
// Furthermore there are functions to compute the Tanimoto coefficient of
// two CLKs and the faster 128-bit folding-hash estimation of the Tanimoto
// coefficient.

class CLK {
protected:

  char *mId;				// ID string
  WORDTYPE *mArray;			// array for stored bits
  WORDTYPE mHashArray[FOLDED_WORDS];	// 128 Bit folded Hash-Key
  int mLength;				// length of CLK in bits
  int mOrigin;				// origin, eg fileA=1, fileB=2
  int mCardinality;			// cardinality buffer of CLK
  bool mInP;				// indicator for set P (CanopyCLustering)
  static int sCardinalityMap[0x10000];	// static 16-bit cardinality-map
  int wLength; 	// length of array of weights for weigthed RLBF
  vector<int> partsize; 	// vector int for partsizes for weigthed RLBF
  int sum = 0;
  int maxPartSize = 0;
  // calculate length of word-array

  inline int arrayLength() {
    return (mLength - 1) / WORD_LEN + 1;
  }

  // initialize word-array

  inline void allocate() {
    mLength = MAX(mLength, 128);
    mArray = new WORDTYPE[arrayLength()];
  }

  // calculate cardinality of 32bit-word

  inline int cardWord(WORDTYPE word) {
    return sCardinalityMap[word & 0xFFFF] + sCardinalityMap[(word >> 16) & 0xFFFF];
  }

public:
  // constructor for empty CLK

  inline CLK() {
    mId = NULL;
    wLength = 0;
    allocate();
    clear();
    fold();

  }

  // constructor for empty CLK of given bit-length

  inline CLK(int l, int origin) {
    mId = NULL;
    mLength = l;
    mOrigin = origin;
    mInP = true;
    allocate();
    clear();
  }

  // constructor for CLK based on given ascii-string

  inline CLK(char * id, const char *str, int origin) {
    int l = 0;

    mId = id;
    mOrigin = origin;

    while ((str[l] == '0') || (str[l] == '1')) {
      l++;
    }

    mLength = l;
    mInP = true;

    allocate();
    clear();

    for (int i = 0; i < l; i++) {
      if (str[i] != '0') {
        setBit(i);
        mCardinality++;
      }
    }
  }

  // constructor for empty CLK of given bit-length

  inline CLK(int length) {
    mId = NULL;
    mLength = length;
    allocate();
    clear();
    fold();
  }

  // constructor for empty CLK of given id and ascci-string

  inline CLK(char * id, const char *str) {
    int l = 0;

    mId = id;

    while ((str[l] == '0') || (str[l] == '1')) {
      l++;
    }

    mLength = l;

    allocate();
    clear();

    for (int i = 0; i < l; i++) {
      if (str[i] != '0') {
        setBit(i);
      }
    }

    fold();
  }

  // copy-constructor for CLKs

  inline CLK(CLK *clk) {

    // copy id
    if (clk->mId != NULL) {
      mId = new char[strlen(clk->mId) + 1];
      strcpy(mId, clk->mId);
    } else {
      mId = NULL;
    }

    // copy word-array
    mLength = clk->mLength;
    allocate();

    for (int i = 0; i < arrayLength(); i++) {
      mArray[i] = clk->mArray[i];
    }

    // copy hash-key
    for (int i = 0; i < FOLDED_WORDS; i++) {
      mHashArray[i] = clk->mHashArray[i];
    }
  }

  // destructor
  inline ~CLK() {
    if (mId != NULL) {
      delete[] mId;
    }
    delete[] mArray;
  }

  // get id

  inline char *getId() {
    return mId;
  }

  // set id

  inline void setId(char *id) {
    mId = id;
  }

  // get origin

  inline int getOrigin() {
    return mOrigin;
  }

  // clear all bits

  inline void clear() {
    for (int i = 0; i < arrayLength(); i++) {
      mArray[i] = 0;
    }
    mCardinality = 0;
  }

  // check if all bits are cleared

  inline int isEmpty() {
    for (int i = 0; i < arrayLength(); i++) {
      if (mArray[i] != 0) {
        return 0;
      }
    }

    return 1;
  }

  // checks if two clks are equal
  inline int isEqual(CLK *clk) {
    for (int i = 0; i < arrayLength(); i++) {
      if (mArray[i] != clk->mArray[i]) {
        return 0;
      }
    }

    return 1;
  }

  // set length in bits
  inline void setLength(int length) {
    mLength = length;
    clear();
    fold();
  }

  // get length in bits

  inline int getLength() {
    return mLength;
  }

  inline bool isInP() {
    return mInP;
  }

  inline void removeFromP() {
    mInP = false;
  }

  // compute cardinality of intersection of two clks

  inline int intersectCard(CLK *clk) {
    int count = 0;

    for (int i = 0; i < arrayLength(); i++) {
      count += cardWord(mArray[i] & clk->mArray[i]);
    }

    return count;
  }

  // compute hash-key

  inline void fold() {
    int len;

    len = arrayLength();

    for (int i = 0; i < FOLDED_WORDS; i++) {
      mHashArray[i] = mArray[i];
    }

    for (int i = FOLDED_WORDS; i < len; i++) {
      mHashArray[i % FOLDED_WORDS] ^= mArray[i];
    }
  }

  // set bit at position n

  inline void setBit(int n) {
    mArray[n / WORD_LEN] |= (BIT1 << (n % WORD_LEN));
  }

  // get bit at position n

  inline WORDTYPE getBit(int n) {
    return (mArray[n / WORD_LEN] >> (n % WORD_LEN)) & BIT1; //in mArray gibt es 32 bit lange "Wörter" mit Nullen und Einsen, die Position in einem mArray-Ellement wird mit BIT1 =1 UND-verknüpft
  }

  // unset bit at position n

  inline void unsetBit(int n) {
    mArray[n / WORD_LEN] &= 0xFFFF - (BIT1 << (n % WORD_LEN));
  }

  // join one bits of second clk
  inline void join(CLK *clk) {
    for (int i = 0; i < arrayLength(); i++) {
      mArray[i] = mArray[i] | clk->mArray[i];
    }
  }

  // compute tanimoto-index

  inline double tanimoto(CLK *clk) {
    int count_and = 0;
    int count_or = 0;
    int min, len;

    len = arrayLength();
    min = clk->arrayLength();

    // handle different bit-length
    if (min <= len) {
      for (int i = min; i < len; i++) {
        count_or += cardWord(mArray[i]);
      }
    } else {
      for (int i = len; i < min; i++) {
        count_or += cardWord(clk->mArray[i]);
      }
      min = len;
    }

    for (int i = 0; i < min; i++) {
      int a = mArray[i] & clk->mArray[i];
      int o = mArray[i] | clk->mArray[i];
      count_and += cardWord(a);
      count_or += cardWord(o);
    }

    return ((double) count_and) / count_or;
  }

  // generic similarity estimation

  inline double estimatedSimilarity(CLK *clk) {
    return tanimoto(clk);
  }

  //Xor

  inline void ClkXOR(CLK *clk) {
    for (int i = 0; i < arrayLength(); i++) {
      mArray[i] = mArray[i] ^ clk->mArray[i];
    }
  }

  // permutes CLK
  inline void permuteCLK(string password, int n) {
    // copy word-array
    int l = MIN(mLength, n);
    seed_seq seed(password.begin(), password.end());
    int *a = new int[l];
    for (int i = 0; i < l; i++) {
      if (getBit(i) == 0) {
        a[i] = 0;
      } else {
        a[i] = 1;
      }
    }
    shuffle(&a[0], &a[l], default_random_engine(seed));
    for (int i = 0; i < l; i++) {
      if (a[i] == 1) {
        setBit(i);
      }
    }
    delete[] a;
  }

  // set size of the parts of CLK
  inline void setPartsize(int lenWeights, float *weigth) {
    wLength = lenWeights;
    partsize.clear();
    for (int w = 0; w < wLength; w++) {
      partsize.push_back((int) floor(mLength * weigth[w]));
      sum += partsize[w];
    }
    if (sum < mLength) {
      for (int i = 0; i < mLength - sum; i++) {
        partsize[i] += 1;
      }
    }
    maxPartSize = *max_element(partsize.begin(), partsize.end());
  }

  // size of the parts of CLK
  inline vector<int> getPartsize() {
    return partsize;
  }

  // maximal size of the parts of CLK
  inline int getMaxPartsize() {
    return maxPartSize;
  }

  // negate CLK
  inline void negate() {
    for (int i = 0; i < mLength / WORD_LEN; i++)
      mArray[i] = ~mArray[i];
  }

  // count set bits

  inline int cardinality() {
    int count = 0;

    for (int i = 0; i < arrayLength(); i++) {
      count += cardWord(mArray[i]);
    }

    return count;
  }

  // generic sortkey calculation

  inline int getCardinality() {
    return(mCardinality);
  }


  // static class function to initialise cardinality-map

  static void init();

  // convert CLK to ascii-string of size n

  void copyToString(char *str, int n);


  // convert CLK to int array of size n
  void copyToInt(int* out, int n);

  void copyFromString(char * id, const char *str);
  void copy(CLK *clk);
};

class CLKXOR: public CLK {
protected:

  WORDTYPE mHashArray[FOLDED_WORDS];	// 128 Bit folded Hash-Key

public:

  inline CLKXOR(int length, int origin) :
    CLK(length, origin) {
    fold();
  }

  inline CLKXOR(char * id, const char *str, int origin) :
    CLK(id, str, origin) {
    fold();
  }
  // constructor for CLK based on given ascii-string

  inline CLKXOR(char * id, const char *str) :
    CLK(id, str) {
    fold();
  }

  // compute hash-key

  inline void fold() {
    int len;

    len = arrayLength();

    for (int i = 0; i < FOLDED_WORDS; i++) {
      mHashArray[i] = mArray[i];
    }

    for (int i = FOLDED_WORDS; i < len; i++) {
      mHashArray[i % FOLDED_WORDS] ^= mArray[i];
    }
  }

  // compute tanimoto estimation on hash-keys

  inline double tanimotoXOR(CLKXOR *clk, int AB) {

    int xorCount = cardWord(mHashArray[0] ^ clk->mHashArray[0]) + cardWord(mHashArray[1] ^ clk->mHashArray[1])
    + cardWord(mHashArray[2] ^ clk->mHashArray[2]) + cardWord(mHashArray[3] ^ clk->mHashArray[3]);

    return ((double) (AB - xorCount)) / (AB + xorCount);
  }

  // generic similarity estimation

  inline double estimatedSimilarity(CLKXOR *clk) {
    return tanimotoXOR(clk, mCardinality + clk->mCardinality);
  }
};

// specialized CLK class for computing Tanimoto similarity
class CLKTanimoto: public CLK {
public:

  inline CLKTanimoto(int length) :
    CLK(length) {
  }
  inline CLKTanimoto(char * id, const char *str) :
    CLK(id, str) {
  }
  inline CLKTanimoto(CLKTanimoto *clk) :
    CLK((CLK*) clk) {
  }
  inline int isEqual(CLK *clk) {
    return CLK::isEqual((CLK*) clk);
  }
  inline void join(CLK *clk) {
    CLK::join((CLK*) clk);
  }
  inline int intersectCard(CLK *clk) {
    return CLK::intersectCard((CLK*) clk);
  }

  typedef float S;

  // compute tanimoto similarity
  inline S tanimoto(CLKTanimoto *clk) {
    int count_and = 0;
    int count_or = 0;
    int min, len;

    len = arrayLength();
    min = clk->arrayLength();

    // handle different bit-length
    if (min <= len) {
      for (int i = min; i < len; i++) {
        count_or += cardWord(mArray[i]);
      }
    } else {
      for (int i = len; i < min; i++) {
        count_or += cardWord(clk->mArray[i]);
      }
      min = len;
    }

    for (int i = 0; i < min; i++) {
      int a = mArray[i] & clk->mArray[i];
      int o = mArray[i] | clk->mArray[i];
      count_and += cardWord(a);
      count_or += cardWord(o);
    }
    return ((S) count_and) / count_or;
  }

  // compute tanimoto estimation on hash-keys
  inline S tanimotoXOR(CLKTanimoto *clk, int AB) {

    int xorCount = cardWord(mHashArray[0] ^ clk->mHashArray[0]) + cardWord(mHashArray[1] ^ clk->mHashArray[1])
    + cardWord(mHashArray[2] ^ clk->mHashArray[2]) + cardWord(mHashArray[3] ^ clk->mHashArray[3]);

    return ((S) (AB - xorCount)) / (AB + xorCount);
  }

  // compute lower bound for 1D-Grid
  static inline int lowerBound(S minTanimotoSimilarity, int card) {
    return (int) ceil(minTanimotoSimilarity * card);
  }

  // compute upper bound for 1D-Grid
  static inline int upperBound(S minTanimotoSimilarity, int card) {
    return (int) (1.0 / minTanimotoSimilarity * card) + 1;
  }
};

// specialized CLK class for computing Hamming distance
class CLKHamming: public CLK {
public:

  inline CLKHamming(int length) :
    CLK(length) {
  }
  inline CLKHamming(char * id, const char *str) :
    CLK(id, str) {
  }
  inline CLKHamming(CLKHamming *clk) :
    CLK((CLK*) clk) {
  }
  inline int isEqual(CLK *clk) {
    return CLK::isEqual((CLK*) clk);
  }
  inline void join(CLK *clk) {
    CLK::join((CLK*) clk);
  }
  inline int intersectCard(CLK *clk) {
    return CLK::intersectCard((CLK*) clk);
  }

  typedef int S;

  // compute hamming distance
  inline S hamming(CLKHamming *clk) {
    int count_xor = 0;
    int min, len;

    len = arrayLength();
    min = clk->arrayLength();

    // handle different bit-length
    if (min <= len) {
      for (int i = min; i < len; i++) {
        count_xor += cardWord(mArray[i]);
      }
    } else {
      for (int i = len; i < min; i++) {
        count_xor += cardWord(clk->mArray[i]);
      }
      min = len;
    }

    for (int i = 0; i < min; i++) {
      count_xor += cardWord(mArray[i] ^ clk->mArray[i]);
    }

    return count_xor;
  }

  // compute hamming distance estimation on hash-keys
  inline S hammingXOR(CLKHamming *clk) {

    return (cardWord(mHashArray[0] ^ clk->mHashArray[0]) + cardWord(mHashArray[1] ^ clk->mHashArray[1])
              + cardWord(mHashArray[2] ^ clk->mHashArray[2]) + cardWord(mHashArray[3] ^ clk->mHashArray[3]));
  }

  // compute lower bound for 1D-Grid
  static inline int lowerBound(S maxHammingDistance, int card) {
    return card - maxHammingDistance;
  }

  // compute upper bound for 1D-Grid
  static inline int upperBound(S maxHammingDistance, int card) {
    return card + maxHammingDistance;
  }
};

#endif
