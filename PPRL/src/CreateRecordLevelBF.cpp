//============================================================================
// Name        : CreateRecordLevelBF.cpp
// Author      : Rukasz
// Version     :
// Copyright   : Your copyright notice
// Description : Ansi-style
//============================================================================

#include "CreateRecordLevelBF.h"

using namespace std;

/*creates RLBs as proposed by Elizabeth Ashley Durham in "A FRAMEWORK FOR ACCURATE, EFFICIENT PRIVATE RECORD
LINKAGE"
1rst method Static/Uniform RBFs page 65ff and 73
*/
void CreateRBF(CLK* rlbf, vector<string> vars, int k, vector<int> padding,
             vector<int> lenQgram, unsigned lenBloom, int lenRLBF, vector<string> password) {
  vector<string> bf(vars.size());
  // if partsize is smaller than the BF from which the bits are taken, see Durham Disesrtation p 90 bottom
  int len = rlbf->getMaxPartsize();
  len = max(len, (int)lenBloom);
  //determine the samples of the bloomfilters

  //determine size of the parts of the RLBF
  vector<int> partsize = rlbf->getPartsize();
  //create single Bloom filter
  for (unsigned i = 0; i < vars.size(); i++) {
    bf[i] = CreateBFBigramSeed(vars[i],  k, padding[i],
                              lenQgram[i], lenBloom, password[i]);
  }
  //vector containing the numbers 0 to length of the bloomfilter
  vector<int> numbers(len);
  for (int j = 0; j < len; j++) {
    numbers[j]=j;
  }
  //take samples of every bloomfilter
  int start = 0; //to determine the posotion in the rlbf
  int n; //
  string password2;
  for (unsigned i = 0; i < vars.size(); i++) {
    n=bf[i].size();
    //shuffle(numbers.begin(), numbers.begin()+n, default_random_engine(seed));//shuffle numbers for each BF
    for (int m = 0; m < partsize[i]; m++) {
      if (bf[i][numbers[m]%n]=='1'){// using modulo on the BF in case partsize is bigger than the BF from which the bits are taken, see Durham Disesrtation p 90 bottom
        rlbf->setBit(m + start);
      }
    }
    start += partsize[i];
    password2 = password2+password[i];
  }

  rlbf->permuteCLK(password2,lenRLBF);
}		//end CreateRLBFUniformBitSelection

