//============================================================================
// Name        : Bloomencoder.cpp
// Author      : Rukasz
// Version     :
// Copyright   : Your copyright notice
// Description : Ansi-style
//============================================================================

#include "BloomHelper.h"

vector<string> CreateBloom(vector<string> vars, unsigned lenBloom, int k,
                         vector<int> lenQgram, vector<int> padding,
                         vector<vector<int> > lookupTable, vector<string> uniBiGrams) {


  vector<string> bloomfilter(vars.size());
  unsigned max_size = 0;
  for (int i = 0; i < (int) vars.size(); i++) {
    if (vars[i].size() > max_size) {
      max_size = vars[i].size();
    }
  }
  vector<string> qgrams(max_size);
  for (int i = 0; i < (int) vars.size(); i++) {
    if (padding[i] > 0)
      Padding(vars[i], padding[i]);

    //create the qgrams from every variable
    qgrams = CreateQgrams(vars[i], lenQgram[i]);

    //create and fill bloomfilter for every variable
    bloomfilter[i] = CreateBloomfilter(qgrams, lenBloom, k, lookupTable,
                                     uniBiGrams, i);
  }

  return bloomfilter;
} // end CreateBloom

void Padding(string &var, int n) {
  for (int i = 0; i < n; i++)
    var = " " + var + " ";
} //end Padding




string CreateEmptyBloomfilter(int lenBloom) {
  string bloomfilter;
  //Create empty bloomfilter
  bloomfilter = "0";
  for (int i = 0; i < (int) (lenBloom - 1); i++) {
    bloomfilter.append("0");
  }
  return bloomfilter;
}

string CreateBloomfilter(vector<string> qgrams, unsigned lenBloom, int k,
                       vector<vector<int> > lookupTable, vector<string> uniBiGrams, int i) {
  //Create empty bloomfilter
  string bloomfilter = CreateEmptyBloomfilter(lenBloom);
  //fill bloomfilter from the qgrams
  for (int l = 0; l < (int) qgrams.size(); l++) {
    //find position of qgram in lookuptable

    for (int m = 0; m < (int) uniBiGrams.size(); m++) {

      if ((qgrams[l].compare(uniBiGrams[m])) == 0) {


        //find the k values in the lookup table for that specific uni-/bigram
        for (int j = 0; j < k; j++) {
          //replace put k "1" in the bloomfilter at the looked up positions, i is the position of the attribute
          bloomfilter.replace(lookupTable[m][k * i + j], 1, "1");
        }
      }
    }
  }


  return bloomfilter;
} //end CreateBloomfilter



vector<string> CreateUniBiGrams(vector<string> uniBiGrams) {
  //inserting unigrams

  for (int i = 0; i < 37; i++) {
    uniBiGrams.push_back(string(ALPHABET, i, 1));
  }
  //inseting bigrams
  for (int i = 0; i < 37; i++) {
    //string b(ALPHABET + i, ALPHABET + 1 + i);
    for (int j = 0; j < 37; j++) {
      //string c(ALPHABET + j, ALPHABET + 1 + j);

      uniBiGrams.push_back(
        (string(ALPHABET, i, 1)).append(string(ALPHABET, j, 1)));
    }
  }

  return uniBiGrams;
}	//end CreateUniBiGrams





vector<vector<int> > LookupTableBloom(int varsSize, string password, int k,
                                      unsigned lenBloom, vector<string> uniBiGrams) {

  //Generating a lookup table for each attribute
  //set seed accourding to the user defined password
  seed_seq seed(password.begin(), password.end());
  //create a random engine
  default_random_engine engine(seed);	// or other engine as std::mt19937
  //vector for uniBiGrams.size() * k random numbers
  vector<int> table(uniBiGrams.size() * k * varsSize);
  //uniform discrete distribution of the random numbers
  // see: http://www.cplusplus.com/reference/random/uniform_real_distribution/
  uniform_int_distribution<int> distr(0, lenBloom - 1);
  //fill vector with uniBiGrams.size() * k random numbers
  generate(begin(table), end(table), [&]() {return distr(engine);});
  //output for testing

  //copy the vector into matrix of size uniBiGrams.size()*k
  vector<vector<int> > lookupTable;
  int counter = 0;

  for (int i = 0; i < (int) uniBiGrams.size(); i++) {
    vector<int> rowUniBiGram;
    for (int j = 0; j < k * varsSize; j++) {
      rowUniBiGram.push_back(table[counter++]);

    }
    lookupTable.push_back(rowUniBiGram);
  }

  return lookupTable;
} //end LookupTableBloom

