//============================================================================
// Name        : CreateCLKBigramSeed.cpp
// Author      : Rukasz
// Version     :
// Copyright   : Your copyright notice
// Description : Ansi-style
//============================================================================

#include "CreateCLKBigramSeed.h"
#include <Rcpp.h>
using namespace std;

string CreateBFBigramSeed(string vars, int k, int padding,
                      int lenQgram, unsigned lenBloom, string password) {

  //Standardisation
  replaceNonAscii(vars);
  if (padding > 0)
    Padding(vars, padding);

  vector<string> qgrams= CreateQgrams(vars, lenQgram);
  char *CLKout = new char[lenBloom+1];
  //init empty CLK
  CLK* bf = new CLK(lenBloom);
  bf->init();
  for (int l = 0; l < (int) qgrams.size(); l++) {
    //generate password from  first bigram (before padding), password.
    string password2;
    password2 = qgrams[l] + password;
    //generate seed with password above
    //set seed accourding to the user defined password
    seed_seq seed(password2.begin(), password2.end());
    //create a random engine
    default_random_engine engine(seed);	// or other engine as std::mt19937
    //vector for k random numbers
    vector<int> table(k);
    //uniform discrete distribution of the random numbers
    // see: http://www.cplusplus.com/reference/random/uniform_real_distribution/
    uniform_int_distribution<int> distr(0, lenBloom - 1);
    //fill vector with k random numbers
    generate(begin(table), end(table), [&]() {return distr(engine);});

    //find the k values in the lookup table for that specific uni-/bigram
    for (int j = 0; j < k; j++) {
      // put k "1" in the CLK
      bf->setBit(table[j]);
    }

  }
   bf->copyToString (CLKout, lenBloom+11);
    string CLKouts = CLKout;
    delete[] CLKout;
    delete bf;
   return CLKouts;
} //end CreateBFBigramSeed

string CreateCLKBigramSeed(vector<string> vars, int k, vector<int> padding,
                         vector<int> lenQgram, unsigned lenBloom, vector<string> password) {

  //Standardisation
  unsigned max_size = 0;
  for (unsigned i = 0; i < vars.size(); i++) {
    replaceNonAscii(vars[i]);
    if (vars[i].size() > max_size) {
      max_size = vars[i].size();
    }
  }

  char *CLKout = new char[lenBloom+1];
  //init empty CLK
  CLK* clk = new CLK(lenBloom);
  clk->init();
  //for every variable
  for (unsigned i = 0; i < vars.size(); i++) {
    vector<string> qgrams;
    if (padding[i] > 0)
      Padding(vars[i], padding[i]);

    qgrams = CreateQgrams(vars[i], lenQgram[i]);
    for (int l = 0; l < (int) qgrams.size(); l++) {
      //generate password from  first bigram (before padding), password and column no.of the input vector.
      string password2;
      password2 = qgrams[l] + password[i];
      //generate seed with password above
      //set seed accourding to the user defined password
      seed_seq seed(password2.begin(), password2.end());
      //create a random engine
      default_random_engine engine(seed);	// or other engine as std::mt19937
      //vector for k random numbers
      vector<int> table(k);
      //uniform discrete distribution of the random numbers
      // see: http://www.cplusplus.com/reference/random/uniform_real_distribution/
      uniform_int_distribution<int> distr(0, lenBloom - 1);
      //fill vector with k random numbers
      generate(begin(table), end(table), [&]() {return distr(engine);});

      for (int j = 0; j < k; j++) {
        // put k "1" in the CLK
        clk->setBit(table[j]);
      }
    }
  }
  clk->copyToString(CLKout,lenBloom+1);
  string CLKouts = CLKout;
  delete[] CLKout;
  delete clk;
  return CLKouts;
} //end CreateCLKBigramSeed

string CreateMarkovCLKc(vector<string> vars, int k1, int k2, vector<int> padding,
                     vector<int> lenQgram, unsigned lenBloom, vector<string> password,
                     vector<vector<float>> markovTable, vector<string> rownames, vector<string> colnames, bool includeOriginalBigram, bool v)
  {

  //Standardisation
  unsigned max_size = 0;// max length of input strings
  for (unsigned i = 0; i < vars.size(); i++) {
    replaceNonAscii(vars[i]);
    if (vars[i].size() > max_size) {
      max_size = vars[i].size();
    }
  }

  // iterators to find qgrams in markovTable
  vector<string>::iterator itrow;
  vector<string>::iterator itcol;
  vector<double> qgram_weights;
  vector<string> extraQgrams;

  //draw k2 new qgrams to each qgram

  char *CLKout = new char[lenBloom+1];
  //init empty CLK
  CLK* clk = new CLK(lenBloom);
  clk->init();
  //for every variable
  for (unsigned i = 0; i < vars.size(); i++) {
    vector<string> qgrams;
    if (padding[i] > 0)
      Padding(vars[i], padding[i]);
    //create qgrams
    qgrams = CreateQgrams(vars[i], lenQgram[i]);
    for (int l = 0; l < (int) qgrams.size()-1; l++) {
      //generate password from  first bigram (before padding), password and column no.of the input vector.

      //draw k2 new qgrams to each qgram and add to vector of qgrams

       string password2;
      // // Übergangwert aus markovTabel wird in seed aufgenommen
      // // Erstelle Liste aus Übergängen
      // // suche in Markov table den Übergangswert

      itrow=find(rownames.begin(),rownames.end(),qgrams[l]);
      itcol=find(colnames.begin(),colnames.end(),qgrams[l+1]);
      //get the vector of weights to qgrams[l]
      qgram_weights.clear();
      password2 = qgrams[l] + password[i];
      if(itrow != rownames.end() &&itcol != colnames.end() ){
        for (unsigned j = 0; j < colnames.size(); j++){
          qgram_weights.push_back(markovTable[(itrow-rownames.begin())][j]);
        }
       //generate k2 additional qgrams
        extraQgrams = pickWeightedRandomQgrams(colnames, qgram_weights, k2, password2);
        if(includeOriginalBigram){
        extraQgrams.insert(extraQgrams.begin(), qgrams[l]);
          }
      }
      //for each of these extra Qgram bits are drawn
      for (unsigned k = 0; k < extraQgrams.size(); k++){
      password2 = extraQgrams[k] + password[i];// + offset;
      //für christian
      if (v)
        Rcpp::Rcout << extraQgrams[k] << " ";
      //generate seed with password above
      //set seed accourding to the user defined password
      seed_seq seed(password2.begin(), password2.end());
      //create a random engine
      default_random_engine engine(seed);	// or other engine as std::mt19937
      //vector for k random numbers
      vector<int> table(k1);
      //uniform discrete distribution of the random numbers
      // see: http://www.cplusplus.com/reference/random/uniform_real_distribution/
      uniform_int_distribution<int> distr(0, lenBloom - 1);
      //fill vector with k random numbers
      generate(begin(table), end(table), [&]() {return distr(engine);});

      //find the k values in the lookup table for that specific uni-/bigram
      for (int j = 0; j < k1; j++) {
        // put k "1" in the CLK
        clk->setBit(table[j]);
      }
    }
      //für christian
      if (v)
        Rcpp::Rcout << endl;
    }
  }
  clk->copyToString(CLKout, lenBloom+1);
  string CLKouts = CLKout;
  delete[] CLKout;
  delete clk;
  return CLKouts;
} //end CreateCLKMarkovc


string CreateMarkovCLKOldc(vector<string> vars, int k, vector<int> padding,
                        vector<int> lenQgram, unsigned lenBloom, vector<string> password,
                        vector<vector<float>> markovTable, vector<string> rownames, vector<string> colnames) {
  //Standardisation
  unsigned max_size = 0;// max length of input strings
  for (unsigned i = 0; i < vars.size(); i++) {
    replaceNonAscii(vars[i]);
    if (vars[i].size() > max_size) {
      max_size = vars[i].size();
    }
  }
  // iterators to find qgrams in markovTable
  vector<string>::iterator itrow;
  vector<string>::iterator itcol;
  char *CLKout = new char[lenBloom+1];
  //init empty CLK
  CLK* clk = new CLK(lenBloom);
  clk->init();
  //for every variable
  for (unsigned i = 0; i < vars.size(); i++) {
    vector<string> qgrams;
    if (padding[i] > 0)
      Padding(vars[i], padding[i]);

    qgrams = CreateQgrams(vars[i], lenQgram[i]);
    for (int l = 0; l < (int) qgrams.size()-1; l++) {
      //generate password from  first bigram (before padding), password and column no.of the input vector.
      string password2;
      stringstream ss;
      stringstream ssm;
      ss << i;
      string offset = ss.str();
      string markov;
      // Übergangwert aus markovTabel wird in seed aufgenommen
      // Erstelle Liste aus Übergängen
      // suche in Markov tabe den Übergangswert

      itrow=find(rownames.begin(),rownames.end(),qgrams[l]);
      itcol=find(colnames.begin(),colnames.end(),qgrams[l+1]);

      if(itrow != rownames.end() &&itcol != colnames.end() ){
        ssm << markovTable[(itrow-rownames.begin())][(itcol-colnames.begin())];//get position in markovTable
        markov = ssm.str();

      }
      password2 = markov + password[i];// + offset;
      //generate seed with password above
      //set seed accourding to the user defined password
      seed_seq seed(password2.begin(), password2.end());
      //create a random engine
      default_random_engine engine(seed);	// or other engine as std::mt19937
      //vector for k random numbers
      vector<int> table(k);
      //uniform discrete distribution of the random numbers
      // see: http://www.cplusplus.com/reference/random/uniform_real_distribution/
      uniform_int_distribution<int> distr(0, lenBloom - 1);
      //fill vector with k random numbers
      generate(begin(table), end(table), [&]() {return distr(engine);});
      //find the k values in the lookup table for that specific uni-/bigram
      for (int j = 0; j < k; j++) {
        // put k "1" in the CLK
        clk->setBit(table[j]);
      }
    }
  }
  clk->copyToString(CLKout, lenBloom+1);
  string CLKouts = CLKout;
  delete[] CLKout;
  return CLKouts;
} //end CreateCLKMarkovOld

string CreateEnsembleCLKc(vector<string> vars, int k, int NumberOfCLK, vector<int> padding,
                       vector<int> lenQgram, unsigned lenBloom, vector<string> password) {
  //Standardisation
  unsigned max_size = 0;
  for (unsigned i = 0; i < vars.size(); i++) {
    replaceNonAscii(vars[i]);
    if (vars[i].size() > max_size) {
      max_size = vars[i].size();
    }
  }
  //vector<string> qgrams(max_size);
  //final clk
  CLK* clkres = new CLK(lenBloom);
  clkres->init();
  char *CLKout = new char[lenBloom+1];

  //empty vector of CLKs for majority calculation
  vector<CLK*> clk(NumberOfCLK);
  for (int n = 0; n < NumberOfCLK; n++){
    clk[n]= new CLK(lenBloom);
    clk[n]->init();
  }
  //vector int for to count set bits
  int *majority = new int[lenBloom];
  for (int i = 0; i < lenBloom; i++){
    majority[i] = 0;
  }
  //memset(majority, 0, lenBloom);
  //for every variable
  for (unsigned i = 0; i < vars.size(); i++) {
    vector<string> qgrams;
    if (padding[i] > 0)
      Padding(vars[i], padding[i]);

    qgrams = CreateQgrams(vars[i], lenQgram[i]);
    for (int l = 0; l < (int) qgrams.size(); l++) {
      //generate password from  first bigram (before padding), password and column no.of the input vector.
      stringstream ss;
      ss << i;
      string offset = ss.str();
      string password2;
      password2 = qgrams[l] + password[i] + offset;
      //generate seed with password above
      //set seed accourding to the user defined password
      seed_seq seed(password2.begin(), password2.end());
      //create a random engine
      default_random_engine engine(seed);	// or other engine as std::mt19937
      //vector for k random numbers
      vector<int> table(k*NumberOfCLK);
      //uniform discrete distribution of the random numbers
      // see: http://www.cplusplus.com/reference/random/uniform_real_distribution/
      uniform_int_distribution<int> distr(0, lenBloom - 1);
      //fill vector with k random numbers
      generate(begin(table), end(table), [&]() {return distr(engine);});

      for (int n = 0; n < NumberOfCLK; n++){
        for (int j = 0; j < k; j++) {
          // put k "1" in the CLK
          clk[n]->setBit(table[j+k*n]);
        }
      }
    }
  }

  //second loop since a bit can be set twice, so clks have to be build completely first

  for (unsigned l = 0; l < lenBloom; l++){
    for (int n = 0; n < NumberOfCLK; n++){
      if (clk[n]->getBit(l)==1){
        majority[l]++;
      }
    }
    if (majority[l]>
          NumberOfCLK/2){
      clkres->setBit(l);
    }
  }
  clkres->copyToString(CLKout,lenBloom-1);
  string CLKouts = CLKout;
  delete[] CLKout;
  delete[] majority;
  for (int n = 0; n < NumberOfCLK; n++){
    delete clk[n];
  }
  delete clkres;
  return CLKouts;
} //end CreateCLKMajority


/**
 * picks from a vector of qgram k2 qgram acoording to their own weight
 *
 * @param qgrams
 * @param items_weight each qgram has its own weight
 * @param k2 number of new qgram to be drawn
 * @param password to draw qgrams
 * @return
 */
vector<string> pickWeightedRandomQgrams(vector<string> qgrams, vector<double> qgram_weight, unsigned k2, string password){
  vector<string> res;
  //set seed accourding to the user defined password
  seed_seq seed(password.begin(), password.end());
  //create a random engine
  default_random_engine engine(seed);
  //create interval
  vector<double> interval;
  for (double i = 1.0; i <= (double)qgrams.size(); i++  ){
    interval.push_back(i);
  }
  std::piecewise_constant_distribution<> dist(interval.begin(), interval.end(), qgram_weight.begin());

  // Demonstrate with k2 randomly generated numbers
  for (unsigned i = 0; i < k2; ++i)
  {
    // Generate random number using gen, distributed according to dist
    unsigned r = static_cast<unsigned>(dist(engine));
    // Save r
    res.push_back(qgrams[r-1]);

  }
  return res;
}
