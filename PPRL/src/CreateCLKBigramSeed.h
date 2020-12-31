//============================================================================
// Name        : CreateCLKBigramSeed.h
// Author      : Rukasz
// Version     :
// Copyright   : Your copyright notice
// Description : Ansi-style
//============================================================================


#ifndef CREATECLKBIGRAMSEED_H
#define CREATECLKBIGRAMSEED_H
#include <sstream>
#include "BloomHelper.h"
#include "CLK.h"
#include <random>
#include <iterator>
#include <type_traits>

using namespace std;

#define NONCE_SIZE 8


// void CreateCLKBigramSeedBit(CLK* clk, vector<string> vars, int k, vector<int> padding,
//                          vector<int> lenQgram, unsigned lenBloom, vector<string> password);

string CreateCLKBigramSeed(vector<string> vars, int k, vector<int> padding,
                         vector<int> lenQgram, unsigned lenBloom, vector<string> password);

//string CreateCLKSalsa20(vector<string> vars, int k, vector<int> padding,
//                         vector<int> lenQgram, unsigned lenBloom, string password, Rcpp::RawVector nonce);
string CreateBFBigramSeed(string vars, int k, int padding,
                        int lenQgram, unsigned lenBloom, string password);

string CreateMarkovCLKc(vector<string> vars, int k1, int k2, vector<int> padding,
                     vector<int> lenQgram, unsigned lenBloom, vector<string> password,
                     vector<vector<float>> markovTable, vector<string> rownames, vector<string> colnames, bool includeOriginalBigram, bool v);

string CreateEnsembleCLKc(vector<string> vars, int k, int NumberOfCLK, vector<int> padding,
                       vector<int> lenQgram, unsigned lenBloom, vector<string> password);
/**
 * helper function for markov
 * picks from a vector of qgram k2 qgram acoording to their own weight
 *
 * @param qgrams
 * @param items_weight each qgram has its own weight
 * @param k2 number of new qgram to be drawn
 * @param password to draw qgrams
 * @return
 */
vector<string> pickWeightedRandomQgrams(vector<string> qgrams, vector<double> qgram_weight, unsigned k2, string password);

#endif
