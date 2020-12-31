/*
 * Ngram.cpp
 *
 *  Created on: Feb 20, 2017
 *      Author: rukasz
 */
#include "MTB_Ngram.h"

vector<string> CreateQgrams(string var, unsigned lenQgram) {
  //Only uni-, bi-, tri- and fourgram are allowed (for now)
  //if length of variable is smaller than length qgram then
  if (var.size()<=lenQgram){
    vector<string> qgrams = {var};
    return qgrams;
  }
  int vectorSize;
  if (lenQgram > 4||lenQgram < 1) {
       Rcpp::Rcerr
    << "Only q-grams of length between 1 and 4 are allowed. Q-grams of length more than 4 or less than 1 are not supported. Length of the qgram will be set to 2."
    << endl;

    vectorSize = var.size()-lenQgram+1;
  }
  else if(lenQgram == 1){

    vectorSize = var.size();

  }
  else{ // if lenQgram == 2
    vectorSize = var.size()-lenQgram+1;
  }
  string rest = var;
  vector<string> qgrams(vectorSize);
  int i = 0;
  while(rest.length() >=lenQgram){
    qgrams[i] = rest.substr(0, lenQgram);
    rest = rest.substr(1,rest.length());
    i++;
  }
  return qgrams;
} //end CreateQgrams

double NgramDistance(vector<string> qgrams1, vector<string> qgrams2){
	double res=0.0;
	// get ngrams
			if ((qgrams1.size() == 0) || (qgrams2.size() == 0))
			{
				return res;
			}
			// sort both
			sort(qgrams1.begin(), qgrams1.end());
			sort(qgrams2.begin(), qgrams2.end());
			// walk through and count identical
			int identicalNgrams = 0;
			for (unsigned pos1 = 0, pos2 = 0; (pos1 < qgrams1.size()) && (pos2 < qgrams2.size());)
			{
				int c = qgrams1[pos1].compare(qgrams2[pos2]);
				if (c < 0)
				{
					pos1++;
				}
				else if (c > 0)
				{
					pos2++;
				}
				else
				{
					identicalNgrams++;
					pos1++;
					pos2++;
				}
			}

			return (double)(2 * identicalNgrams) / (double)(qgrams1.size() + qgrams2.size());
		}

string MTB_NgramDistanceAlgorithm::getName(){
	return "NgramDistance";
}


double MTB_NgramDistanceAlgorithm::getRelativeValue(string var1, string var2)
	{
		return NgramDistance(CreateQgrams(var1,this->l1), CreateQgrams(var2, this->l2));
	}

double MTB_NgramDistanceAlgorithm::getAbsoluteValue(string var1, string var2)
	{
		return getRelativeValue(var1,var2);
	}
