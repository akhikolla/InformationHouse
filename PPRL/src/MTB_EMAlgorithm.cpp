/*
 * EmAlgorithm.cpp
 *
 *  Created on: 30.01.2017
 *      Author: schnell-42
 */

#include "MTB_EMAlgorithm.h"

/**
 * frequencies = varsA.size()*varsB.size()
 * m initiale m values for each method
 * u initiale u values for each method
 * p the probability that a possible match is in fact correct
 * epsilon
 */

MTB_EMAlgorithm::MTB_EMAlgorithm(vector<int> frequencies, vector<double> m,
		vector<double> u, double p, double epsilon) {
//	MAX_ITERATIONS = 500;
//	tol = 0.0001;
//	maxM = 0.9999;
//	minU =  0.0001;
	this->patternFrequencies = frequencies; //in some cases frequencies are the product of the number of records of data1 * number of records data2
	this->epsilon = epsilon; // if difference of m and u of previous iteration smaller epsilon, stop
	// algorithm
	this->patternCount = this->patternFrequencies.size();
	//variableCount
	this->variableCount = (int) round(log2(this->patternCount));
	for (int i = 0; i < this->patternCount; i++) {
		this->recordPairs += frequencies[i];
	}

	this->p.resize(MAX_ITERATIONS);
	this->p[0] = p;

	// create m and u and init
	////this->m = new double[MAX_ITERATIONS][this->variableCount];//to be removed
	//Grow rows of m by MAX_ITERATIONS
	this->m.resize(MAX_ITERATIONS);
	for (int i = 0; i < MAX_ITERATIONS; ++i) {
//		//Grow Columns by this->variableCount
		this->m[i].resize(this->variableCount);
	}
//
//	//this->u = new double[MAX_ITERATIONS][this->variableCount];//to be removed
//	//Grow rows of n by MAX_ITERATIONS
	this->u.resize(MAX_ITERATIONS);
	for (int i = 0; i < MAX_ITERATIONS; ++i) {
		//Grow Columns by this->variableCount
		this->u[i].resize(this->variableCount);
	}

	for (int i = 0; i < this->variableCount; i++) {
		this->m[0][i] = m[i];
		this->u[0][i] = u[i];
	}
//	// create gm, gu
//	//this->gm = new double[MAX_ITERATIONS][this->patternCount];//to be removed
	this->gm.resize(MAX_ITERATIONS);
	for (int i = 0; i < MAX_ITERATIONS; ++i) {
		//Grow Columns by this->patternCount
		this->gm[i].resize(this->patternCount);
	}
//	//this->gu = new double[MAX_ITERATIONS][this->patternCount];//to be removed
	this->gu.resize(MAX_ITERATIONS);
	for (int i = 0; i < MAX_ITERATIONS; ++i) {
		//Grow Columns by this->patternCount
		this->gu[i].resize(this->patternCount);
	}
	this->patternMatrix = getPatternMatrix(this->variableCount);

}
MTB_EMAlgorithm::~MTB_EMAlgorithm() {
}

bool MTB_EMAlgorithm::calculate() {
	this->iterationCount = 0;
	do {
		expectation(this->iterationCount + 1);
		maximization(this->iterationCount + 1);
		this->iterationCount++;
	} while ((delta() > this->epsilon)
			&& (this->iterationCount < MAX_ITERATIONS));
	return checkValidity();
}

string MTB_EMAlgorithm::getException() {
	return this->exception;
}

int MTB_EMAlgorithm::getIterations() {
	return this->iterationCount;
}

vector<double> MTB_EMAlgorithm::getMArray() {
	return this->m[this->iterationCount - 1];
}

bool MTB_EMAlgorithm::checkValidity() {
	if (delta() > this->epsilon) {
		this->exception = EX_CONVERGENCY;
		return false;
	}
	vector<double> mArray = getMArray();
	for (double element : mArray) {
//			if (Double->isNaN(element)) //TODO
//			{
//				this->exception = EX_M_NAN;
//				return false;
//			}
		if (element == 0) {
			this->exception = EX_M_0;
			return false;
		}
	}
	return true;
}

double MTB_EMAlgorithm::delta() {
	double delta = 0.0;
	for (unsigned i = 0; i < this->m[this->iterationCount - 1].size(); i++) {
		delta +=  abs(this->m[this->iterationCount][i]- this->m[this->iterationCount - 1][i])
				+ abs(this->u[this->iterationCount][i]- this->u[this->iterationCount - 1][i]);
	}
	return delta;
}

void MTB_EMAlgorithm::expectation(int i) {
	// --------- E-Step --------- //
	for (int patCount = 0; patCount < this->patternCount; patCount++) // PATTERNS
	{
		double mproduct = 1;
		double uproduct = 1;
		for (int varCount = 0; varCount < this->variableCount; varCount++) // VARIABLES
				{
			mproduct *= (pow(this->m[i - 1][varCount],//see Winkler "Using the EM Algorithm for Weight Computation in the Felligi Sunter
					this->patternMatrix[patCount][varCount])
					* pow((1 - (this->m[i - 1][varCount])),
							(1 - (this->patternMatrix[patCount][varCount]))));
			uproduct *= (pow(this->u[i - 1][varCount],
					this->patternMatrix[patCount][varCount])
					* pow((1 - (this->u[i - 1][varCount])),
							(1 - (this->patternMatrix[patCount][varCount]))));
		}
		this->gm[i][patCount] = (this->p[i - 1] * mproduct)	/ ((this->p[i - 1] * mproduct)+ ((1 - this->p[i - 1]) * uproduct));
		this->gu[i][patCount] = ((1 - this->p[i - 1]) * uproduct)/ ((this->p[i - 1] * mproduct)	+ ((1 - this->p[i - 1]) * uproduct));
	}
}

/**
 * Agreement Pattern
 *
 * Inits the patternMatrix assuming that the given pf has the following form (in case of pf = 3): freq[0] = kein agreement; freq[1] = agr var 1; freq[2] = agr var 2; freq[3] =
 * agr var 1+2; freq[4] = agr var 3; freq[5] = agr var 1+3; freq[6] = agr var 2+3; freq[7] = agr var 1+2+3;
 *
 * @param pf = number of merge (comparision) functions (Configurations )
 *
 * see Yancey 2004 "Improving EM Algorithm Estimates for
 * Record Linkage Parameters"
 * There are 2^pf possible agreement patterns in the comparison space Gamma, so the pattern Matrix contains 2^pf
 *
 *
 */

vector<vector<int>> MTB_EMAlgorithm::getPatternMatrix(int pf) {
	int pn = (int) round(pow(2, pf));
	vector<vector<int>> mat(pn, vector<int>(pf));
	for (int pat = 0; pat < pn; pat++) {
		// get binary representation of pattern number
		char a = pat;
		string tmp = bitset<32>(a).to_string();
		for (int b = 0; b < pf; b++) { // get the pf last binary numbers
			mat[pat][b] = stoi((tmp.substr(tmp.length() - b - 1, 1)));
		}
	}


	return mat;
}

void MTB_EMAlgorithm::maximization(int i) {
	for (int varCount = 0; varCount < this->variableCount; varCount++) // VARIALBES
			{
		double gmprod = 0;
		double guprod = 0;
		double gmprodnofreq = 0;
		double guprodnofreq = 0;
		for (int patCount = 0; patCount < this->patternCount; patCount++) // PATTERNS
				{
			gmprod += (this->gm[i][patCount]
					* this->patternMatrix[patCount][varCount]
					* this->patternFrequencies[patCount]);
			gmprodnofreq += (this->gm[i][patCount]
					* this->patternFrequencies[patCount]);
			guprod += (this->gu[i][patCount]
					* this->patternMatrix[patCount][varCount]
					* this->patternFrequencies[patCount]);
			guprodnofreq += (this->gu[i][patCount]
					* this->patternFrequencies[patCount]);
		}
		this->m[i][varCount] = (gmprod) / (gmprodnofreq);
		this->u[i][varCount] = (guprod) / (guprodnofreq);
		this->p[i] = gmprodnofreq / this->recordPairs;
	}
}

