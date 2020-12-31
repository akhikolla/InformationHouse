#include <stdlib.h>
#include <stdio.h> // for console input/output
#include <math.h> // for exp and log functions
#include <vector>
using namespace std;

#include <Rcpp.h>
using namespace Rcpp;

#define LAMBDA_MIN_RATIO 0.000001
#define NUM_LAMBDA 100

static int my_min (int a, int b){
	/*
	* Returns the minimum of two integers
	*/ 
	if (a < b) return a;
	return b;
}
	
static int my_max(int a, int b){
	/*
	* Returns the maximum of two integers
	*/
	if (a > b) return a;
	return b;
}

static void getMeansInStratum (const vector<double>& X, int k, int p, int stratumPos, int stratumSize, double **xMean){
	/*
	*	Calculates the mean of each predictor in stratum k and stores the result in the appropriate part of xMean
	*	'stratumPos' is the number of elements in previous strata (so we know where we are in X and xMean)
	*	'stratumSize' is the current size of the stratum
	*/
	double currentSum;
	for (int j = 0; j < p; j++){ // loop over predictors
		currentSum = 0;
		for (int i = 0; i < stratumSize; i++){
			currentSum += X[stratumPos*p + i*p + j];
		}
		xMean[k][j] = currentSum/stratumSize; // store mean for current stratum
	}
}

static void getMeans (const vector<double>& X, int K, int p, const vector<int>& nVec, double **xMean){
	/*
	*	Computes the mean of each of the p predictors in each of the K strata and stores the result in xMean
	*/
	int sumN = 0;
	for (int k = 0; k < K; k++){
		getMeansInStratum (X, k, p, sumN, nVec[k], xMean);
		sumN += nVec[k];
	}
}

static void getCaseSumsInStratum (const vector<int>& y, const vector<double>& X, int k, int p, int stratumPos, int stratumSize, double **caseSums, double **xMean){
	/*
	*	Calculates, for every predictor (mean centered), the sum over the cases within stratum k
	*	Stores the result at the appropraite place in 'caseSums'
	*/
	double currentSum;
	
	for (int j = 0; j < p; j++){
		currentSum = 0;
		for (int i = 0; i < stratumSize; i++){
			if (y[stratumPos + i] == 1)	currentSum += (X[stratumPos*p + i*p + j] - xMean[k][j]); // case, so add contribution
		}
		caseSums[k][j] = currentSum;
	}
}

static void getCaseSums(const vector<int>& y, const vector<double>& X, int K, int p, const vector<int>& nVec, double **xMean, double **caseSums){
	/*
	*	Calculates the case sums for each predictor within each stratum
	*/
	int sumN = 0;
	for (int k = 0; k < K; k++){
		getCaseSumsInStratum(y, X, k, p, sumN, nVec[k], caseSums, xMean);
		sumN += nVec[k];
	}
}

static void allocateLinearPredictors(double ***linPred, int K, const vector<int>& nVec, double **maxLinPred, double **sumLinPredCases){
	/*
	*	Allocates memory space to the linear predictors (a ragged array with backbone of length K, each slot carrying just enough slots for a linear predictor for each observation in the stratum)
	*	All linear predictors are initially 0
	*	Also allocates space for the arrays storing the maximum linear predictor value and the sum of the linear predictors of the cases within each stratum
	*/
	*linPred = (double**) calloc (K, sizeof (double*));
	*maxLinPred = (double*) calloc (K, sizeof(double));
	*sumLinPredCases = (double*) calloc(K, sizeof(double));
	for (int k = 0; k < K; k++){
		(*linPred)[k] = (double*) calloc(nVec[k], sizeof(double)); // allocate space for the current stratum
	}
}

static void getLinearPredictorsInStratum (const vector<int>& y, const vector<double>& X, int k, double *beta, int nonZeroBeta, vector<int>& nonZeroBetaInd, int p, int stratumPos, int stratumSize, double *linPred, double *maxLinPred, double *sumLinPredCases){
	/*
	*	Calculates the linear precdictors within stratum k, keeping track of the maximum predictor value
	*	Places the linear predictors in 'linPred' and the maximum one in 'maxLinPred'
	*	Also computes the sum of the linear predictors of the cases and stores it in 'sumLinPredCases'
	*/
	int first = 1;
	double maxSoFar = 0, currentLinPred;
	
	*sumLinPredCases = 0;
	for (int i = 0; i < stratumSize; i++){
		currentLinPred = 0;
		for (int j = 0; j < p; j++){
			currentLinPred += beta[j]*X[stratumPos*p + i*p + j];
		}
		linPred[i] = currentLinPred; // save current linear predictor value
		if (y[stratumPos + i] == 1) *sumLinPredCases += currentLinPred; // if case, add to case sum
		if (first == 1 || currentLinPred > maxSoFar) maxSoFar = currentLinPred;
		first = 0;
	}
	*maxLinPred = maxSoFar;
}

static void getLinearPredictors(const vector<int>& y, const vector<double>& X, int K, double *beta, int nonZeroBeta, vector<int>& nonZeroBetaInd, int p, const vector<int>& nVec, double **linPred, double *maxLinPred, double *sumLinPredCases){
	/*
	*	Computes the linear predictors of all observations in all strata, using the provided beta vector
	*	Results are stored in 'linPred' - a ragged array with a backbone having a slot for each stratum
	*	Maximum linear predictors for every stratum are computed as we go along and stored in 'maxLinPred'
	*	Similarly, the sum of the linear predictors over all cases in each stratum is computed and stored in 'sumLinPredCases'
	*/
	int sumN = 0;
	for (int k = 0; k < K; k++){
		getLinearPredictorsInStratum(y, X, k, beta, nonZeroBeta, nonZeroBetaInd, p, sumN, nVec[k], linPred[k], &(maxLinPred[k]), &(sumLinPredCases[k]));
		sumN += nVec[k];
	}
}

static double getLikelihoodInStratum (vector<int>& whichVars, int numVars, int n, int m, double *linPred, double maxLinPred, double sumLinPredCases){
	/*
	*	Computes the likelihood and score for 'whichVars' in stratum k
	*/
	vector<double> normConst(m+1);// = (double*) calloc (m + 1, sizeof(double)); // will store the results of the recursive computations building up to the value of the normalising constant
	normConst[0] = 1;
	
	// recursive computations
	for (int i = 0; i < n; i++){ // loop over each observation in the stratum
		double u = exp(linPred[i]-maxLinPred);
		for (int j = my_min (m, i+1); j > my_max(0, m+i-n); j--){ // loop over column in the (n+1)x(m+1) matrix of recursive computations
			normConst[j] += normConst[j-1]*u; // update the normalising constant
		}
	}
	
	// add results to tallies
	return sumLinPredCases - m*maxLinPred - log(normConst[m]);
	//free(normConst);
}

static void getLikelihood (vector<int>& whichVars, int numVars, int K, const vector<int>& nVec, const vector<int>& mVec, double **linPred, double *maxLinPred, double *sumLinPredCases, double *loglik){
	/*
	*	Computes the likelihood and score (need to recursive calculations from the former for the latter) for all the variables the indices of which are listed in 'whichVars' (there are 'numVars' of them)
	* 	Arguments X, K, nVec, mVec, linPred, maxLinPred, sumLinPredCases, xMean and caseSums are all pre computed and ready for use in the recursive computation
	* 	The results will be stored in 'loglik' and 'score'
	*/
	int k;
	double tempSum = 0;
	//*loglik = 0;
	for (k = 0; k < K; k++){
		tempSum += getLikelihoodInStratum(whichVars, numVars, nVec[k], mVec[k], linPred[k], maxLinPred[k], sumLinPredCases[k]);
	}
	*loglik = tempSum;
}

static void getLikelihoodAndScoreInStratum (vector<int>& whichVars, int numVars, const vector<double>& X, int p, int stratumPos, int n, int m, double *linPred, double maxLinPred, double sumLinPredCases, const double *xMean, double *caseSums, double *loglik, double *score){
	/*
	*	Computes the likelihood and score for 'whichVars' in stratum k
	*/
	double *normConst = (double*) calloc (m + 1, sizeof(double)); // will store the results of the recursive computations building up to the value of the normalising constant
	double *normConstDeriv = (double*) calloc (numVars*(m+1), sizeof(double)); // will store the results of the recursive computations building up to the value of the derivative of the normalising constant
	normConst[0] = 1;
	
	// recursive computations
	for (int i = 0; i < n; i++){ // loop over each observation in the stratum
		double u = exp(linPred[i] - maxLinPred);
		for (int j = my_min (m, i+1); j > my_max(0, m+i-n); j--){ // loop over column in the (n+1)x(m+1) matrix of recursive computations
			for (int l = 0; l < numVars; l++){ // for each of the variables for which we require the score
				int current_l = whichVars[l];
				double *normConstDeriv_l = normConstDeriv + l*(m+1); // move to correct palce in array
				normConstDeriv_l[j] += u*(normConstDeriv_l[j-1] + normConst[j-1]*(X[stratumPos*p + i*p + current_l] - xMean[current_l]));
			}
			normConst[j] += normConst[j-1]*u; // update the normalising constant
		}
	}
	
	// add results to tallies
	if (stratumPos == 0) *loglik = 0;
	*loglik += (sumLinPredCases - m*maxLinPred - log(normConst[m]));
	for (int l = 0; l < numVars; l++){
		double *normConstDeriv_l = normConstDeriv + l*(m+1);
		int current_l = whichVars[l];
		if (stratumPos == 0) score[current_l] = 0;
		score[current_l] += (caseSums[current_l] - normConstDeriv_l[m]/normConst[m]);
	}
	
	free(normConst);
	free(normConstDeriv);
}

static void getLikelihoodAndScore (vector<int>& whichVars, int numVars, const vector<double>& X, int p, int K, const vector<int>& nVec, const vector<int>& mVec, double **linPred, double *maxLinPred, double *sumLinPredCases, double **xMean, double **caseSums, double *loglik, double *score){
	/*
	*	Computes the likelihood and score (need to recursive calculations from the former for the latter) for all the variables the indices of which are listed in 'whichVars' (there are 'numVars' of them)
	* 	Arguments X, K, nVec, mVec, linPred, maxLinPred, sumLinPredCases, xMean and caseSums are all pre computed and ready for use in the recursive computation
	* 	The results will be stored in 'loglik' and 'score'
	*/
	int sumN = 0;
	for (int k = 0; k < K; k++){
		getLikelihoodAndScoreInStratum(whichVars, numVars, X, p, sumN, nVec[k], mVec[k], linPred[k], maxLinPred[k], sumLinPredCases[k], xMean[k], caseSums[k], loglik, score);
		sumN += nVec[k]; // keep track of which row we are at in the X matrix
	}
}

static void getLikelihoodScoreAndHessianInStratum (vector<int>& whichVars, int numVars, const vector<double>& X, int p, int stratumPos, int n, int m, double *linPred, double maxLinPred, double sumLinPredCases, double *xMean, double *caseSums, double *loglik, double *score, double **hessian){
	/*
	*	Computes the likelihood and score for 'whichVars' in stratum k
	*/
	double *normConst = (double*) calloc (m + 1, sizeof(double)); // will store the results of the recursive computations building up to the value of the normalising constant
	double *normConstDeriv = (double*) calloc (numVars*(m+1), sizeof(double)); // will store the results of the recursive computations building up to the value of the derivative of the normalising constant
	double *normConstDeriv2 = (double*) calloc (numVars*numVars*(m+1), sizeof(double));
	normConst[0] = 1;
	
	// recursive computations
	for (int i = 0; i < n; i++){ // loop over each observation in the stratum
		double u = exp(linPred[i] - maxLinPred);
		for (int j = my_min (m, i+1); j > my_max(0, m+i-n); j--){ // loop over column in the (n+1)x(m+1) matrix of recursive computations
			for (int s = 0; s < numVars; s++){ // for each of the variables for which we require the score
				int current_s = whichVars[s]; // current variable number
				double *normConstDeriv_s = normConstDeriv + s*(m+1); // move to correct place in array
				for (int t = s; t < numVars; t++){
					int current_t = whichVars[t];
					double *normConstDeriv_t = normConstDeriv + t*(m+1);
					double *normConstDeriv2_st = normConstDeriv2 + (m+1)*(s + numVars*t);
					normConstDeriv2_st[j] += u*(normConstDeriv2_st[j-1] + 
										normConstDeriv_s[j-1]*(X[stratumPos*p + i*p + current_t] - xMean[current_t]) + 
										normConstDeriv_t[j-1]*(X[stratumPos*p + i*p + current_s] - xMean[current_s]) + 
										normConst[j-1]*(X[stratumPos*p + i*p + current_s] - xMean[current_s])*(X[stratumPos*p + i*p + current_t] - xMean[current_t]));
				}
				normConstDeriv_s[j] += u*(normConstDeriv_s[j-1] + normConst[j-1]*(X[stratumPos*p + i*p + current_s] - xMean[current_s]));
			}
			normConst[j] += normConst[j-1]*u; // update the normalising constant
		}
	}
	
	// add results to tallies
	if (stratumPos == 0) *loglik = 0;
	*loglik += (sumLinPredCases - m*maxLinPred - log(normConst[m]));
	for (int s = 0; s < numVars; s++){
		int current_s = whichVars[s];
		double *normConstDeriv_s = normConstDeriv + s*(m+1);
		if (stratumPos == 0) score[whichVars[s]] = 0;
		score[whichVars[s]] += (caseSums[whichVars[s]] - normConstDeriv_s[m]/normConst[m]);
		for (int t = s; t < numVars; t++){
			int current_t = whichVars[t];
			int minVar = my_min(current_s, current_t);
			double *normConstDeriv_t = normConstDeriv + t*(m+1);
			double *normConstDeriv2_st = normConstDeriv2 + (m+1)*(s + numVars*t);
			if (stratumPos == 0) hessian [minVar][abs(current_s - current_t)] = 0;
			hessian[minVar][abs(current_s - current_t)] += (normConstDeriv2_st[m]/normConst[m] - (normConstDeriv_t[m]/normConst[m])*(normConstDeriv_s[m]/normConst[m]));
		}
	}
	
	free(normConst);
	free(normConstDeriv);
	free (normConstDeriv2);
}

static void getLikelihoodScoreAndHessian (vector<int>& whichVars, int numVars, const vector<double>& X, int p, int K, const vector<int>& nVec, const vector<int>& mVec, double **linPred, double *maxLinPred, double *sumLinPredCases, double **xMean, double **caseSums, double *loglik, double *score, double **hessian){
	/*
	*	Computes the likelihood, score and hessian(need to recursive calculations from the former for the latter) for all the variables the indices of which are listed in 'whichVars' (there are 'numVars' of them)
	* 	Arguments X, K, nVec, mVec, linPred, maxLinPred, sumLinPredCases, xMean and caseSums are all pre computed and ready for use in the recursive computation
	* 	The results will be stored in 'loglik', 'score' and hessian
	*/
	int sumN = 0;
	for (int k = 0; k < K; k++){
		getLikelihoodScoreAndHessianInStratum(whichVars, numVars, X, p, sumN, nVec[k], mVec[k], linPred[k], maxLinPred[k], sumLinPredCases[k], xMean[k], caseSums[k], loglik, score, hessian);
		sumN += nVec[k]; // keep track of which row we are at in the X matrix
	}
}

static double getMaximumLambda (int p, double *score){
	/*
	*	The largest score element gives the smallest lambda for which all betas are 0 (so we assume that 'score' is computed at beta = 0)
	*	This is obtained from the coordinate descent update for the conditional likelihood model, plugging in beta = 0 on the right hand side of the update equation
	*/
	double maxL = 0;
	int first = 1;
	
	for (int i = 0; i < p; i++){
		if (first || fabs(score[i]) > maxL) maxL = fabs(score[i]);
		first = 0;
	}
	return maxL;
}

static double softThreshold(double x, double lambda){
	/*
	*	Returns a soft thresholding on 'x' at threshold 'lambda'
	*/
	
	if (x > lambda) return x - lambda;
	if (x < -lambda) return x + lambda;
	return 0;
}


static void allocateMemoryForLikelihoodComputations(int K, int p, const vector<int>& mVec, const vector<int>& nVec, double **score, double ***hessian, double ***xMean, double ***caseSums){
	/*
	*	Allocates memory to the various ragged arrays storing the important initial and intermediate computations required during the computation of all relevant
	*	likelihood quantities
	*/
	*score = (double*) calloc(p, sizeof(double));
	*hessian = (double**) calloc (p, sizeof(double*)); // backbone of hessian
	for (int i = 0; i < p; i++){
		(*hessian)[i] = (double*) calloc(p-i, sizeof(double)); // arrays containing hessian entries - store only upper triangle of hessian
	}
	//printf ("Score and hessian allocated\n");
	
	*xMean = (double**) calloc(K, sizeof(double*)); // a different mean vector for each stratum
	*caseSums = (double**) calloc(K, sizeof(double*)); // same for case sums
	for (int k = 0; k < K; k++){
		(*xMean)[k] = (double*) calloc(p, sizeof(double)); // one stratum mean for each predictor
		(*caseSums)[k] = (double*) calloc(p, sizeof(double));
	}
	//printf ("Means and case sums allocated\n");
}

static void deallocateMemoryForLikelihoodComputations(int K, int p, const vector<int>& mVec, double **score, double ***hessian, double ***xMean, double ***caseSums, double ***linPred, double **maxLinPred, double **sumLinPredCases){
	/*
	*	Deallocates all the memory previously allocated for use in computing likelihood tallies
	*/
	free (*score);
	for (int i = 0; i < p; i++){
		free ((*hessian)[i]);
	}
	free (*hessian);
	//printf ("Deallocated score and hessian\n");
	
	for (int k = 0; k < K; k++){
		free((*xMean)[k]);
		free((*caseSums)[k]);
		free((*linPred)[k]);
		//printf ("Startum %d: freed means and case sums\n", k);
	}
	free(*xMean);
	free(*caseSums);
	free(*linPred);
	free(*maxLinPred);
	free(*sumLinPredCases);
}

static double getPenaltyTerm(const double *beta, int nonZeroBeta, const vector<int>& nonZeroBetaInd, double lambda, double alpha){
	/*
	*	Computes the l1 penalty term associated with the current beta and lambda values
	*/
	double pen = 0;
	for (int i = 0; i < nonZeroBeta; i++){
		pen += (alpha*fabs(beta[nonZeroBetaInd[i]]) + (1-alpha)*pow(beta[nonZeroBetaInd[i]], 2)/2);
	}
	pen *= lambda;
	return pen;
}

static void coordinateDescent (vector<int>& whichVars, int numVars, double lambda, double *topOfLoopBeta, double *currentBeta, double *score, double **hessian, double initPenalty, double eps, double alpha){
	/*
	*	Runs one epoch of coordinate descent, using the 'score', 'hessian', 'beta' and 'lambda' in successive cyclic univariate updates to construct a new beta
	*	Visits only those variables in 'whichVars'
	*   By the time we have converged (convergence is measured in changes in the quadratic approximation to the loglikelihood), 'currentBeta' will contain the new betas
	*/
	int first = 1;
	double oldQuadApprox, oldPenalty, offDiagonalHessian, newBeta;
	double quadApprox = 0; // quadratic approximation relative to current likelihood
	double penalty = initPenalty;
	while (first || fabs((-quadApprox + oldQuadApprox + penalty - oldPenalty)/(-quadApprox + penalty)) > eps){
		oldQuadApprox = quadApprox; 
		oldPenalty = penalty; // save the current values
		penalty = 0;
		quadApprox = 0;
		for (int j = 0; j < numVars; j++){ // visit each of the variables in 'whichVars'
			int current_j = whichVars[j];
			offDiagonalHessian = 0;
			for (int i = 0; i < numVars; i++){ // off diagonal hessian contribution
				int current_i = whichVars[i];
				if (current_i != current_j){
					int minVar = my_min(current_i, current_j);
					offDiagonalHessian += (topOfLoopBeta[current_i] - currentBeta[current_i])*hessian[minVar][abs(current_i - current_j)];
				}
			}
			newBeta = softThreshold(topOfLoopBeta[current_j]*hessian[current_j][0] + score[current_j] + offDiagonalHessian, lambda*alpha)/(hessian[current_j][0] + lambda*(1-alpha)); // beta update
			currentBeta[current_j] = newBeta;
				
			penalty += (alpha*fabs(currentBeta[current_j]) + (1-alpha)*pow(currentBeta[current_j], 2)/2);
			quadApprox += score[current_j]*(currentBeta[current_j] - topOfLoopBeta[current_j]) - offDiagonalHessian - 0.5*hessian[current_j][0]*pow(currentBeta[current_j] - topOfLoopBeta[current_j], 2);
		}
		penalty *= lambda;
		first = 0;
	}
}

static int getStrongSet (const double *score, double currentLambda, double previousLambda, int p, vector<int>& strongSet, vector<int>& notStrongSet, double alpha){
	/*
	*	Applies the sequential strong rule to determine which variables will feature in the next round of cyclic coordinate descent steps
	*	Returns the number of elements in the strong set and populates the first few elements of 'strongSet' with the appropriate indices.
	*	The indices of the all variables that fail the sequential strong rule are placed in 'notStrongSet'
	*/
	int count = 0, notCount = 0;
	//beta(l, p + 7) = 2*currentLambda - previousLambda;
	for (int i = 0; i < p; i++){
		//beta(l, p+8+i) = score[i];
		if (fabs(score[i]) < alpha*(2*currentLambda - previousLambda)){ // fails sequential strong rule
			notStrongSet[notCount] = i;
			notCount++;
		}
		else{
			strongSet[count] = i;
			count++;
		}
	}
	return count;
}

NumericMatrix fit_cloglik (const vector<int>& y, const vector<double>& X, int K, int p, const vector<int> &mVec, const vector<int> &nVec, int numLambda, int switchIter, double minLambda, double alpha){
	/*
	*	Function for fitting the conditional loglikelihood using coordinate descent
	*		input:
	*			- y = Pointer to vector marking cases (1) and controls (0)
	*			- X = Pointer to predictor matrix 
	*			- K = Number of strata
	*			- p = Number of predictors (columns of X)
	*			- mVec = Vector of length K containing the number of cases per stratum
	*			- nVec = Vector of length K containing the number of observations per stratum
	*/
	
	int firstMiddle, count, numViolations, nonZeroBeta;
	double loglik, oldloglik, penalty, oldpenalty, nullDev;
	double *score, **hessian, **xMean, **caseSums, **linPred, *maxLinPred, *sumLinPredCases;
	double eps = 0.0001;
	
	// create vectors for strong set bookeeping
	vector<int> strongSet(p);
	vector<int> notStrongSet(p);
	int numStrongSet = 0;
	for (int i = 0; i < p; i++) notStrongSet[i] = i;
	
	int totalObs = 0;
	for (int k = 0; k < K; k++) totalObs += nVec[k];
	
	// declare matrix that will house the betas
	NumericMatrix beta = NumericMatrix (numLambda + 1, p + 7);
	
	// initial computations
	allocateMemoryForLikelihoodComputations(K, p, mVec, nVec, &score, &hessian, &xMean, &caseSums);
	getMeans(X, K, p, nVec, xMean);
	getCaseSums(y, X, K, p, nVec, xMean, caseSums);
	allocateLinearPredictors (&linPred, K, nVec, &maxLinPred, &sumLinPredCases);
	//getLinearPredictors(y, X, K, beta[0], 0, strongSet, p, nVec, linPred, maxLinPred, sumLinPredCases);
	getLikelihoodAndScore (notStrongSet, p-numStrongSet, X, p, K, nVec, mVec, linPred, maxLinPred, sumLinPredCases, xMean, caseSums, &loglik, score);
	penalty = 0;
	nullDev = -2*loglik; // maximised likelihood of saturated model is 0
	
	// lambda bookkeeping
	vector<double> lambdaVec(numLambda + 1);
	lambdaVec[0] = getMaximumLambda(p, score)/alpha;
	//if (switchIter != 100) minLambdaRatio = LAMBDA_MIN_RATIO; // can go to 0 only if we have linear grid
	
	// make space for the beta working variables
	double *topOfLoopBeta = (double*) calloc (p, sizeof(double)); // 'topOfLoopBeta' always points to an array containing the initial betas in the current quadratic approximation
	double *currentBeta = (double *) calloc(p, sizeof(double));; // current beta contains the betas currently undergoing cyclic coordinate descent
	beta(0, p) = lambdaVec[0];
	
	for (int l = 1; l <= numLambda; l++){ // loop over lambdas
		for (int i = 0; i < p; i++){
			topOfLoopBeta[i] = beta(l-1, i); // initialise top of loop betas
		}
		if (l <= switchIter) lambdaVec[l] = lambdaVec[l-1] - (lambdaVec[0] - lambdaVec[0]*minLambda)/numLambda;
		else lambdaVec[l] = exp(log(lambdaVec[l-1]) + log(minLambda)/numLambda);
		//printf ("Current lambda; maximum lambda: %f\t%f\n", lambdaVec[l], lambdaVec[0]);
		
		// find the strong set
		numStrongSet = getStrongSet(score, lambdaVec[l], lambdaVec[l-1], p, strongSet, notStrongSet, alpha); // get the strong set
		
		numViolations = 1; // just make sure we break into this loop
		while (numViolations > 0){
			firstMiddle = 1; // about to have the first middle loop iteration for this lambda value
			count = 0;
			//countUp = 0;
			while (firstMiddle || fabs((-loglik+oldloglik+penalty-oldpenalty)/(-loglik + penalty)) > eps){ // quadratic approximation to the loglikelihood - middle loop
				count++;
				//printf ("\tMiddle iteration %d\n", count);
				oldloglik = loglik;
				oldpenalty = penalty;
				getLikelihoodScoreAndHessian (strongSet, numStrongSet, X, p, K, nVec, mVec, linPred, maxLinPred, sumLinPredCases, xMean, caseSums, &loglik, score, hessian);
			
				for (int i = 0; i < p; i++) currentBeta[i] = topOfLoopBeta[i]; // copy top of loop betas into current beta
			
				coordinateDescent(strongSet, numStrongSet, lambdaVec[l], topOfLoopBeta, currentBeta, score, hessian, penalty, eps, alpha);
			
				getLinearPredictors (y, X, K, currentBeta, numStrongSet, strongSet, p, nVec, linPred, maxLinPred, sumLinPredCases);
				getLikelihood (strongSet, numStrongSet, K, nVec, mVec, linPred, maxLinPred, sumLinPredCases, &loglik);
				penalty = getPenaltyTerm(currentBeta, numStrongSet, strongSet, lambdaVec[l], alpha);
				//if (-loglik+oldloglik+penalty-oldpenalty > 0) countUp++; // increase by 1 if we do not decrease the objective
				
				// make sure top of loop beta now points to current beta (avoiding having to copy current beta into top of loop beta
				double *temp = topOfLoopBeta;
				topOfLoopBeta = currentBeta;
				currentBeta = temp;
				
				firstMiddle = 0;
				if (count == 100) break;
			}
		
			// check KKT conditions
			if (count != 100){
				getLikelihoodAndScore (strongSet, numStrongSet, X, p, K, nVec, mVec, linPred, maxLinPred, sumLinPredCases, xMean, caseSums, &loglik, score);
				getLikelihoodAndScore (notStrongSet, p-numStrongSet, X, p, K, nVec, mVec, linPred, maxLinPred, sumLinPredCases, xMean, caseSums, &loglik, score); // get the score of all the other variables
				numViolations = 0;
				for (int i = 0; i < p-numStrongSet; i++){ // loop over predictors that are not in the strong set
					if (fabs(score[notStrongSet[i]] - lambdaVec[l]*(1-alpha)*currentBeta[notStrongSet[i]]) > lambdaVec[l]) {
						beta(l, p+6) = 1; // did we have any violations here?
						numViolations++;
					
						strongSet[numStrongSet] = notStrongSet[i];
						numStrongSet++; // put into strong set
					
						// remove the offending index from the nonstrong set
						for (int j = i; j < p-numStrongSet-1;j++){
							notStrongSet[j] = notStrongSet[j+1];
						}
						//printf ("Violation!");
					}
				}
			}
		}
		// copy the final betas for this lambda into our NumericMatrix
		nonZeroBeta = 0;
		for (int i = 0; i < p; i++){
			beta(l, i) = topOfLoopBeta[i];
			if (fabs(beta(l, i)) > 0) nonZeroBeta++;
			//printf ("%f ", beta(l, i));
		}
		
		beta (l, p) = lambdaVec[l]; // save lambda value
		beta (l, p+1) = nonZeroBeta; // save number of non-zero betas
		beta (l, p+2) = numStrongSet; // save number of strong set betas
		beta (l, p+3) = count; // save the iteration count (middle
		double dev = -2*loglik;
		beta (l, p+4) = 1 - dev/nullDev; //save the proportion of deviance explained
		if (nullDev-dev >= 0.99*nullDev) break; // for wide datasets - see Simon et al (2011)
		//printf ("Active betas: %d\n", numStrongSet);
	}
	//printf ("\nOut of for loop\n");
	deallocateMemoryForLikelihoodComputations(K, p, mVec, &score, &hessian, &xMean, &caseSums, &linPred, &maxLinPred, &sumLinPredCases);
	free (currentBeta);
	free (topOfLoopBeta);
	return beta;
}

static void copyData (NumericVector yvec, NumericMatrix xmat, vector<int>& y, vector<double>& X){
	/*
	*	Copies elements of 'yvec' into 'y' and elements of 'xmat' into 'X'
	*/
	int rows = xmat.nrow(), cols = xmat.ncol();
	for (int i = 0; i < rows; i++){
		for (int j = 0; j < cols; j++){
			X[i*cols + j] = xmat(i,j);
		}
		y[i] = (int)yvec[i];
	}
}

// [[Rcpp::export]]
NumericMatrix clogitl1_c (NumericVector n, NumericVector m, NumericMatrix Xmat, NumericVector yvec, int switchIter, int numLambda, double minLambda, double alpha){
	// find the number of strata and initialise mVec and nVec, copying the relevant values into them
	int K = n.size();
	vector<int> mVec(K);// = (int*) calloc (K, sizeof(int));
	vector<int> nVec(K);// = (int*) calloc (K, sizeof(int));
	int nSum = 0, mSum = 0;
	for (int k = 0; k < K; k++){
		nVec[k] = (int)n[k];
		mVec[k] = (int)m[k];
		nSum += nVec[k];
		mSum += mVec[k];
	}
	int p = Xmat.ncol();
	
	// allocate space for y and X, and copy into them
	vector<int> y(nSum); // = (int*) calloc (nSum, sizeof(int));
	vector<double> X(nSum*p);// = (double*) calloc (nSum*p, sizeof(double));
	copyData(yvec, Xmat, y, X);
	//printf ("Data copied\n");
	
	// fit the model
	NumericMatrix beta = fit_cloglik(y, X, K, p, mVec, nVec, numLambda, switchIter, minLambda, alpha);

	//free (mVec);
	//free (nVec);
	//printf ("Freed m and n vecs\n");
	//free (y);
	//free (X); // deallocate space
	//printf ("Freed y and X\n");
	return beta;
}
