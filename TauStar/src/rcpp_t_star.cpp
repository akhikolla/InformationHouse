/***
 * Copyright (C) 2015 Luca Weihs
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/***
 * A collection of functions that implement the main functionality of the
 * algorithms described in
 *
 * Weihs, Luca, Mathias Drton, and Dennis Leung. "Efficient Computation of the
 * Bergsma-Dassios Sign Covariance." arXiv preprint arXiv:1504.00964 (2015).
 */

#include<RcppArmadilloExtensions/sample.h>
#include"red_black_tree.h"
#include<stdio.h>
#include<ctype.h>
#include<cmath>
#include<algorithm>

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::interfaces(r, cpp)]]

// The following functions are used as input to initalize a red black tree whose
// nodes contain doubles. Knowing why these functions exist is likely
// unimportant to you, if otherwise then please read the red_black_tree.cpp file
// to see where they are used.
void DoubDest(void* a) { free((double*)a); }
int DoubComp(const void* a, const void* b) {
	if (*(double*)a > *(double*)b) {
    return(1);
	}
	if (*(double*)a < *(double*)b) {
    return(-1);
	}
	return(0);
}
void DoubPrint(const void* a) { Rprintf("%f", *(double*)a); }
void InfoPrint(void* a) { }
void InfoDest(void *a) { }

/***
 * A simple bubble sort implementation, this is only ever used on a vec of size
 * 4 in the following functions and hence its inefficiency is unimportant.
 * @param vec the array to sort
 * @param n the number of elements in the vec array
 */
void bubbleSort(long double * vec, int n) {
	int i,j,k;
	long double tmp;
	for(i = n - 1; i > 0; i--) {
		for(j = 0; j < i; j++) {
			if(vec[j + 1] < vec[j]) {
				tmp = vec[j + 1];
				vec[j + 1] = vec[j];
				vec[j] = tmp;
			}
		}
	}
}

/***
 * For an input vector x = (x_1,...,x_n) returns a new vector r of length n
 * with r[i] equal to the rank of x_i in x. So if x == (4,2,2,3,2,3) then
 * r == (3, 1, 1, 2, 1, 2).
 * @param x an arma::vec to get the ranks of
 * @return the vector of ranks
 */
arma::uvec vecToRanks(const arma::vec& x) {
  if (x.n_elem == 0) {
    return(arma::zeros<arma::uvec>(0));
  }
  arma::uvec xSortedInds = arma::sort_index(x);
  arma::uvec ranks = arma::zeros<arma::uvec>(x.n_elem);
  arma::uvec repeats = arma::zeros<arma::uvec>(x.n_elem);
  int k = 1;
  double last = x[xSortedInds[0]];
  for (int i = 0; i < x.n_elem; i++) {
    double current = x[xSortedInds[i]];
    if (current != last) {
      k++;
      last = current;
    }
    repeats[i] = k;
  }
  for (int i = 0; i < xSortedInds.n_elem; i++) {
    ranks[xSortedInds[i]] = repeats[i];
  }
  return ranks;
}

/***
 * For two input vectors of ranks xRanks, yRanks (i.e. taking values in the
 * positive integers and including all integers between 1 and their maximum
 * value) returns the matrix A where A[i,j] is the number times and element of
 * xRanks, and the corresponding element of yRanks, are less than or equal to i
 * and j respectively. Hence (xRanks[s],yRanks[s]) contributes to A[i,j] if and
 * only if x[s] <= i and y[s] <= i.
 * @param xRanks a vector of ranks.
 * @param yRanks a vector of ranks (must be the same length as x).
 * @return the "less than or equal" to matrix.
 */
arma::umat ranksToLeqMat(const arma::uvec& xRanks, const arma::uvec& yRanks) {
  int xMax = xRanks.max();
  int yMax = yRanks.max();
  arma::umat leqMat = arma::zeros<arma::umat>(xMax + 1, yMax + 1);
  for (int i = 0; i < xRanks.n_elem; i++) {
    leqMat(xRanks[i], yRanks[i]) += 1;
  }
  for (int i = 1; i <= xMax; i++) {
    for (int j = 1; j <= yMax; j++) {
      leqMat(i,j) =
        leqMat(i,j - 1) + leqMat(i - 1,j) - leqMat(i - 1, j - 1) + leqMat(i,j);
    }
  }
  return leqMat;
}

/***
 * Computes a version of the matrix defined at the top of page 5 of Heller and
 * Heller (2016) (arxiv:1605.08732). In particular, we let the rank by always
 * <= rather than < in the x case. This makes the notation nice in
 * TStarHellerAndHellerRCPP.
 * @param leqMat a leqMat as created by ranksToLeqMat.
 * @return the unique count matrix.
 */
arma::umat leqMatToUniqueCountMat(const arma::umat& leqMat) {
  arma::umat uCountMat = arma::zeros<arma::umat>(leqMat.n_rows, leqMat.n_cols);
  for (int i = 1; i < leqMat.n_rows; i++) {
    for (int j = 1; j < leqMat.n_cols; j++) {
      int numEqYAndLessOrEqX = leqMat(i, j) - leqMat(i, j - 1);
      uCountMat(i,j) = uCountMat(i, j - 1) +
        (numEqYAndLessOrEqX * (numEqYAndLessOrEqX - 1)) / 2;
    }
  }
  return uCountMat;
}

/***
 * Given an arma::uvec x and a set of indices inds, creates a new vector whose
 * entries are those of x indexed by inds. So if
 * x = (1,2,3)
 * inds = (0,0,2)
 * then we would return the vector (1,1,3).
 * @param x vector be be indexed into
 * @param inds the indices
 * @return a new vector.
 */
arma::uvec indexUvec(const arma::uvec& x, const arma::uvec& inds) {
  arma::uvec newVec = arma::zeros<arma::uvec>(inds.n_elem);
  for (int i = 0; i < newVec.n_elem; i++) {
    newVec[i] = x[inds[i]];
  }
  return newVec;
}

/***
 * Function implementing the O(n^2) computation of the t* U-statistic
 * described in
 *
 * Heller, Yair and Heller, Ruth. "Computing the Bergsma Dassios sign-
 * covariance." arXiv preprint arXiv:1605.08732 (2016).
 *
 * @param x a arma::vec of values
 * @param y a arma::vec of values of the same length as x
 * @return the U-statistic computed upon the two input vectors.
 */
// [[Rcpp::export]]
double TStarHellerAndHellerRCPP(const arma::vec& x, const arma::vec& y) {
  arma::uvec xRanks = vecToRanks(x);
  arma::uvec yRanks = vecToRanks(y);
  arma::umat leqMat = ranksToLeqMat(xRanks, yRanks);
  arma::uvec xOrder = arma::sort_index(xRanks);
  xRanks = indexUvec(xRanks, xOrder);
  yRanks = indexUvec(yRanks, xOrder);

  arma::umat uCountMat = leqMatToUniqueCountMat(leqMat);

  double numCon = 0;
  double numDis = 0;
  for (int i = 0; i < xRanks.n_elem - 1; i++) {
    for (int j = i + 1; j < xRanks.n_elem; j++) {
      int xRankMin = xRanks[i];
      int yRankMin = std::min(yRanks[i], yRanks[j]);
      int yRankMax = std::max(yRanks[i], yRanks[j]);

      int bot = leqMat(xRankMin - 1, yRankMin - 1);
      int mid = (yRankMin == yRankMax) ? 0 :
        leqMat(xRankMin - 1, yRankMax - 1) - leqMat(xRankMin - 1, yRankMin);
      int top = leqMat(xRankMin - 1, leqMat.n_cols - 1) -
                 leqMat(xRankMin - 1, yRankMax);
      int eqMin = leqMat(xRankMin - 1, yRankMin) -
        leqMat(xRankMin - 1, yRankMin - 1);
      int eqMax = leqMat(xRankMin - 1, yRankMax) -
        leqMat(xRankMin - 1, yRankMax - 1);

      numCon += top * (top - 1) / 2.0 + bot * (bot - 1) / 2.0;

      if (yRankMin != yRankMax) {
        numDis += top * (mid + eqMin + bot) + bot * (mid + eqMax) +
          eqMin * (mid + eqMax) +
          eqMax * mid + mid * (mid - 1) / 2.0;
        numDis -= uCountMat(xRankMin - 1, yRankMax - 1) -
          uCountMat(xRankMin - 1, yRankMin);
      }
    }
  }

  int n = xRanks.n_elem;
  double c = 16 * numCon - 8 * numDis;
  double d = (c < 0) ? -1 : 1;
  return d * expl(logl(d * c) - (logl(n) + logl(n - 1) +
                               logl(n - 2) + logl(n - 3)));
}

/***
 * Function implementing the O(n^2) computation of the t* V-statistic
 * using a modified version of the algorithm used in TStarHellerAndHellerRCPP.
 *
 * Heller, Yair and Heller, Ruth. "Computing the Bergsma Dassios sign-
 * covariance." arXiv preprint arXiv:1605.08732 (2016).
 *
 * @param x a arma::vec of values
 * @param y a arma::vec of values of the same length as x
 * @return the V-statistic computed upon the two input vectors.
 */
// [[Rcpp::export]]
double VTStarHellerAndHellerRCPP(const arma::vec& x, const arma::vec& y) {
  arma::uvec xRanks = vecToRanks(x);
  arma::uvec yRanks = vecToRanks(y);
  arma::umat leqMat = ranksToLeqMat(xRanks, yRanks);
  arma::uvec xOrder = arma::sort_index(xRanks);
  xRanks = indexUvec(xRanks, xOrder);
  yRanks = indexUvec(yRanks, xOrder);

  arma::umat uCountMat = leqMatToUniqueCountMat(leqMat);

  double numCon = 0;
  double numDis = 0;
  for (int i = 0; i < xRanks.n_elem; i++) {
    int top = leqMat(xRanks[i] - 1, leqMat.n_cols - 1) -
      leqMat(xRanks[i] - 1, yRanks[i]);
    int bot = leqMat(xRanks[i] - 1, yRanks[i] - 1);
    numCon += 0.5 * (top * (top - 1) / 2.0 + bot * (bot - 1) / 2.0) +
      0.25 * (top + bot);

    for (int j = i + 1; j < xRanks.n_elem; j++) {
      int xRankMin = xRanks[i];
      int yRankMin = std::min(yRanks[i], yRanks[j]);
      int yRankMax = std::max(yRanks[i], yRanks[j]);

      int bot = leqMat(xRankMin - 1, yRankMin - 1);
      int mid = (yRankMin == yRankMax) ? 0 :
        leqMat(xRankMin - 1, yRankMax - 1) - leqMat(xRankMin - 1, yRankMin);
      int top = leqMat(xRankMin - 1, leqMat.n_cols - 1) -
        leqMat(xRankMin - 1, yRankMax);
      int eqMin = leqMat(xRankMin - 1, yRankMin) -
        leqMat(xRankMin - 1, yRankMin - 1);
      int eqMax = leqMat(xRankMin - 1, yRankMax) -
        leqMat(xRankMin - 1, yRankMax - 1);

      numCon += top * (top - 1) / 2.0 + bot * (bot - 1) / 2.0 + 0.5 * (top + bot);

      if (yRankMin != yRankMax) {
        numDis += top * (mid + eqMin + bot) + bot * (mid + eqMax) +
          eqMin * (mid + eqMax) +
          eqMax * mid + mid * (mid - 1) / 2.0;
        numDis -= uCountMat(xRankMin - 1, yRankMax - 1) -
          uCountMat(xRankMin - 1, yRankMin);
      }
    }
  }

  int n = xRanks.n_elem;
  double c = 16 * numCon - 8 * numDis;
  double d = (c < 0) ? -1 : 1;
  return d * expl(logl(d * c) - 4 * logl(n));
}

/***
 * The main function implementing the O(n^2*log(n)) computation of the t*
 * U-statistic described in
 *
 * Weihs, Luca, Mathias Drton, and Dennis Leung. "Efficient Computation of the
 * Bergsma-Dassios Sign Covariance." arXiv preprint arXiv:1504.00964 (2015).
 *
 * NOTE: Assumes that the input xNumeric vector has been sorted from least to
 * greatest and that the yNumeric vector values have been rearranged to match.
 *
 * @param xNumeric a NumericVector of values
 * @param yNumeric a NumericVector of values of the same length as xNumeric
 * @return the U-statistic computed upon the two input vectors.
 */
// [[Rcpp::export]]
Rcpp::NumericVector TStarWeihsEtAlRCPP(NumericVector xNumeric,
                                      NumericVector yNumeric) {
	// Copy input data to something understood by the rbtree functions
  int n = xNumeric.size();
  double *x = (double*) malloc(sizeof(double) * n);
  double *y = (double*) malloc(sizeof(double) * n);
	for(int i = 0; i < n; i++) {
		x[i] = xNumeric[i];
		y[i] = yNumeric[i];
	}

	// Set up variables / data structures that will be needed during forward
  // loop over all data
	rb_red_blk_tree* tree =
    RBTreeCreate(DoubComp, DoubDest, InfoDest, DoubPrint, InfoPrint);
  int *savedYsInds = (int*) malloc(sizeof(int) * n); // Used to store y values
                                                     // corresponding to
                                                     // duplicate x values
  int numSavedYs = 0; // Number of y values actually saved in above. Note that
                      // the algorithm in the paper uses a list instead of the
                      // above two variables but the running time is uneffected
                      // by this change.
  long double concord = 0; // running count of # concordant quadruples
  long double discord = 0; // running count of # discordant quadruples
  int totalYInTree = 0;
	double lastX = 0; // the last x value seen, used to determine duplicates

  // Begin forward loop
	for (int i = 0; i < n - 1; i++) {
		if (lastX == x[i] && i != 0) {
      // If new x equal to last x then save the most recent y value
			savedYsInds[numSavedYs] = i;
			numSavedYs += 1;
		} else {
      // Otherwise insert all previously saved y values into the tree, delete
      // the old saved values, and save the most recent y value.
			for (int k = 0; k < numSavedYs; k++) {
				totalYInTree += 1;
				RBTreeInsert(tree, &(y[savedYsInds[k]]), 0);
			}
			lastX = x[i];
			numSavedYs = 1;
			savedYsInds[0] = i;
		}

    // Loop over the x,y values with indices greater than i and use the
    // equations from the paper to record the number of concordant/discordant
    // quadruples.
		for (int j = i + 1; j < n; j++) {
			double minY = std::min(y[i], y[j]);
			double maxY = std::max(y[i], y[j]);

			int numLessMax = RBNumLessThan(tree, &maxY);
			int numGreaterMin = RBNumGreaterThan(tree, &minY);
			int top = RBNumGreaterThan(tree, &maxY);
			int bottom = RBNumLessThan(tree, &minY);

      int middle = 0;
			if(minY != maxY) {
				middle = numGreaterMin + numLessMax - totalYInTree;
			}

			concord += ((bottom * (bottom - 1)) / 2) + ((top * (top - 1)) / 2);

			if(minY != maxY) {
				discord += ((middle * (middle - 1)) / 2) + top * middle + top * bottom +
          middle * bottom;
				int numYMin = totalYInTree - numGreaterMin - bottom;
				int numYMax = totalYInTree - top - numLessMax;
				discord += numYMin * numGreaterMin + numYMax * (numLessMax - numYMin);
			}
		}
	}

  // Now run everything in reverse to get rid of quadruples over-counted
  // as discordant
  rb_red_blk_tree* revTree = RBTreeCreate(DoubComp, DoubDest, InfoDest, DoubPrint, InfoPrint);
	numSavedYs = 0;
	lastX = 0;
	for(int i = n - 1; i > 0; i--) {
		if(lastX == x[i] && i != n - 1) {
			savedYsInds[numSavedYs] = i;
			numSavedYs += 1;
		} else {
			for(int k = 0; k < numSavedYs; k++) {
				RBTreeInsert(revTree, &(y[savedYsInds[k]]), 0);
			}
			lastX = x[i];
			numSavedYs = 1;
			savedYsInds[0] = i;
		}
		for(int j = i - 1; j >= 0; j--) {
			double minY = std::min(y[i], y[j]);
			double maxY = std::max(y[i], y[j]);
			if(minY == maxY) {
				discord -= RBNumGreaterThan(revTree, &maxY) *
          RBNumLessThan(revTree, &maxY);
			}
		}
	}

  // Clean up used memory
	RBTreeDestroy(tree);
	RBTreeDestroy(revTree);
	free(savedYsInds);
	free(x);
	free(y);

  // Use counts of concordant/discordant quadruples to find value of statistic
	long double c = 16 * concord - 8 * discord;
	long double d = (c < 0) ? -1 : 1;
	return NumericVector::create(d * expl(logl(d * c) - (logl(n) + logl(n - 1) +
    logl(n - 2) + logl(n - 3))));
}

/***
 * The same as TStarWeihsEtAlRCPP except computes the V-statistic form of the
 * statistic.
 *
 * @param xNumeric a NumericVector of values
 * @param yNumeric a NumericVector of values of the same length as xNumeric
 * @return the V-statistic computed upon the two input vectors.
 */
// [[Rcpp::export]]
Rcpp::NumericVector VTStarWeihsEtAlRCPP(NumericVector xNumeric,
                                       NumericVector yNumeric) {
  /***
   * Note, there are very few inline comments here as they would be almost
   * identical to those comments in the TStarWeihsEtAlRCPP function. See the
   * paper for differences in the two algorithms.
   */
	// Copy input data to something understood by the rbtree functions
  int n = xNumeric.size();
  double *x = (double*) malloc(sizeof(double) * n);
  double *y = (double*) malloc(sizeof(double) * n);
	for(int i = 0; i < n; i++) {
		x[i] = xNumeric[i];
		y[i] = yNumeric[i];
	}

  rb_red_blk_tree* tree =
    RBTreeCreate(DoubComp, DoubDest, InfoDest, DoubPrint, InfoPrint);
  long double concord = 0;
  long double discord = 0;
	int *savedYsInds = (int*) malloc(sizeof(int) * n);
  int totalYInTree = 0;
  int numSavedYs = 0;
  double lastX = 0;

	for(int i = 0; i < n; i++) {
		if(lastX == x[i] && i != 0) {
			savedYsInds[numSavedYs] = i;
			numSavedYs += 1;
		} else {
			for(int k = 0; k < numSavedYs; k++) {
				totalYInTree += 1;
				RBTreeInsert(tree, &(y[savedYsInds[k]]), 0);
			}
			lastX = x[i];
			numSavedYs = 1;
			savedYsInds[0] = i;
		}
		int top = RBNumGreaterThan(tree, &(y[i]));
		int bottom = RBNumLessThan(tree, &(y[i]));
		concord += .5 * ((bottom * (bottom - 1)) / 2) + .5 * ((top * (top - 1)) / 2)
      + .25 * (top + bottom);
		for(int j = i + 1; j < n; j++) {
			double minY = std::min(y[i], y[j]);
			double maxY = std::max(y[i], y[j]);

			int numLessMax = RBNumLessThan(tree, &maxY);
			int numGreaterMin = RBNumGreaterThan(tree, &minY);
			int top = RBNumGreaterThan(tree, &maxY);
			int bottom = RBNumLessThan(tree, &minY);

      int middle = 0;
			if(minY != maxY) {
				middle = numGreaterMin + numLessMax - totalYInTree;
			}

			concord += ((bottom * (bottom - 1)) / 2) + ((top * (top - 1)) / 2) +
        .5 * (top + bottom);

			if(minY != maxY) {
				discord += ((middle * (middle - 1)) / 2) + top * middle + top * bottom +
          middle * bottom;
				int numYMin = totalYInTree - numGreaterMin - bottom;
				int numYMax = totalYInTree - top - numLessMax;
				discord += numYMin * numGreaterMin + numYMax * (numLessMax - numYMin);
			}
		}
	}

	// Now run everything in reverse to get rid of
	// quadruples incorrectly counted as discordant
  rb_red_blk_tree* revTree = RBTreeCreate(DoubComp, DoubDest, InfoDest, DoubPrint, InfoPrint);
	numSavedYs = 0;
	lastX = 0;

	for(int i = n - 1; i > 0; i--) {
		if(lastX == x[i] && i != n - 1) {
			savedYsInds[numSavedYs] = i;
			numSavedYs += 1;
		} else {
			for(int k = 0; k < numSavedYs; k++) {
				RBTreeInsert(revTree, &(y[savedYsInds[k]]), 0);
			}
			lastX = x[i];
			numSavedYs = 1;
			savedYsInds[0] = i;
		}
		for(int j = i - 1; j >= 0; j--) {
			double minY = std::min(y[i], y[j]);
			double maxY = std::max(y[i], y[j]);
			if(minY == maxY) {
				discord -= RBNumGreaterThan(revTree, &maxY) * RBNumLessThan(revTree, &maxY);
			}
		}
	}

	RBTreeDestroy(tree);
	RBTreeDestroy(revTree);
	free(savedYsInds);
	free(x);
	free(y);

	long double c = 16 * concord - 8 * discord;
	long double d = (c < 0) ? -1 : 1;
	return NumericVector::create(d * expl(logl(d * c) - 4 * logl(n)));
}

/***
 * A helper function for TStarFastResample.
 */
Rcpp::NumericVector getSubset(NumericVector x, IntegerVector inds) {
  NumericVector subset(inds.size());
  for(int i = 0; i < inds.size(); i++) {
    subset[i] = x[inds[i]];
  }
  return subset;
}

/***
 * A helper function for TStarFastResample.
 */
IntegerVector integerSort(IntegerVector x) {
   IntegerVector y = clone(x);
   std::sort(y.begin(), y.end());
   return y;
}

/***
 * Estimates tau* by iteratively resampling small subsets of the data, computing
 * t* (see VTStarWeihsEtAlRCPP) across these subsets, and then averaging the
 * result.
 *
 * @param xNumeric a NumericVector of values
 * @param yNumeric a NumericVector of values of the same length as xNumeric
 * @param numResamples the number of subsets to draw from the data
 * @param sampleSize the size of subsets to draw from the data, should be
 *        less than the size of xNumeric (and thus also yNumeric)
 * @return the V-statistic computed upon the two input vectors.
 */
// [[Rcpp::export]]
Rcpp::NumericVector TStarFastResampleRCPP(NumericVector xNumeric,
                                          NumericVector yNumeric,
                                          int numResamples,
                                          int sampleSize) {
  IntegerVector allIndices(xNumeric.size());
  for(int i = 0; i < xNumeric.size(); i++) {
    allIndices[i] = i;
  }

  NumericVector toReturn = NumericVector::create(0);
  IntegerVector inds;
  for(int i = 0; i < numResamples; i++) {
      inds = integerSort(RcppArmadillo::sample(allIndices, sampleSize, false, NumericVector::create()));
      toReturn[0] += TStarHellerAndHellerRCPP(
        getSubset(xNumeric, inds),
        getSubset(yNumeric, inds)
      );
  }
  toReturn[0] = toReturn[0] / numResamples;
  return toReturn;
}

/***
 * A helper function that computes the 'a' function from
 *
 * Bergsma, Wicher; Dassios, Angelos. A consistent test of independence based on
 * a sign covariance related to Kendall's tau. Bernoulli 20 (2014), no. 2,
 * 1006–1028.
 */
int bergDassAFunc(double z1, double z2, double z3, double z4) {
  int returnVal = 0;
  returnVal += (std::max(z1, z3) < std::min(z2, z4)) ? 1 : 0;
  returnVal += (std::max(z2, z4) < std::min(z1, z3)) ? 1 : 0;
  returnVal += (std::max(z1, z2) < std::min(z3, z4)) ? -1 : 0;
  returnVal += (std::max(z3, z4) < std::min(z1, z2)) ? -1 : 0;
  return returnVal;
}

/***
 * Similar to TStarWeihsEtAlRCPP but implements a slow, naive algorithm for the
 * computation of the U/V-statistics. This allows for the comparision of the old
 * algorithm to the new.
 *
 * @param x see xNumeric from TStarWeihsEtAlRCPP
 * @param y see yNumeric from TStarWeihsEtAlRCPP
 * @param vStat a boolean that, if true, will result in the computation of the
 *        V-statistic form of t* and, if false, will result in the computation
 *        of the U-statistic form of t*
 */
// [[Rcpp::export]]
Rcpp::NumericVector TStarNaiveRCPP(NumericVector x, NumericVector y, bool vStat) {
  int i,j,k,l;
  long double sumVal = 0;
	long double z[4];
	int a,b,c,d,p;
	int n = x.size();

  if (!vStat) {
    for(i = 0; i < n - 3; i++) {
    	for(j = i + 1; j < n - 2; j++) {
  			for(k = j + 1; k < n - 1 ; k++) {
  				for(l = k + 1; l < n; l++) {
  					a = i; b = j; c = k; d = l;
  					if(x[a] > x[c]) {
  						p = a; a = c; c = p;
  					}
  					if(x[b] > x[c]) {
  						p = b; b = c; c = p;
  					}
  					if(x[a] > x[d]) {
  						p = a; a = d; d = p;
  					}
  					if(x[b] > x[d]) {
              p = b; b = d; d = p;
  					}

  					z[0] = y[i]; z[1] = y[j]; z[2] = y[k]; z[3] = y[l];
  					bubbleSort(z, 4);

  					if ((std::max(x[a], x[b]) < std::min(x[c], x[d])) &&
  						 ((std::max(y[a], y[b]) < std::min(y[c], y[d])) ||
  							(std::min(y[a], y[b]) > std::max(y[c], y[d])))) {
  						sumVal += 16;
  					} else if ((std::max(x[a], x[b]) < std::min(x[c], x[d])) && (z[1] < z[2])) {
  						sumVal -= 8;
  					}
  				}
  			}
  		}
  	}
  } else {
    for(i = 0; i < n; i++) {
      for(j = 0; j < n; j++) {
  			for(k = 0; k < n ; k++) {
  				for(l = 0; l < n; l++) {
            sumVal += bergDassAFunc(x[i], x[j], x[k], x[l]) *
              bergDassAFunc(y[i], y[j], y[k], y[l]);
  				}
  			}
  		}
  	}
  }

  int bar = (sumVal < 0) ? -1 : 1;
  if(vStat) {
    return NumericVector::create(bar * expl(
      logl(bar * sumVal) - (4.0 * logl(n))
    ));
  } else {
    return NumericVector::create(bar * expl(
      logl(bar * sumVal) - (logl(n) + logl(n - 1) + logl(n - 2) + logl(n - 3))
    ));
  }
}
