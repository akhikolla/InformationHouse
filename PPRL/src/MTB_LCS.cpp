/*
 * LCS.cpp
 *
 *  Longest common subsequence
 *
 *  Created on: Feb 17, 2017
 *      Author: rukasz
 *  source: LCS (Hirschberg, 1977)
 */
#include "MTB_LCS.h"

/* Returns length of LCS for X[0..m-1], Y[0..n-1] */
int lcsHelper(string X, string Y) {
	int m = X.size();
	int n = Y.size();

	vector<vector<int>> L(m + 1, vector<int>(n + 1));

	/* Following steps build L[m+1][n+1] in bottom up fashion. Note
	 that L[i][j] contains length of LCS of X[0..i-1] and Y[0..j-1] */
	for (int i = 0; i <= m; i++) {
		for (int j = 0; j <= n; j++) {
			if (i == 0 || j == 0)
				L[i][j] = 0;
			else if (X[i - 1] == Y[j - 1])
				L[i][j] = L[i - 1][j - 1] + 1;
			else
				L[i][j] = max(L[i - 1][j], L[i][j - 1]);
		}
	}
	/* L[m][n] contains length of LCS for X[0..n-1] and Y[0..m-1] */
	return L[m][n];
}

/* Utility function to get max of 2 integers */
int max(int a, int b) {
	return (a > b) ? a : b;
}


double calculateCoeff(string str1, string str2, double result)
	{
		if ((str1.length() > 0) && (str2.length() > 0))
		{
			return ((2 * result) / (str1.length() + str2.length()));
		}
		if (str1.length() > 0)
		{
			return (result / str1.length());
		}
		return (result / str2.length());
	}

double MTB_LongestCommonSubsequenceAlgorithm::getRelativeValue(string o1, string o2){

	if ((o1.length() > 0) && (o2.length() > 0))
			{
				return calculateCoeff(o1, o2, (double)lcsHelper(o1, o2));
			}

	return (double)lcsHelper(o1, o2);
}

double MTB_LongestCommonSubsequenceAlgorithm::getAbsoluteValue(string o1, string o2){

	return (double)lcsHelper(o1, o2);
}

string MTB_LongestCommonSubsequenceAlgorithm::getName(){return "Longest common subsequence";}
