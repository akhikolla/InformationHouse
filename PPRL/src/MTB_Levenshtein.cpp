/*
 * levenshtein.cpp
 *
 *  Created on: 30.01.2017
 *      Author: rukasz
 */
#include "MTB_Levenshtein.h"

/**
 * Calculates Levenshtein distance
 *
 * @param s1
 * @param s2
 * @return
 */
int levenshteinDistance(const std::string &s1, const std::string &s2)
{

	// To change the type this function manipulates and returns, change
	// the return type and the types of the two variables below.
	int s1len = s1.size();
	int s2len = s2.size();

	auto column_start = (decltype(s1len))1;
	auto column = new decltype(s1len)[s1len + 1];
	std::iota(column + column_start, column + s1len + 1, column_start);

	for (auto x = column_start; x <= s2len; x++) {
		column[0] = x;
		auto last_diagonal = x - column_start;
		for (auto y = column_start; y <= s1len; y++) {
			auto old_diagonal = column[y];
			auto possibilities = { // deletion, insertion, substitution
				column[y] + 1,
				column[y - 1] + 1,
				last_diagonal + (s1[y - 1] == s2[x - 1]? 0 : 1)
			};
			column[y] = std::min(possibilities);
			last_diagonal = old_diagonal;
		}
	}
	auto result = column[s1len];
	delete[] column;
	return result;
}

string MTB_LevenshteinAlgorithm::getName(){
	return "Levenshtein";
}

double MTB_LevenshteinAlgorithm::getRelativeValue(string o1, string o2){
			return calcFromDistanceToInterval(o1, o2, levenshteinDistance(o1,o2));
	}

double MTB_LevenshteinAlgorithm::getAbsoluteValue(string o1, string o2){
			return levenshteinDistance(o1,o2);
	}


/**
	 * used for Levenshtein and Damerau-Levensthein
	 *
	 *
	 * @param str1
	 * @param str2
	 * @param result
	 * @return
	 */
	double calcFromDistanceToInterval(string str1, string str2, double result)
	{
		double avg = averageLength(str1, str2);
		if (avg > result)
		{
			return (1 - (result / avg));

		}
		return 0;
	}

		/**
	 * Returns average length of 2 given Strings
	 *
	 * @param str1
	 *            String 1
	 * @param str2
	 *            String 2
	 * @return Average length
	 */
	double averageLength(string str1, string str2)
	{
		vector<string> strvec;
		strvec.push_back(str1);
		strvec.push_back(str2);
		return averageLengthVec(strvec);
	}

	/**
	 * Returns average length of the Strings in the given String Array
	 *
	 * @param strings
	 *            Vector of strings
	 * @return Average length
	 */
	double averageLengthVec(vector<string> strings)
	{
		if (strings.size() > 0)
		{
			double result = 0;
			for (string string : strings)
			{
				if (string.length()>0 )
				{
					result += string.length();
				}
			}
			return result / strings.size();
		}
		return 0;
	}
