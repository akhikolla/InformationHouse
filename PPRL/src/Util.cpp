/*
 * Util.cpp
 *
 *  Created on: 05.07.2017
 *      Author: schnell-42
 */

#include "Util.h"

// void replaceAll(string& str, const string& from, const string& to) {
// 	size_t start_pos = 0;
// 	while ((start_pos = str.find(from, start_pos)) != string::npos) {
// 		str.replace(start_pos, from.length(), to);
// 		start_pos += to.length(); // Handles case where 'to' is a substring of 'from'
// 	}
// 	//return str;
// }

// void replaceSharpS(string& str) {
// 	//if (str.length() == 0) {
// 	//return "";
// 	//}
// 	replaceAll(str, "\u00df", "ss");
// 	//return str;
// }

// void replaceUmlauts(string& str) {
//
// 	replaceAll(str, "\u00E4", "ae");
// 	replaceAll(str, "\u00F6", "oe");
// 	replaceAll(str, "\u00FC", "ue");
// 	replaceAll(str, "\u00C4", "Ae");
// 	replaceAll(str, "\u00D6", "Oe");
// 	replaceAll(str, "\u00DC", "Ue");
//
// }


string replaceEnding(string fullString, string ending, char repl) {
	bool endsWith = false;
	unsigned i = 0;
	if (fullString.length() >= ending.length()) {
		while ((fullString[fullString.length() - ending.length() + i]
				== ending[i])
				&& (fullString.length() - ending.length() + i
						<= fullString.length())) {
			endsWith = true;
			i++;
		}
		if (endsWith && (i - 1 == ending.length()))
			return fullString.substr(0, fullString.length() - ending.length())
					+ repl;
	}
	return fullString;
}

vector<int> convertToSingleFrequencies(vector<int> patternFrequencies) {
	// return null if null or no valid number of patterns

	int singleFrequNo = (int) round(log2(patternFrequencies.size()));
	vector<int> singleFrequencies(singleFrequNo);
	if (pow(2, singleFrequNo) != patternFrequencies.size()) {
		return singleFrequencies;
	}
	for (int i = 0; i < singleFrequNo; i++){singleFrequencies[i]= 0;}
	for (int i = patternFrequencies.size(); i > 1; i--) {
		int pos = patternFrequencies.size() / 2;
		int j = i - 1;
		while ((j > 0) && (pos >= 1)) {
			if ((j - pos) > -1) {
				int sfPos = (int) round(log2(pos));
				singleFrequencies[sfPos] += patternFrequencies[i - 1];
				j -= pos;
			}
			pos = pos > 1 ? pos / 2 : 0;
		}
	}
	return singleFrequencies;
}

 int sum(vector<int> values)
	{
		int sum = 0;
		if (values.size() != 0)
		{
			for (int value : values)
			{
				sum += value;
			}
		}
		return sum;
	}

