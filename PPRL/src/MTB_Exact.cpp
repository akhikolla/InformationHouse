/*
 * MTB_Exact.cpp
 *
 *  Created on: 04.04.2017
 *      Author: schnell-42
 */

#include "MTB_Exact.h"

/**
 *
 * @param var1
 * @param var2
 * @return
 */
bool exact(string var1, string var2){
	return var1.compare(var2);
}

string MTB_ExactAlgorithm::getName(){
	return "Exact";
}


double MTB_ExactAlgorithm::getRelativeValue(string var1, string var2)
	{
		if (var1.compare(var2)==0)
		{
			return 1.0;
		}
		return 0.0;
	}

double MTB_ExactAlgorithm::getAbsoluteValue(string var1, string var2)
	{
				return getRelativeValue(var1,var2);
	}
/**
 *
 * @param var1
 * @param var2
 * @return
 */
//bool exactCapitalLetters(string var1, string var2){
//	return StringToUpper(var1).compare(StringToUpper(var2));
//		}

string MTB_ExactCapitalLettersAlgorithm::getName(){
	return "ExactCapitlalLetter";
}

double MTB_ExactCapitalLettersAlgorithm::getRelativeValue(string var1, string var2)
	{
		if (StringToUpper(var1).compare(StringToUpper(var2))==0)
		{
			return 1.0;
		}
		return 0.0;
	}

double MTB_ExactCapitalLettersAlgorithm::getAbsoluteValue(string var1, string var2)
	{

		return getRelativeValue( var1,  var2);
	}

/**
 *
 * @param var
 * @return
 */
string StringToUpper(string var)
{
    std::transform(var.begin(), var.end(), var.begin(), ::toupper);

    return var;
}
