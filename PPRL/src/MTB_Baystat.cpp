/*
 * Baystat.cpp
 *
 *  Created on: Feb 17, 2017
 *      Author: rukasz
 */

#include "MTB_Baystat.h"

int maxLength(string str1, string str2)
	{
		if (str1.length() > str2.length())
		{
			return str1.length();
		}
		return str2.length();
	}
string repeatedString(string str, int count)
	{
		for (int i = 0; i < count; i++)
		{
			str = str.append(str);
		}
		return str;
	}

string MTB_Baystat::getName()
	{
		return "Baystat";
	}


double MTB_Baystat::getRelativeValue(string string1, string string2)
	{
	int minSL, rangeLeft, rangeRight, maxL;

		if ((string1.size() == 0) || (string2.size() == 0))
		{
			return 0;
		}
		if ((maxL = maxLength(string1, string2)) >= 7)
		{
			minSL = this->minSubstringLengthLongerEqual7;
			rangeLeft = this->rangeLeftLongerEqual7;
			rangeRight = this->rangeRightLongerEqual7;
		}
		else
		{
			minSL = this->minSubstringLengthSmaller7;
			rangeLeft = this->rangeLeftSmaller7;
			rangeRight = this->rangeRightSmaller7;
		}
		// match
		unsigned index = 0;
		int lengthCS = 0;
		string substring = "";
		string range = "";
		do
		{
			if ((index + minSL) > string1.length())
			{
				return (lengthCS / (double) maxL);
			}
			int rangeBegin = 0;
			int rangeLeftLength = 0;
			int rangeRightLength = 0;
			int hitLength = 0;
			int i = 1;
			substring = string1.substr(index, minSL);
			rangeBegin = index - rangeLeft;
			if (rangeBegin < 0)
			{
				rangeBegin = 0;
				rangeLeftLength = index;
			}
			else
			{
				rangeLeftLength = rangeLeft;
			}
			if ((unsigned)(index + minSL + rangeRight) >= string2.length())
			{
				rangeRightLength = string2.length() - index - minSL;
			}
			else
			{
				rangeRightLength = rangeRight;
			}
			if ((rangeBegin + rangeLeftLength + minSL + rangeRightLength) > rangeBegin)
			{
				range = string2.substr(rangeBegin, rangeLeftLength + minSL + rangeRightLength);
			}
			else
			{
				range = "";
			}
			string flaggedString2 = "";
			while (((signed)range.find(substring) > -1) && ((index + i) <= string1.length()))
			{
				flaggedString2 = string2;
				string replaceWith = repeatedString("#", substring.length());
				ReplaceAllSubstr(flaggedString2, substring, replaceWith);
				hitLength = substring.length();
				if ((unsigned)(index + minSL + i) <= string1.length())
				{
					substring = string1.substr(index,  minSL + i);
				}
				if (!(unsigned)((index + minSL + rangeRightLength + 1) > string2.length()))
				{
					rangeRightLength += 1;
				}
				if ((unsigned)(rangeBegin + rangeLeftLength + minSL + rangeRightLength) <= string2.length())
				{
					range = string2.substr(rangeBegin, rangeLeftLength + minSL + rangeRightLength);
				}
				i++;
			}
			if (hitLength > 0)
			{
				string2 = flaggedString2;
			}
			lengthCS += hitLength;
			index += i;
		}
		while (true);
	}

double MTB_Baystat::getAbsoluteValue(string string1, string string2){
	return getRelativeValue(string1,string2);
}
