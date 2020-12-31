/*
 * MTB_Jaro.cpp
 *
 *  Created on: 10.04.2017
 *      Author: schnell-42
 */

#include "MTB_Jaro.h"

#include <algorithm>
double jaroDistance(const std::string& a, const std::string& b)
{
    // Register strings length.
    int aLength(a.size());
    int bLength(b.size());

    // If one string has null length, we return 0.
    if (aLength == 0 || bLength == 0)
    {
        return 0.0;
    }

    // Calculate max length range.
    int maxRange(std::max(0, std::max(aLength, bLength) / 2 - 1));

    // Creates 2 vectors of integers.
    std::vector<bool> aMatch(aLength, false);
    std::vector<bool> bMatch(bLength, false);

    // Calculate matching characters.
    int matchingCharacters(0);
    for (int aIndex(0); aIndex < aLength; ++aIndex)
    {
        // Calculate window test limits (limit inferior to 0 and superior to bLength).
        int minIndex(std::max(aIndex - maxRange, 0));
        int maxIndex(std::min(aIndex + maxRange + 1, bLength));

        if (minIndex >= maxIndex)
        {
            // No more common character because we don't have characters in b to test with characters in a.
            break;
        }

        for (int bIndex(minIndex); bIndex < maxIndex; ++bIndex)
        {
            if (!bMatch.at(bIndex) && a.at(aIndex) == b.at(bIndex))
            {
                // Found some new match.
                aMatch[aIndex] = true;
                bMatch[bIndex] = true;
                ++matchingCharacters;
                break;
            }
        }
    }

    // If no matching characters, we return 0.
    if (matchingCharacters == 0)
    {
        return 0.0;
    }

    // Calculate character transpositions.
    std::vector<int> aPosition(matchingCharacters, 0);
    std::vector<int> bPosition(matchingCharacters, 0);

    for (int aIndex(0), positionIndex(0); aIndex < aLength; ++aIndex)
    {
        if (aMatch.at(aIndex))
        {
            aPosition[positionIndex] = aIndex;
            ++positionIndex;
        }
    }

    for (int bIndex(0), positionIndex(0); bIndex < bLength; ++bIndex)
    {
        if (bMatch.at(bIndex))
        {
            bPosition[positionIndex] = bIndex;
            ++positionIndex;
        }
    }

    // Counting half-transpositions.
    int transpositions(0);
    for (int index(0); index < matchingCharacters; ++index)
    {
        if (a.at(aPosition.at(index)) != b.at(bPosition.at(index)))
        {
            ++transpositions;
        }
    }

    // Calculate Jaro distance.
    return (
        JARO_WEIGHT_STRING_A * matchingCharacters / aLength +
        JARO_WEIGHT_STRING_B * matchingCharacters / bLength +
        JARO_WEIGHT_TRANSPOSITIONS * (matchingCharacters - transpositions / 2) / matchingCharacters
    );
}

string MTB_JaroAlgorithm::getName(){
	return "Jaro";
}

double MTB_JaroAlgorithm::getRelativeValue(string var1, string var2)
	{
		return jaroDistance(var1, var2);
	}

double MTB_JaroAlgorithm::getAbsoluteValue(string var1, string var2)
	{
		return jaroDistance(var1, var2);
	}


double jaroWinklerDistance(const std::string&a, const std::string& b)
{
    // Calculate Jaro distance.
    double distance(jaroDistance(a, b));

    if (distance > JARO_WINKLER_BOOST_THRESHOLD)
    {
        // Calculate common string prefix.
        int commonPrefix(0);
        for (int index(0), indexEnd(min((unsigned long)min(a.size(), b.size()),
                                    (unsigned long)JARO_WINKLER_PREFIX_SIZE)); index < indexEnd; ++index)
        {
            if (a.at(index) == b.at(index))
            {
                ++commonPrefix;
            }
            else
            {
                break;
            }
        }

        // Calculate Jaro-Winkler distance.
        distance += JARO_WINKLER_SCALING_FACTOR * commonPrefix * (1.0 - distance);
    }

    return distance;
}

string MTB_JaroWinklerAlgorithm::getName(){
	return "JaroWinkler";
}

double MTB_JaroWinklerAlgorithm::getRelativeValue(string var1, string var2)
	{
		return jaroWinklerDistance(var1, var2);
	}

double MTB_JaroWinklerAlgorithm::getAbsoluteValue(string var1, string var2)
	{
		return jaroWinklerDistance(var1, var2);
	}
