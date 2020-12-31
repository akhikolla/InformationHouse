/*Copyright (C) Kristian Nybo, Jarkko Venna
 *
 *This software is released under the GNU Lesser General Public
 *License. See the included file LICENSE for details.*/

#include "dynamicdouble.h"

namespace dredviz {
DynamicDouble::DynamicDouble (double initialValue, double finalValue):
    originalValue (initialValue),
    finalValue (finalValue),
    currentValue (initialValue)
{
}


void
DynamicDouble::update (size_t currentRound, size_t totalRounds)
{
  currentValue = originalValue +
                 ((finalValue - originalValue) / (double) totalRounds)
                     * (double) currentRound;
}
}
