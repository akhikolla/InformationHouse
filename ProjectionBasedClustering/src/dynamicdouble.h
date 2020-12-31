/*Copyright (C) Kristian Nybo, Jarkko Venna
 *
 *This software is released under the GNU Lesser General Public
 *License. See the included file LICENSE for details.*/

#ifndef DYNAMICDOUBLE_HH
#define DYNAMICDOUBLE_HH

#include <cstddef>

namespace dredviz {
class DynamicDouble
{
protected:
  const double originalValue;
  const double finalValue;
  double currentValue;
public:
  virtual ~ DynamicDouble ()
  {
  }
  DynamicDouble (double initialValue, double finalValue);

  double value () const
  {
    return currentValue;
  }
  virtual void update (size_t currentRound, size_t totalRounds);
};
}

#endif
