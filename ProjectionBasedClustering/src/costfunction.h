/*Copyright (C) Kristian Nybo, Jarkko Venna
 *
 *This software is released under the GNU Lesser General Public
 *License. See the included file LICENSE for details.*/

#ifndef COSTFUNCTION_HH
#define COSTFUNCTION_HH

#include "datamatrix.h"
#include <string>

namespace dredviz {
class CostFunction
{
protected:
  //  Matrix gradient;
public:

  /* Evaluates the gradient at data and returns the result in target. */

  virtual double getGradient (const DataMatrix & data, DataMatrix & target) =
    0;
  /* Evaluates the cost function at data. */

  virtual double evaluate (const DataMatrix & data) = 0;

  virtual ~ CostFunction ()
  {
  }

  /* Reports the parameters of the cost function as a string. */
  virtual void reportParameters(std::string& target) = 0;

  /*Updates any parameters that change as the optimization progresses */

  virtual void updateDynamicParameters (size_t currentRound,
                                        size_t totalRounds,
                                        const DataMatrix & projData) = 0;

  virtual void updateDataRepresentation (const DataMatrix & newData) = 0;

};
}

#endif
