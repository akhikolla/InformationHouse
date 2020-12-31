/*Copyright (C) Kristian Nybo, Jarkko Venna
 *
 *This software is released under the GNU Lesser General Public
 *License. See the included file LICENSE for details.*/

#ifndef OPTIMIZATIONSTEPBATCH_HH
#define OPTIMIZATIONSTEPBATCH_HH

#include "optimizationstep.h"

namespace dredviz {
class OptimizationStepBatch:public OptimizationStep
{
public:
  virtual void setIterationsPerStep (size_t number) = 0;
};

}
#endif
