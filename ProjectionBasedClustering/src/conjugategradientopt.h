/*Copyright (C) Kristian Nybo, Jarkko Venna
 *
 *This software is released under the GNU Lesser General Public
 *License. See the included file LICENSE for details.*/

#ifndef CONJUGATEGRADIENTOPT_HH
#define CONJUGATEGRADIENTOPT_HH

#include "optimizationstepbatch.h"
#include "distancematrix.h"
#include "dataset.h"
#include "costfunction.h"
#include "linesearch.h"
#include "recorder.h"

#include <iostream>



namespace dredviz {
#define DEFAULT_ITERATIONS 5
class ConjugateGradientOpt:public OptimizationStepBatch
{
private:
  size_t iterationsPerStep;
  CostFunction & costFunc;
  LineSearch & linesearch;

  std::ostream & feedback;

  double previousStepSize;

  bool record; // export current state of data after each step?

  Recorder recorder;

public:
  ConjugateGradientOpt (CostFunction & costFunc, LineSearch & lineSearch,
                        std::ostream & feedback);

  ConjugateGradientOpt (CostFunction & costFunc, LineSearch & lineSearch,
                        std::ostream & feedback, std::string filename_stem);

  void perform (DataMatrix & projData);

  void updateDynamicParameters (size_t currentRound, size_t totalRounds);

  void setIterationsPerStep (size_t number);
};
}

#endif
