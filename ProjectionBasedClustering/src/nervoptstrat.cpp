/*Copyright (C) Kristian Nybo, Jarkko Venna
 *
 *This software is released under the GNU Lesser General Public
 *License. See the included file LICENSE for details.*/

#include "nervoptstrat.h"

#include <sstream>

namespace dredviz {
void
NeRVOptStrat::optimize (const DistanceMatrix & origDist,
                        DataMatrix & initialProjData,
                        OptimizationStepBatch & optStep,
                        CostFunction & costFunc, std::ostream & feedback)
{

  optStep.setIterationsPerStep (STEPS_PER_ROUND);
  costFunc.updateDynamicParameters (TOTAL_ROUNDS, TOTAL_ROUNDS,
                                    initialProjData);

  feedback << "Initial cost: " << costFunc.evaluate (initialProjData)
  << std::endl;

  for (size_t roundsLeft = TOTAL_ROUNDS; roundsLeft > 0; roundsLeft--)
  {
    optStep.updateDynamicParameters (TOTAL_ROUNDS - roundsLeft, TOTAL_ROUNDS);

    costFunc.updateDynamicParameters (TOTAL_ROUNDS - roundsLeft, TOTAL_ROUNDS,
                                      initialProjData);

    feedback << "Starting round " << TOTAL_ROUNDS - roundsLeft << "...\n";

    optStep.perform (initialProjData);

    feedback << "Done.\n" << std::endl;
  }

  feedback << "Starting final round, performing " << STEPS_ON_LAST_ROUND
  << " optimization steps.\n";

  optStep.updateDynamicParameters (TOTAL_ROUNDS, TOTAL_ROUNDS);

  costFunc.updateDynamicParameters (TOTAL_ROUNDS, TOTAL_ROUNDS,
                                    initialProjData);

  optStep.setIterationsPerStep (STEPS_ON_LAST_ROUND);
  optStep.perform (initialProjData);

  feedback << "Done.\n" << std::endl;
}
}
