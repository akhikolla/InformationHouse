/*Copyright (C) Kristian Nybo, Jarkko Venna
 *
 *This software is released under the GNU Lesser General Public
 *License. See the included file LICENSE for details.*/

#ifndef RANDOMDATAGENERATOR_HH
#define RANDOMDATAGENERATOR_HH

#include "dataloader.h"
#include "dataset.h"
#include "datamatrix.h"

namespace dredviz {
class RandomDataGenerator:public DataLoader
{
  double range;
  DataMatrix matrix;

public:

  RandomDataGenerator (size_t dimension, size_t numberOfPoints,
                       double range);

  void loadData (DataSet & target);
};
}

#endif
