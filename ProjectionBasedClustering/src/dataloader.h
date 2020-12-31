/*Copyright (C) Kristian Nybo, Jarkko Venna
 *
 *This software is released under the GNU Lesser General Public
 *License. See the included file LICENSE for details.*/

#ifndef DATALOADER_HH
#define DATALOADER_HH

#include "dataset.h"
#include "exception.h"

namespace dredviz {
class DataLoader
{
public:

  // Returns a DataSet object constructed from the data loaded by this loader.

  virtual void loadData (DataSet & target) = 0;

  virtual ~ DataLoader ()
  {
  };

  // Exceptions:

class InputError:public Exception
  {
  public:
    InputError (std::string errMsg):Exception (errMsg)
    {
    };
  };
};
}

#endif
