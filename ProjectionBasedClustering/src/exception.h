/*Copyright (C) Kristian Nybo, Jarkko Venna
 *
 *This software is released under the GNU Lesser General Public
 *License. See the included file LICENSE for details.*/

#ifndef EXCEPTION_HH
#define EXCEPTION_HH

#include <string>

namespace dredviz {
class Exception
{
private:
  std::string errMsg;

public:
  Exception (std::string errMsg);
  virtual ~ Exception ()
  {
  };

  const std::string & getErrMsg () const;
};

}
#endif
