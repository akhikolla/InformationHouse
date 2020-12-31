#include <string>
#include "datamatrix.h"

namespace dredviz {
class Recorder
{
  private:
    std::string filename_stem;
    size_t counter;

  public:
    void record(const DataMatrix& data);

    Recorder(std::string filename_stem)
      : filename_stem(filename_stem), counter(0)
    {};
};
}
