
#include "NDArrayUtils.h"

int32_t maxAbsElement(const std::vector<int32_t>& r)
{
  int32_t m = 0;
  for (size_t i = 0; i < r.size(); ++i)
  {
    m = std::max(m, abs(r[i]));
  }
  return m;
}

std::vector<int32_t> diff(const std::vector<uint32_t>& x, const std::vector<uint32_t>& y)
{
  size_t size = x.size();
  assert(size == y.size());

  std::vector<int32_t> result(size);

  for (size_t i = 0; i < x.size(); ++i)
  {
    result[i] = x[i] - y[i];
  }
  return result;
}




