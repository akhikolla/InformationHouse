#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
std::vector<int> subv_int(std::vector<int> &v, std::vector<int> idx) {
  std::vector<int> subv;
  for (int i = 0; i < idx.size(); i++) {
    subv.push_back(std::move(v[idx[i]]));
  }
  return subv;
}

// [[Rcpp::export]]
bool does_intersect_vi(std::vector<int> x, std::vector<int> &y) {
  bool result = false;
  unsigned int x_size = x.size();
  unsigned int y_size = y.size();
  for (unsigned int i = 0; i < x_size; i++) {
    for (unsigned int j = 0; j < y_size; j++) {
      if (x[i] == y[j]) {
        result = true;
        j = y_size;
        i = x_size;
      }
    }
  }
  return result;
}

// [[Rcpp::export]]
std::vector<int> noc_cpp(std::vector<std::vector<int>> x) {
  int x_size = x.size();
  std::vector<int> remain_idx(x_size);
  for (int k = 0; k < x_size; k++) {
    remain_idx[k] = k;
  }
  int i = 0;
  std::vector<int> u(1);
  u[0] = 0;

  std::vector<int> inter;
  while (i < x_size) {
    inter.clear();
    for (int j = 0; j < remain_idx.size(); j++) {
      if (!does_intersect_vi(x[i], x[remain_idx[j]])) {
        inter.push_back(j);
      }
    }
    // remain_idx = subv_int(remain_idx, inter);
    if (inter.size() > 0) {
      remain_idx = subv_int(remain_idx, inter);
      i = *std::min_element(remain_idx.begin(), remain_idx.end());
      u.push_back(i);
    } else {
      i = x_size;
    }
  }
  return u;
}
