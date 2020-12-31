#ifndef __read_boundaries_file_h__
#define __read_boundaries_file_h__

#include <vector>
#include <string>
#include <utility>

std::pair<std::vector<double>, std::vector<double> > read_boundaries_file(const std::string& filename);
void verify_boundary_is_valid(const std::vector<double>& steps);

#endif
