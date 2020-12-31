#include "read_boundaries_file.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <cassert>
#include <limits>

#include "string_utils.h"

using namespace std;

static bool is_monotone_increasing(const vector<double>& v)
{
    double prev = -numeric_limits<double>::infinity();
    for (int i = 0; i < (int)v.size(); ++i) {
        if (v[i] < prev) {
            return false;
        }
        prev = v[i];
    }
    return true;
}

pair<vector<double>, vector<double> > read_boundaries_file(const string& filename)
{
    string line;
    ifstream f(filename.c_str());
    if (!f.is_open()) {
        throw runtime_error("Unable to read input file '" + filename + "'");
    }
    f.exceptions(ifstream::failbit | ifstream::badbit);

    getline(f, line);
    vector<double> lower_bound_steps = read_comma_delimited_doubles(line);

    getline(f, line);
    vector<double> upper_bound_steps = read_comma_delimited_doubles(line);

    return pair<vector<double>, vector<double> >(lower_bound_steps, upper_bound_steps);
}

void verify_boundary_is_valid(const vector<double>& steps)
{
    if (!is_monotone_increasing(steps)) {
        throw runtime_error("Bound steps are not monotone increasing.");
    }
    if ((steps.size() > 0) && ((steps.front() < 0.0) || (steps.back() > 1.0))) {
        throw runtime_error("Steps must be in the range 0 to 1.");
    }
}

