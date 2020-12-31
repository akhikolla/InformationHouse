#include <iomanip>
#include <iostream>
#include <fstream>
#include <vector>
#include <stdexcept>
#include <cassert>
#include <algorithm>
#include <numeric>

#include "stdio.h"
#include <iterator>
#include <fstream>
#include <string>
#include <math.h>
#include <set>
#include <cmath>
#include <random>


#include "two_sided_noncrossing_probability.h"
#include "read_boundaries_file.h"
#include "string_utils.h"

using namespace std;

/*========================================================================*/
/******************************************************************************/
double cont_ks_distribution(long n){


    pair<vector<double>, vector<double> > bounds = read_boundaries_file("Boundary_Crossing_Time.txt");
    const vector<double>& lower_bound_steps = bounds.first;
    const vector<double>& upper_bound_steps = bounds.second;

    bool use_fft = true;

    return 1.0 - ecdf_noncrossing_probability(n, lower_bound_steps, upper_bound_steps, use_fft);

}


