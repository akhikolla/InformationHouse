#include <algorithm>
#include <cmath>
#include <cassert>
#include <iostream>
#include <stdexcept>
#include <numeric>
#include <algorithm>
#include <sstream>
#include "fftwconvolver.h"
#include "aligned_mem.h"
#include "common.h"
#include "poisson_pmf.h"
#include "string_utils.h"
/*************/
#include <Rcpp.h>
using namespace Rcpp;
/*************/

using namespace std;

enum BoundType {H_STEP, G_STEP, END};

struct Bound {
    double location;
    BoundType tag;
};

static bool operator<(Bound b0, Bound b1)
{
    return (b0.location < b1.location);
}

static vector<Bound> join_all_bounds(const vector<double>& h_steps, const vector<double>& g_steps)
{
    assert(h_steps.size() >= g_steps.size());

    vector<Bound> bounds;
    bounds.reserve(h_steps.size()+g_steps.size()+1);
    Bound b;

    for (int i = 0; i < (int)h_steps.size(); ++i) {
        b.location = h_steps[i];
        b.tag = H_STEP;
        bounds.push_back(b);
    }

    for (int i = 0; i < (int)g_steps.size(); ++i) {
        b.location = g_steps[i];
        b.tag = G_STEP;
        bounds.push_back(b);
    }

    sort(bounds.begin(), bounds.end());

    b.location = 1.0;
    b.tag = END;
    bounds.push_back(b);

    return bounds;
}

static bool lower_and_upper_boundaries_cross(const vector<double>& g_steps, const vector<double>& h_steps)
{
    if (g_steps.size() > h_steps.size()) {
        Rcpp::Rcout << "The lower and upper boundaries cross: g(1) > h(1).\n";
        return true;
    }
    for (size_t i = 0; i < g_steps.size(); ++i) {
        if (g_steps[i] < h_steps[i]) {
            Rcpp::Rcout << "The lower and upper boundaries cross! i=" << i << ".\n";
            return true;
        }
    }
    return false;
}


vector<double> poisson_process_noncrossing_probability(double intensity, const vector<double>& g_steps, const vector<double>& h_steps, bool use_fft)
{
    if (lower_and_upper_boundaries_cross(g_steps, h_steps)) {
        return vector<double>();
    }

    vector<Bound> bounds = join_all_bounds(h_steps, g_steps);
    int n = h_steps.size();
    DoubleBuffer<double> buffers(n+1, 0.0);
    buffers.get_src()[0] = 1.0;
    FFTWConvolver fftconvolver(n+1);
    PoissonPMFGenerator pmfgen(n+1);
    double* pmf = allocate_aligned_doubles(n+1);
    int h_step_count = 0;
    int g_step_count = 0;
    double prev_location = 0.0;

    for (unsigned int i = 0; i < bounds.size(); ++i) {
        int cur_size = h_step_count - g_step_count + 1;

        double lambda = intensity*(bounds[i].location-prev_location);
        pmfgen.compute_pmf(cur_size, lambda, pmf);

        if (use_fft) {
            fftconvolver.convolve_same_size(cur_size, pmf, &buffers.get_src()[g_step_count], &buffers.get_dest()[g_step_count]);
        } else {
            convolve_same_size(cur_size, pmf, &buffers.get_src()[g_step_count], &buffers.get_dest()[g_step_count]);
        }

        BoundType tag = bounds[i].tag;
        if (tag == H_STEP) {
            ++h_step_count;
            buffers.get_dest()[h_step_count] = 0.0;
            buffers.get_src()[h_step_count] = 0.0; // Not strictly necessary. This just keeps the arrays cleaner when printing.
        } else if (tag == G_STEP) {
            buffers.get_dest()[g_step_count] = 0.0;
            buffers.get_src()[g_step_count] = 0.0; // Not strictly necessary. This just keeps the arrays cleaner when printing.
            ++g_step_count;
        } else {
            assert(tag == END);
            break;
        }
        prev_location = bounds[i].location;
        buffers.flip();
    }

    free_aligned_mem(pmf);
    return buffers.get_dest();
}



double ecdf_noncrossing_probability(int n, const vector<double>& g_steps, const vector<double>& h_steps, bool use_fft)
{
    if ((int)g_steps.size() > n) {
        stringstream ss;
        ss << "Empirical CDF must cross lower boundary g(t) since g(1)==" << g_steps.size() << " > n and the number of samples is n." << endl;
        throw runtime_error(ss.str());
    }
    vector<double> processed_h_steps(n, 0.0);
    if (h_steps.size() == 0) {
        // Special case, only the lower bound is specified.
        // We treat this as an implicit upper bound satisfying h(t) = n for all t.
    } else {
        if (lower_and_upper_boundaries_cross(g_steps, h_steps)) {
            return 0.0;
        }
        if ((int)h_steps.size() < n) {
            stringstream ss;
            ss << "Empirical CDF must cross lower boundary g(t) since h(1)==" << h_steps.size() << " > n and the number of samples is n. h_steps:" << endl;
            throw runtime_error(ss.str() + vector_to_string(h_steps));
        }
        copy(h_steps.begin(), h_steps.begin() + n, processed_h_steps.begin());
    }
    vector<double> poisson_nocross_probs = poisson_process_noncrossing_probability(n, g_steps, processed_h_steps, use_fft);
    return poisson_nocross_probs[n] / poisson_pmf(n, n);
}

