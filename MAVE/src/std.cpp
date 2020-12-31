/*
 * Student License - for use by students to meet course requirements and
 * perform academic research at degree granting institutions only.  Not
 * for government, commercial, or other organizational use.
 *
 * std.cpp
 *
 * Code generation for function 'std'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "MAVEfast.h"
#include "std.h"
#include "MAVEfast_emxutil.h"

/* Function Definitions */
void b_std(const emxArray_real_T *varargin_1, emxArray_real_T *y)
{
  int n;
  int d;
  int ix;
  int nx;
  double xbar;
  int k;
  double r;
  double b_y;
  n = varargin_1->size[0];
  if (varargin_1->size[0] > 1) {
    d = varargin_1->size[0] - 1;
  } else {
    d = varargin_1->size[0];
  }

  ix = y->size[0] * y->size[1];
  y->size[0] = 1;
  y->size[1] = varargin_1->size[1];
  emxEnsureCapacity((emxArray__common *)y, ix, sizeof(double));
  if ((varargin_1->size[0] == 0) || (varargin_1->size[1] == 0)) {
    ix = y->size[0] * y->size[1];
    y->size[0] = 1;
    emxEnsureCapacity((emxArray__common *)y, ix, sizeof(double));
    nx = y->size[1];
    for (ix = 0; ix < nx; ix++) {
      y->data[y->size[0] * ix] = 0.0;
    }
  } else {
    for (nx = 0; nx + 1 <= varargin_1->size[1]; nx++) {
      ix = 0;
      xbar = varargin_1->data[varargin_1->size[0] * nx];
      for (k = 2; k <= n; k++) {
        ix++;
        xbar += varargin_1->data[ix + varargin_1->size[0] * nx];
      }

      xbar /= (double)n;
      ix = 0;
      r = varargin_1->data[varargin_1->size[0] * nx] - xbar;
      b_y = r * r;
      for (k = 2; k <= n; k++) {
        ix++;
        r = varargin_1->data[ix + varargin_1->size[0] * nx] - xbar;
        b_y += r * r;
      }

      b_y /= (double)d;
      y->data[nx] = b_y;
    }
  }

  nx = y->size[1];
  for (k = 0; k + 1 <= nx; k++) {
    y->data[k] = std::sqrt(y->data[k]);
  }
}

/* End of code generation (std.cpp) */
