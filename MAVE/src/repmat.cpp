/*
 * Student License - for use by students to meet course requirements and
 * perform academic research at degree granting institutions only.  Not
 * for government, commercial, or other organizational use.
 *
 * repmat.cpp
 *
 * Code generation for function 'repmat'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "MAVEfast.h"
#include "repmat.h"
#include "MAVEfast_emxutil.h"

/* Function Definitions */
void b_repmat(const emxArray_real_T *a, double varargin_2, emxArray_real_T *b)
{
  int outsize_idx_0;
  int jtilecol;
  boolean_T p;
  int ibtile;
  int k;
  outsize_idx_0 = a->size[0];
  jtilecol = b->size[0] * b->size[1];
  b->size[0] = outsize_idx_0;
  b->size[1] = (int)varargin_2;
  emxEnsureCapacity((emxArray__common *)b, jtilecol, sizeof(double));
  if (!(a->size[0] == 0)) {
    if (outsize_idx_0 == 0) {
      p = true;
    } else if ((int)varargin_2 == 0) {
      p = true;
    } else {
      p = false;
    }

    if (!p) {
      outsize_idx_0 = a->size[0];
      for (jtilecol = 1; jtilecol <= (int)varargin_2; jtilecol++) {
        ibtile = (jtilecol - 1) * outsize_idx_0;
        for (k = 0; k + 1 <= outsize_idx_0; k++) {
          b->data[ibtile + k] = a->data[k];
        }
      }
    }
  }
}

void repmat(const emxArray_real_T *a, double varargin_1, emxArray_real_T *b)
{
  int outsize_idx_1;
  int ibmat;
  boolean_T p;
  int itilerow;
  outsize_idx_1 = a->size[1];
  ibmat = b->size[0] * b->size[1];
  b->size[0] = (int)varargin_1;
  b->size[1] = outsize_idx_1;
  emxEnsureCapacity((emxArray__common *)b, ibmat, sizeof(double));
  if (!(a->size[1] == 0)) {
    if ((int)varargin_1 == 0) {
      p = true;
    } else if (outsize_idx_1 == 0) {
      p = true;
    } else {
      p = false;
    }

    if (!p) {
      for (outsize_idx_1 = 0; outsize_idx_1 + 1 <= a->size[1]; outsize_idx_1++)
      {
        ibmat = outsize_idx_1 * (int)varargin_1;
        for (itilerow = 1; itilerow <= (int)varargin_1; itilerow++) {
          b->data[(ibmat + itilerow) - 1] = a->data[outsize_idx_1];
        }
      }
    }
  }
}

/* End of code generation (repmat.cpp) */
