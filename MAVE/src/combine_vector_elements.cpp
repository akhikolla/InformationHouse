/*
 * Student License - for use by students to meet course requirements and
 * perform academic research at degree granting institutions only.  Not
 * for government, commercial, or other organizational use.
 *
 * combine_vector_elements.cpp
 *
 * Code generation for function 'combine_vector_elements'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "MAVEfast.h"
#include "combine_vector_elements.h"
#include "MAVEfast_emxutil.h"

/* Function Definitions */
void b_combine_vector_elements(const emxArray_real_T *x, emxArray_real_T *y)
{
  int vlen;
  int i;
  int xoffset;
  double s;
  int k;
  vlen = y->size[0] * y->size[1];
  y->size[0] = 1;
  y->size[1] = x->size[1];
  emxEnsureCapacity((emxArray__common *)y, vlen, sizeof(double));
  if ((x->size[0] == 0) || (x->size[1] == 0)) {
    vlen = y->size[0] * y->size[1];
    y->size[0] = 1;
    emxEnsureCapacity((emxArray__common *)y, vlen, sizeof(double));
    i = y->size[1];
    for (vlen = 0; vlen < i; vlen++) {
      y->data[y->size[0] * vlen] = 0.0;
    }
  } else {
    vlen = x->size[0];
    for (i = 0; i + 1 <= x->size[1]; i++) {
      xoffset = i * vlen;
      s = x->data[xoffset];
      for (k = 2; k <= vlen; k++) {
        s += x->data[(xoffset + k) - 1];
      }

      y->data[i] = s;
    }
  }
}

double c_combine_vector_elements(const emxArray_real_T *x)
{
  double y;
  int k;
  if (x->size[1] == 0) {
    y = 0.0;
  } else {
    y = x->data[0];
    for (k = 2; k <= x->size[1]; k++) {
      y += x->data[k - 1];
    }
  }

  return y;
}

void combine_vector_elements(const emxArray_real_T *x, emxArray_real_T *y)
{
  int vstride;
  int j;
  double s;
  int k;
  vstride = y->size[0];
  y->size[0] = x->size[0];
  emxEnsureCapacity((emxArray__common *)y, vstride, sizeof(double));
  if ((x->size[0] == 0) || (x->size[1] == 0)) {
    j = y->size[0];
    vstride = y->size[0];
    y->size[0] = j;
    emxEnsureCapacity((emxArray__common *)y, vstride, sizeof(double));
    for (vstride = 0; vstride < j; vstride++) {
      y->data[vstride] = 0.0;
    }
  } else {
    vstride = x->size[0];
    for (j = 0; j + 1 <= vstride; j++) {
      s = x->data[j];
      for (k = 2; k <= x->size[1]; k++) {
        s += x->data[j + (k - 1) * vstride];
      }

      y->data[j] = s;
    }
  }
}

/* End of code generation (combine_vector_elements.cpp) */
