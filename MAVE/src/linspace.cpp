/*
 * Student License - for use by students to meet course requirements and
 * perform academic research at degree granting institutions only.  Not
 * for government, commercial, or other organizational use.
 *
 * linspace.cpp
 *
 * Code generation for function 'linspace'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "MAVEfast.h"
#include "linspace.h"
#include "MAVEfast_emxutil.h"

/* Function Definitions */
void linspace(double d1, double d2, double n1, emxArray_real_T *y)
{
  int i10;
  double delta1;
  int k;
  i10 = y->size[0] * y->size[1];
  y->size[0] = 1;
  y->size[1] = (int)std::floor(n1);
  emxEnsureCapacity((emxArray__common *)y, i10, sizeof(double));
  if (y->size[1] >= 1) {
    y->data[y->size[1] - 1] = d2;
    if (y->size[1] >= 2) {
      y->data[0] = d1;
      if (y->size[1] >= 3) {
        delta1 = (d2 - d1) / ((double)y->size[1] - 1.0);
        i10 = y->size[1];
        for (k = 0; k <= i10 - 3; k++) {
          y->data[k + 1] = d1 + (1.0 + (double)k) * delta1;
        }
      }
    }
  }
}

/* End of code generation (linspace.cpp) */
