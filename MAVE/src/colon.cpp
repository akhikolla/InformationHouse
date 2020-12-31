/*
 * Student License - for use by students to meet course requirements and
 * perform academic research at degree granting institutions only.  Not
 * for government, commercial, or other organizational use.
 *
 * colon.cpp
 *
 * Code generation for function 'colon'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "MAVEfast.h"
#include "colon.h"
#include "MAVEfast_emxutil.h"

/* Function Definitions */
void eml_signed_integer_colon(int b, emxArray_int32_T *y)
{
  int n;
  int yk;
  int k;
  if (b < 1) {
    n = 0;
  } else {
    n = b;
  }

  yk = y->size[0] * y->size[1];
  y->size[0] = 1;
  y->size[1] = n;
  emxEnsureCapacity((emxArray__common *)y, yk, sizeof(int));
  if (n > 0) {
    y->data[0] = 1;
    yk = 1;
    for (k = 2; k <= n; k++) {
      yk++;
      y->data[k - 1] = yk;
    }
  }
}

/* End of code generation (colon.cpp) */
