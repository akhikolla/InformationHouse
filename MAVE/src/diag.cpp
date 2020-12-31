/*
 * Student License - for use by students to meet course requirements and
 * perform academic research at degree granting institutions only.  Not
 * for government, commercial, or other organizational use.
 *
 * diag.cpp
 *
 * Code generation for function 'diag'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "MAVEfast.h"
#include "diag.h"
#include "MAVEfast_emxutil.h"

/* Function Definitions */
void b_diag(const emxArray_real_T *v, emxArray_real_T *d)
{
  int unnamed_idx_0;
  int unnamed_idx_1;
  int i7;
  unnamed_idx_0 = v->size[0];
  unnamed_idx_1 = v->size[0];
  i7 = d->size[0] * d->size[1];
  d->size[0] = unnamed_idx_0;
  d->size[1] = unnamed_idx_1;
  emxEnsureCapacity((emxArray__common *)d, i7, sizeof(double));
  unnamed_idx_0 *= unnamed_idx_1;
  for (i7 = 0; i7 < unnamed_idx_0; i7++) {
    d->data[i7] = 0.0;
  }

  for (unnamed_idx_0 = 0; unnamed_idx_0 + 1 <= v->size[0]; unnamed_idx_0++) {
    d->data[unnamed_idx_0 + d->size[0] * unnamed_idx_0] = v->data[unnamed_idx_0];
  }
}

void c_diag(const emxArray_creal_T *v, emxArray_creal_T *d)
{
  int u0;
  int u1;
  int stride;
  if ((v->size[0] == 1) && (v->size[1] == 1)) {
    u0 = d->size[0];
    d->size[0] = 1;
    emxEnsureCapacity((emxArray__common *)d, u0, sizeof(creal_T));
    d->data[0] = v->data[0];
  } else {
    if (0 < v->size[1]) {
      u0 = v->size[0];
      u1 = v->size[1];
      if (u0 < u1) {
        u1 = u0;
      }

      stride = v->size[0] + 1;
    } else {
      u1 = 0;
      stride = 0;
    }

    u0 = d->size[0];
    d->size[0] = u1;
    emxEnsureCapacity((emxArray__common *)d, u0, sizeof(creal_T));
    for (u0 = 0; u0 + 1 <= u1; u0++) {
      d->data[u0] = v->data[u0 * stride];
    }
  }
}

void diag(const emxArray_real_T *v, emxArray_real_T *d)
{
  int u0;
  int u1;
  int stride;
  if ((v->size[0] == 1) && (v->size[1] == 1)) {
    u0 = d->size[0];
    d->size[0] = 1;
    emxEnsureCapacity((emxArray__common *)d, u0, sizeof(double));
    d->data[0] = v->data[0];
  } else {
    if (0 < v->size[1]) {
      u0 = v->size[0];
      u1 = v->size[1];
      if (u0 < u1) {
        u1 = u0;
      }

      stride = v->size[0] + 1;
    } else {
      u1 = 0;
      stride = 0;
    }

    u0 = d->size[0];
    d->size[0] = u1;
    emxEnsureCapacity((emxArray__common *)d, u0, sizeof(double));
    for (u0 = 0; u0 + 1 <= u1; u0++) {
      d->data[u0] = v->data[u0 * stride];
    }
  }
}

/* End of code generation (diag.cpp) */
