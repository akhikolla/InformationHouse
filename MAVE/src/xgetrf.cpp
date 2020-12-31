/*
 * Student License - for use by students to meet course requirements and
 * perform academic research at degree granting institutions only.  Not
 * for government, commercial, or other organizational use.
 *
 * xgetrf.cpp
 *
 * Code generation for function 'xgetrf'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "MAVEfast.h"
#include "xgetrf.h"
#include "colon.h"

/* Function Definitions */
void xgetrf(int m, int n, emxArray_real_T *A, int lda, emxArray_int32_T *ipiv,
            int *info)
{
  int iy;
  int b_info;
  int u0;
  int j;
  int mmj;
  int c;
  int ix;
  double smax;
  int jA;
  int i22;
  int jy;
  double s;
  int b_j;
  int ijA;
  if (m < n) {
    iy = m;
  } else {
    iy = n;
  }

  eml_signed_integer_colon(iy, ipiv);
  b_info = 0;
  if ((m < 1) || (n < 1)) {
  } else {
    u0 = m - 1;
    if (!(u0 < n)) {
      u0 = n;
    }

    for (j = 0; j + 1 <= u0; j++) {
      mmj = m - j;
      c = j * (lda + 1);
      if (mmj < 1) {
        iy = -1;
      } else {
        iy = 0;
        if (mmj > 1) {
          ix = c;
          smax = std::abs(A->data[c]);
          for (jA = 1; jA + 1 <= mmj; jA++) {
            ix++;
            s = std::abs(A->data[ix]);
            if (s > smax) {
              iy = jA;
              smax = s;
            }
          }
        }
      }

      if (A->data[c + iy] != 0.0) {
        if (iy != 0) {
          ipiv->data[j] = (j + iy) + 1;
          ix = j;
          iy += j;
          for (jA = 1; jA <= n; jA++) {
            smax = A->data[ix];
            A->data[ix] = A->data[iy];
            A->data[iy] = smax;
            ix += lda;
            iy += lda;
          }
        }

        i22 = c + mmj;
        for (iy = c + 1; iy + 1 <= i22; iy++) {
          A->data[iy] /= A->data[c];
        }
      } else {
        b_info = j + 1;
      }

      iy = (n - j) - 1;
      jA = (c + lda) + 1;
      jy = c + lda;
      for (b_j = 1; b_j <= iy; b_j++) {
        smax = A->data[jy];
        if (A->data[jy] != 0.0) {
          ix = c + 1;
          i22 = mmj + jA;
          for (ijA = jA; ijA + 1 < i22; ijA++) {
            A->data[ijA] += A->data[ix] * -smax;
            ix++;
          }
        }

        jy += lda;
        jA += lda;
      }
    }

    if ((b_info == 0) && (m <= n) && (!(A->data[(m + A->size[0] * (m - 1)) - 1]
          != 0.0))) {
      b_info = m;
    }
  }

  *info = b_info;
}

/* End of code generation (xgetrf.cpp) */
