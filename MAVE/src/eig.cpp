/*
 * Student License - for use by students to meet course requirements and
 * perform academic research at degree granting institutions only.  Not
 * for government, commercial, or other organizational use.
 *
 * eig.cpp
 *
 * Code generation for function 'eig'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "MAVEfast.h"
#include "eig.h"
#include "MAVEfast_emxutil.h"
#include "schur.h"
#include "xzlarfg.h"
#include "xzggev.h"
#include "MAVEfast_rtwutil.h"

/* Function Definitions */
void eig(const emxArray_real_T *A, emxArray_creal_T *V, emxArray_creal_T *D)
{
  int nx;
  int i5;
  boolean_T p;
  int k;
  int j;
  unsigned int uv2[2];
  boolean_T exitg2;
  emxArray_creal_T *At;
  int n;
  double absxk;
  int exitg1;
  emxArray_creal_T *alpha1;
  emxArray_creal_T *beta1;
  int c;
  int coltop;
  double colnorm;
  double scale;
  double t;
  double alpha1_re;
  double alpha1_im;
  if ((A->size[0] == 0) || (A->size[1] == 0)) {
    i5 = V->size[0] * V->size[1];
    V->size[0] = A->size[0];
    V->size[1] = A->size[1];
    emxEnsureCapacity((emxArray__common *)V, i5, sizeof(creal_T));
    nx = A->size[0] * A->size[1];
    for (i5 = 0; i5 < nx; i5++) {
      V->data[i5].re = A->data[i5];
      V->data[i5].im = 0.0;
    }

    i5 = D->size[0] * D->size[1];
    D->size[0] = A->size[0];
    D->size[1] = A->size[1];
    emxEnsureCapacity((emxArray__common *)D, i5, sizeof(creal_T));
    nx = A->size[0] * A->size[1];
    for (i5 = 0; i5 < nx; i5++) {
      D->data[i5].re = A->data[i5];
      D->data[i5].im = 0.0;
    }
  } else {
    nx = A->size[0] * A->size[1];
    p = false;
    for (k = 0; k + 1 <= nx; k++) {
      if (p || rtIsInf(A->data[k]) || rtIsNaN(A->data[k])) {
        p = true;
      } else {
        p = false;
      }
    }

    if (p) {
      if ((A->size[0] == 1) && (A->size[1] == 1)) {
        i5 = V->size[0] * V->size[1];
        V->size[0] = 1;
        V->size[1] = 1;
        emxEnsureCapacity((emxArray__common *)V, i5, sizeof(creal_T));
        for (i5 = 0; i5 < 1; i5++) {
          V->data[0].re = 1.0;
          V->data[0].im = 0.0;
        }

        i5 = D->size[0] * D->size[1];
        D->size[0] = A->size[0];
        D->size[1] = A->size[1];
        emxEnsureCapacity((emxArray__common *)D, i5, sizeof(creal_T));
        nx = A->size[0] * A->size[1];
        for (i5 = 0; i5 < nx; i5++) {
          D->data[i5].re = A->data[i5];
          D->data[i5].im = 0.0;
        }

        i5 = V->size[0] * V->size[1];
        emxEnsureCapacity((emxArray__common *)V, i5, sizeof(creal_T));
        nx = V->size[1];
        for (i5 = 0; i5 < nx; i5++) {
          k = V->size[0];
          for (j = 0; j < k; j++) {
            V->data[j + V->size[0] * i5].re = rtNaN;
            V->data[j + V->size[0] * i5].im = 0.0;
          }
        }

        i5 = D->size[0] * D->size[1];
        emxEnsureCapacity((emxArray__common *)D, i5, sizeof(creal_T));
        nx = D->size[1];
        for (i5 = 0; i5 < nx; i5++) {
          k = D->size[0];
          for (j = 0; j < k; j++) {
            D->data[j + D->size[0] * i5].re = rtNaN;
            D->data[j + D->size[0] * i5].im = 0.0;
          }
        }
      } else {
        for (i5 = 0; i5 < 2; i5++) {
          uv2[i5] = (unsigned int)A->size[i5];
        }

        i5 = V->size[0] * V->size[1];
        V->size[0] = (int)uv2[0];
        V->size[1] = (int)uv2[1];
        emxEnsureCapacity((emxArray__common *)V, i5, sizeof(creal_T));
        nx = (int)uv2[0] * (int)uv2[1];
        for (i5 = 0; i5 < nx; i5++) {
          V->data[i5].re = rtNaN;
          V->data[i5].im = 0.0;
        }

        for (i5 = 0; i5 < 2; i5++) {
          uv2[i5] = (unsigned int)A->size[i5];
        }

        i5 = D->size[0] * D->size[1];
        D->size[0] = (int)uv2[0];
        D->size[1] = (int)uv2[1];
        emxEnsureCapacity((emxArray__common *)D, i5, sizeof(creal_T));
        nx = (int)uv2[0] * (int)uv2[1];
        for (i5 = 0; i5 < nx; i5++) {
          D->data[i5].re = 0.0;
          D->data[i5].im = 0.0;
        }

        for (k = 0; k + 1 <= (int)uv2[0]; k++) {
          D->data[k + D->size[0] * k].re = rtNaN;
          D->data[k + D->size[0] * k].im = 0.0;
        }
      }
    } else if ((A->size[0] == 1) && (A->size[1] == 1)) {
      i5 = V->size[0] * V->size[1];
      V->size[0] = 1;
      V->size[1] = 1;
      emxEnsureCapacity((emxArray__common *)V, i5, sizeof(creal_T));
      for (i5 = 0; i5 < 1; i5++) {
        V->data[0].re = 1.0;
        V->data[0].im = 0.0;
      }

      i5 = D->size[0] * D->size[1];
      D->size[0] = A->size[0];
      D->size[1] = A->size[1];
      emxEnsureCapacity((emxArray__common *)D, i5, sizeof(creal_T));
      nx = A->size[0] * A->size[1];
      for (i5 = 0; i5 < nx; i5++) {
        D->data[i5].re = A->data[i5];
        D->data[i5].im = 0.0;
      }
    } else {
      p = (A->size[0] == A->size[1]);
      if (p) {
        j = 0;
        exitg2 = false;
        while ((!exitg2) && (j <= A->size[1] - 1)) {
          nx = 0;
          do {
            exitg1 = 0;
            if (nx <= j) {
              if (!(A->data[nx + A->size[0] * j] == A->data[j + A->size[0] * nx]))
              {
                p = false;
                exitg1 = 1;
              } else {
                nx++;
              }
            } else {
              j++;
              exitg1 = 2;
            }
          } while (exitg1 == 0);

          if (exitg1 == 1) {
            exitg2 = true;
          }
        }
      }

      if (p) {
        schur(A, V, D);
        n = D->size[0];
        absxk = D->data[0].re;
        D->data[0].re = absxk;
        D->data[0].im = 0.0;
        for (j = 1; j + 1 <= n; j++) {
          absxk = D->data[j + D->size[0] * j].re;
          D->data[j + D->size[0] * j].re = absxk;
          D->data[j + D->size[0] * j].im = 0.0;
          D->data[j + D->size[0] * (j - 1)].re = 0.0;
          D->data[j + D->size[0] * (j - 1)].im = 0.0;
          for (nx = 1; nx <= j; nx++) {
            D->data[(nx + D->size[0] * j) - 1].re = 0.0;
            D->data[(nx + D->size[0] * j) - 1].im = 0.0;
          }
        }
      } else {
        emxInit_creal_T(&At, 2);
        i5 = At->size[0] * At->size[1];
        At->size[0] = A->size[0];
        At->size[1] = A->size[1];
        emxEnsureCapacity((emxArray__common *)At, i5, sizeof(creal_T));
        nx = A->size[0] * A->size[1];
        for (i5 = 0; i5 < nx; i5++) {
          At->data[i5].re = A->data[i5];
          At->data[i5].im = 0.0;
        }

        emxInit_creal_T1(&alpha1, 1);
        emxInit_creal_T1(&beta1, 1);
        xzggev(At, &nx, alpha1, beta1, V);
        n = A->size[0];
        c = (A->size[0] - 1) * A->size[0];
        coltop = 0;
        emxFree_creal_T(&At);
        while (coltop + 1 <= c + 1) {
          colnorm = 0.0;
          if (n == 1) {
            colnorm = rt_hypotd_snf(V->data[coltop].re, V->data[coltop].im);
          } else {
            scale = 2.2250738585072014E-308;
            nx = coltop + n;
            for (k = coltop; k + 1 <= nx; k++) {
              absxk = std::abs(V->data[k].re);
              if (absxk > scale) {
                t = scale / absxk;
                colnorm = 1.0 + colnorm * t * t;
                scale = absxk;
              } else {
                t = absxk / scale;
                colnorm += t * t;
              }

              absxk = std::abs(V->data[k].im);
              if (absxk > scale) {
                t = scale / absxk;
                colnorm = 1.0 + colnorm * t * t;
                scale = absxk;
              } else {
                t = absxk / scale;
                colnorm += t * t;
              }
            }

            colnorm = scale * std::sqrt(colnorm);
          }

          i5 = coltop + n;
          for (j = coltop; j + 1 <= i5; j++) {
            absxk = V->data[j].re;
            scale = V->data[j].im;
            if (scale == 0.0) {
              V->data[j].re = absxk / colnorm;
              V->data[j].im = 0.0;
            } else if (absxk == 0.0) {
              V->data[j].re = 0.0;
              V->data[j].im = scale / colnorm;
            } else {
              V->data[j].re = absxk / colnorm;
              V->data[j].im = scale / colnorm;
            }
          }

          coltop += n;
        }

        i5 = D->size[0] * D->size[1];
        D->size[0] = alpha1->size[0];
        D->size[1] = alpha1->size[0];
        emxEnsureCapacity((emxArray__common *)D, i5, sizeof(creal_T));
        nx = alpha1->size[0] * alpha1->size[0];
        for (i5 = 0; i5 < nx; i5++) {
          D->data[i5].re = 0.0;
          D->data[i5].im = 0.0;
        }

        for (k = 0; k < alpha1->size[0]; k++) {
          alpha1_re = alpha1->data[k].re;
          alpha1_im = alpha1->data[k].im;
          absxk = beta1->data[k].re;
          t = beta1->data[k].im;
          if (t == 0.0) {
            if (alpha1_im == 0.0) {
              D->data[k + D->size[0] * k].re = alpha1_re / absxk;
              D->data[k + D->size[0] * k].im = 0.0;
            } else if (alpha1_re == 0.0) {
              D->data[k + D->size[0] * k].re = 0.0;
              D->data[k + D->size[0] * k].im = alpha1_im / absxk;
            } else {
              D->data[k + D->size[0] * k].re = alpha1_re / absxk;
              D->data[k + D->size[0] * k].im = alpha1_im / absxk;
            }
          } else if (absxk == 0.0) {
            if (alpha1_re == 0.0) {
              D->data[k + D->size[0] * k].re = alpha1_im / t;
              D->data[k + D->size[0] * k].im = 0.0;
            } else if (alpha1_im == 0.0) {
              D->data[k + D->size[0] * k].re = 0.0;
              D->data[k + D->size[0] * k].im = -(alpha1_re / t);
            } else {
              D->data[k + D->size[0] * k].re = alpha1_im / t;
              D->data[k + D->size[0] * k].im = -(alpha1_re / t);
            }
          } else {
            colnorm = std::abs(absxk);
            scale = std::abs(t);
            if (colnorm > scale) {
              scale = t / absxk;
              absxk += scale * t;
              D->data[k + D->size[0] * k].re = (alpha1_re + scale * alpha1_im) /
                absxk;
              D->data[k + D->size[0] * k].im = (alpha1_im - scale * alpha1_re) /
                absxk;
            } else if (scale == colnorm) {
              if (absxk > 0.0) {
                absxk = 0.5;
              } else {
                absxk = -0.5;
              }

              if (t > 0.0) {
                scale = 0.5;
              } else {
                scale = -0.5;
              }

              D->data[k + D->size[0] * k].re = (alpha1_re * absxk + alpha1_im *
                scale) / colnorm;
              D->data[k + D->size[0] * k].im = (alpha1_im * absxk - alpha1_re *
                scale) / colnorm;
            } else {
              scale = absxk / t;
              absxk = t + scale * absxk;
              D->data[k + D->size[0] * k].re = (scale * alpha1_re + alpha1_im) /
                absxk;
              D->data[k + D->size[0] * k].im = (scale * alpha1_im - alpha1_re) /
                absxk;
            }
          }
        }

        emxFree_creal_T(&beta1);
        emxFree_creal_T(&alpha1);
      }
    }
  }
}

/* End of code generation (eig.cpp) */
