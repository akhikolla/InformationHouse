/*
 * Student License - for use by students to meet course requirements and
 * perform academic research at degree granting institutions only.  Not
 * for government, commercial, or other organizational use.
 *
 * schur.cpp
 *
 * Code generation for function 'schur'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "MAVEfast.h"
#include "schur.h"
#include "MAVEfast_emxutil.h"
#include "xscal.h"
#include "xzlarf.h"
#include "xzlarfg.h"
#include "xdlanv2.h"
#include "xdhseqr.h"
#include "xgehrd.h"
#include "MAVEfast_rtwutil.h"

/* Function Definitions */
void schur(const emxArray_real_T *A, emxArray_creal_T *V, emxArray_creal_T *T)
{
  int nx;
  boolean_T p;
  int iajm1;
  emxArray_real_T *b_A;
  int n;
  int varargin_1[2];
  emxArray_real_T *Vr;
  int itau;
  emxArray_real_T *tau;
  int j;
  int i;
  int nh;
  emxArray_real_T *work;
  double r;
  double b;
  double c;
  double d;
  double s;
  double rt1i;
  double t1_re;
  double t1_im;
  double mu1_im;
  double mu1_re;
  nx = A->size[0] * A->size[1];
  p = false;
  for (iajm1 = 0; iajm1 + 1 <= nx; iajm1++) {
    if (p || rtIsInf(A->data[iajm1]) || rtIsNaN(A->data[iajm1])) {
      p = true;
    } else {
      p = false;
    }
  }

  if (p) {
    for (iajm1 = 0; iajm1 < 2; iajm1++) {
      varargin_1[iajm1] = A->size[iajm1];
    }

    iajm1 = V->size[0] * V->size[1];
    V->size[0] = varargin_1[0];
    V->size[1] = varargin_1[1];
    emxEnsureCapacity((emxArray__common *)V, iajm1, sizeof(creal_T));
    nx = varargin_1[0] * varargin_1[1];
    for (iajm1 = 0; iajm1 < nx; iajm1++) {
      V->data[iajm1].re = rtNaN;
      V->data[iajm1].im = 0.0;
    }

    itau = V->size[0];
    if ((V->size[0] == 0) || (V->size[1] == 0) || (2 >= V->size[0])) {
    } else {
      nx = 3;
      if (V->size[0] - 3 < V->size[1] - 1) {
        iajm1 = V->size[0] - 2;
      } else {
        iajm1 = V->size[1];
      }

      for (j = 1; j <= iajm1; j++) {
        for (i = nx; i <= itau; i++) {
          V->data[(i + V->size[0] * (j - 1)) - 1].re = 0.0;
          V->data[(i + V->size[0] * (j - 1)) - 1].im = 0.0;
        }

        nx++;
      }
    }

    for (iajm1 = 0; iajm1 < 2; iajm1++) {
      varargin_1[iajm1] = A->size[iajm1];
    }

    iajm1 = T->size[0] * T->size[1];
    T->size[0] = varargin_1[0];
    T->size[1] = varargin_1[1];
    emxEnsureCapacity((emxArray__common *)T, iajm1, sizeof(creal_T));
    nx = varargin_1[0] * varargin_1[1];
    for (iajm1 = 0; iajm1 < nx; iajm1++) {
      T->data[iajm1].re = rtNaN;
      T->data[iajm1].im = 0.0;
    }
  } else {
    emxInit_real_T1(&b_A, 2);
    n = A->size[0];
    iajm1 = b_A->size[0] * b_A->size[1];
    b_A->size[0] = A->size[0];
    b_A->size[1] = A->size[1];
    emxEnsureCapacity((emxArray__common *)b_A, iajm1, sizeof(double));
    nx = A->size[0] * A->size[1];
    for (iajm1 = 0; iajm1 < nx; iajm1++) {
      b_A->data[iajm1] = A->data[iajm1];
    }

    emxInit_real_T1(&Vr, 2);
    emxInit_real_T(&tau, 1);
    xgehrd(b_A, tau);
    iajm1 = Vr->size[0] * Vr->size[1];
    Vr->size[0] = b_A->size[0];
    Vr->size[1] = b_A->size[1];
    emxEnsureCapacity((emxArray__common *)Vr, iajm1, sizeof(double));
    nx = b_A->size[0] * b_A->size[1];
    for (iajm1 = 0; iajm1 < nx; iajm1++) {
      Vr->data[iajm1] = b_A->data[iajm1];
    }

    if (A->size[0] != 0) {
      nh = A->size[0] - 1;
      for (j = A->size[0]; j > 1; j--) {
        nx = (j - 1) * n;
        for (i = 1; i < j; i++) {
          Vr->data[(nx + i) - 1] = 0.0;
        }

        iajm1 = nx - n;
        for (i = j; i + 1 <= n; i++) {
          Vr->data[nx + i] = Vr->data[iajm1 + i];
        }

        for (i = n; i + 1 <= n; i++) {
          Vr->data[nx + i] = 0.0;
        }
      }

      for (i = 1; i <= n; i++) {
        Vr->data[i - 1] = 0.0;
      }

      Vr->data[0] = 1.0;
      for (j = A->size[0]; j + 1 <= n; j++) {
        nx = j * n;
        for (i = 1; i <= n; i++) {
          Vr->data[(nx + i) - 1] = 0.0;
        }

        Vr->data[nx + j] = 1.0;
      }

      if (!(A->size[0] - 1 < 1)) {
        for (j = A->size[0]; j - 1 < nh; j++) {
          nx = n + (j - 1) * n;
          for (i = 0; i < nh; i++) {
            Vr->data[(nx + i) + 1] = 0.0;
          }

          Vr->data[nx + j] = 1.0;
        }

        emxInit_real_T(&work, 1);
        itau = A->size[0] - 2;
        varargin_1[0] = Vr->size[1];
        iajm1 = work->size[0];
        work->size[0] = varargin_1[0];
        emxEnsureCapacity((emxArray__common *)work, iajm1, sizeof(double));
        nx = varargin_1[0];
        for (iajm1 = 0; iajm1 < nx; iajm1++) {
          work->data[iajm1] = 0.0;
        }

        for (i = A->size[0] - 1; i >= 1; i--) {
          nx = ((A->size[0] + i) + (i - 1) * n) + 1;
          if (i < n - 1) {
            Vr->data[nx - 1] = 1.0;
            xzlarf(n - i, nh - i, nx, tau->data[itau], Vr, nx + n, n, work);
            xscal(nh - i, -tau->data[itau], Vr, nx + 1);
          }

          Vr->data[nx - 1] = 1.0 - tau->data[itau];
          for (j = 1; j < i; j++) {
            Vr->data[(nx - j) - 1] = 0.0;
          }

          itau--;
        }

        emxFree_real_T(&work);
      }
    }

    emxFree_real_T(&tau);
    eml_dlahqr(b_A, Vr);
    itau = b_A->size[0];
    if ((b_A->size[0] == 0) || (b_A->size[1] == 0) || (3 >= b_A->size[0])) {
    } else {
      nx = 4;
      if (b_A->size[0] - 4 < b_A->size[1] - 1) {
        iajm1 = b_A->size[0] - 3;
      } else {
        iajm1 = b_A->size[1];
      }

      for (j = 1; j <= iajm1; j++) {
        for (i = nx; i <= itau; i++) {
          b_A->data[(i + b_A->size[0] * (j - 1)) - 1] = 0.0;
        }

        nx++;
      }
    }

    iajm1 = T->size[0] * T->size[1];
    T->size[0] = b_A->size[0];
    T->size[1] = b_A->size[1];
    emxEnsureCapacity((emxArray__common *)T, iajm1, sizeof(creal_T));
    nx = b_A->size[0] * b_A->size[1];
    for (iajm1 = 0; iajm1 < nx; iajm1++) {
      T->data[iajm1].re = b_A->data[iajm1];
      T->data[iajm1].im = 0.0;
    }

    iajm1 = V->size[0] * V->size[1];
    V->size[0] = Vr->size[0];
    V->size[1] = Vr->size[1];
    emxEnsureCapacity((emxArray__common *)V, iajm1, sizeof(creal_T));
    nx = Vr->size[0] * Vr->size[1];
    for (iajm1 = 0; iajm1 < nx; iajm1++) {
      V->data[iajm1].re = Vr->data[iajm1];
      V->data[iajm1].im = 0.0;
    }

    for (iajm1 = 0; iajm1 < 2; iajm1++) {
      varargin_1[iajm1] = b_A->size[iajm1];
    }

    nx = varargin_1[0];
    if (varargin_1[1] < varargin_1[0]) {
      nx = varargin_1[1];
    }

    for (iajm1 = 0; iajm1 < 2; iajm1++) {
      varargin_1[iajm1] = Vr->size[iajm1];
    }

    emxFree_real_T(&Vr);
    iajm1 = varargin_1[0];
    if (varargin_1[1] < varargin_1[0]) {
      iajm1 = varargin_1[1];
    }

    if (nx < iajm1) {
      iajm1 = nx;
    }

    if (iajm1 != 0) {
      for (itau = iajm1 - 1; itau + 1 >= 2; itau--) {
        if (b_A->data[itau + b_A->size[0] * (itau - 1)] != 0.0) {
          r = b_A->data[(itau + b_A->size[0] * (itau - 1)) - 1];
          b = b_A->data[(itau + b_A->size[0] * itau) - 1];
          c = b_A->data[itau + b_A->size[0] * (itau - 1)];
          d = b_A->data[itau + b_A->size[0] * itau];
          xdlanv2(&r, &b, &c, &d, &s, &rt1i, &t1_re, &t1_im, &mu1_im, &mu1_re);
          mu1_re = s - b_A->data[itau + b_A->size[0] * itau];
          r = rt_hypotd_snf(rt_hypotd_snf(mu1_re, rt1i), b_A->data[itau +
                            b_A->size[0] * (itau - 1)]);
          if (rt1i == 0.0) {
            mu1_re /= r;
            mu1_im = 0.0;
          } else if (mu1_re == 0.0) {
            mu1_re = 0.0;
            mu1_im = rt1i / r;
          } else {
            mu1_re /= r;
            mu1_im = rt1i / r;
          }

          s = b_A->data[itau + b_A->size[0] * (itau - 1)] / r;
          for (j = itau - 1; j + 1 <= iajm1; j++) {
            t1_re = T->data[(itau + T->size[0] * j) - 1].re;
            t1_im = T->data[(itau + T->size[0] * j) - 1].im;
            c = T->data[(itau + T->size[0] * j) - 1].re;
            d = T->data[(itau + T->size[0] * j) - 1].im;
            r = T->data[(itau + T->size[0] * j) - 1].im;
            b = T->data[(itau + T->size[0] * j) - 1].re;
            T->data[(itau + T->size[0] * j) - 1].re = (mu1_re * c + mu1_im * d)
              + s * T->data[itau + T->size[0] * j].re;
            T->data[(itau + T->size[0] * j) - 1].im = (mu1_re * r - mu1_im * b)
              + s * T->data[itau + T->size[0] * j].im;
            r = mu1_re * T->data[itau + T->size[0] * j].re - mu1_im * T->
              data[itau + T->size[0] * j].im;
            b = mu1_re * T->data[itau + T->size[0] * j].im + mu1_im * T->
              data[itau + T->size[0] * j].re;
            T->data[itau + T->size[0] * j].re = r - s * t1_re;
            T->data[itau + T->size[0] * j].im = b - s * t1_im;
          }

          for (i = 0; i + 1 <= itau + 1; i++) {
            t1_re = T->data[i + T->size[0] * (itau - 1)].re;
            t1_im = T->data[i + T->size[0] * (itau - 1)].im;
            r = mu1_re * T->data[i + T->size[0] * (itau - 1)].re - mu1_im *
              T->data[i + T->size[0] * (itau - 1)].im;
            b = mu1_re * T->data[i + T->size[0] * (itau - 1)].im + mu1_im *
              T->data[i + T->size[0] * (itau - 1)].re;
            c = T->data[i + T->size[0] * itau].re;
            d = T->data[i + T->size[0] * itau].im;
            T->data[i + T->size[0] * (itau - 1)].re = r + s * c;
            T->data[i + T->size[0] * (itau - 1)].im = b + s * d;
            c = T->data[i + T->size[0] * itau].re;
            d = T->data[i + T->size[0] * itau].im;
            r = T->data[i + T->size[0] * itau].im;
            b = T->data[i + T->size[0] * itau].re;
            T->data[i + T->size[0] * itau].re = (mu1_re * c + mu1_im * d) - s *
              t1_re;
            T->data[i + T->size[0] * itau].im = (mu1_re * r - mu1_im * b) - s *
              t1_im;
          }

          for (i = 0; i + 1 <= iajm1; i++) {
            t1_re = V->data[i + V->size[0] * (itau - 1)].re;
            t1_im = V->data[i + V->size[0] * (itau - 1)].im;
            r = mu1_re * V->data[i + V->size[0] * (itau - 1)].re - mu1_im *
              V->data[i + V->size[0] * (itau - 1)].im;
            b = mu1_re * V->data[i + V->size[0] * (itau - 1)].im + mu1_im *
              V->data[i + V->size[0] * (itau - 1)].re;
            c = V->data[i + V->size[0] * itau].re;
            d = V->data[i + V->size[0] * itau].im;
            V->data[i + V->size[0] * (itau - 1)].re = r + s * c;
            V->data[i + V->size[0] * (itau - 1)].im = b + s * d;
            c = V->data[i + V->size[0] * itau].re;
            d = V->data[i + V->size[0] * itau].im;
            r = V->data[i + V->size[0] * itau].im;
            b = V->data[i + V->size[0] * itau].re;
            V->data[i + V->size[0] * itau].re = (mu1_re * c + mu1_im * d) - s *
              t1_re;
            V->data[i + V->size[0] * itau].im = (mu1_re * r - mu1_im * b) - s *
              t1_im;
          }

          T->data[itau + T->size[0] * (itau - 1)].re = 0.0;
          T->data[itau + T->size[0] * (itau - 1)].im = 0.0;
        }
      }
    }

    emxFree_real_T(&b_A);
  }
}

/* End of code generation (schur.cpp) */
