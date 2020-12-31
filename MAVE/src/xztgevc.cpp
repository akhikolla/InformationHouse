/*
 * Student License - for use by students to meet course requirements and
 * perform academic research at degree granting institutions only.  Not
 * for government, commercial, or other organizational use.
 *
 * xztgevc.cpp
 *
 * Code generation for function 'xztgevc'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "MAVEfast.h"
#include "xztgevc.h"
#include "MAVEfast_emxutil.h"

/* Function Definitions */
void xztgevc(const emxArray_creal_T *A, emxArray_creal_T *V)
{
  emxArray_creal_T *work1;
  int n;
  int i19;
  int loop_ub;
  emxArray_creal_T *work2;
  emxArray_real_T *rworka;
  double SMALL;
  double BIG;
  double BIGNUM;
  double anorm;
  int j;
  double xmx;
  double ascale;
  double y;
  int je;
  int b_je;
  double temp;
  double scale;
  double salpha_re;
  double salpha_im;
  double acoeff;
  boolean_T lscalea;
  boolean_T lscaleb;
  double acoefa;
  int jr;
  double dmin;
  int b_j;
  double d_re;
  double d_im;
  double work1_im;
  emxInit_creal_T1(&work1, 1);
  n = A->size[0] - 1;
  i19 = work1->size[0];
  work1->size[0] = A->size[0];
  emxEnsureCapacity((emxArray__common *)work1, i19, sizeof(creal_T));
  loop_ub = A->size[0];
  for (i19 = 0; i19 < loop_ub; i19++) {
    work1->data[i19].re = 0.0;
    work1->data[i19].im = 0.0;
  }

  emxInit_creal_T1(&work2, 1);
  i19 = work2->size[0];
  work2->size[0] = A->size[0];
  emxEnsureCapacity((emxArray__common *)work2, i19, sizeof(creal_T));
  loop_ub = A->size[0];
  for (i19 = 0; i19 < loop_ub; i19++) {
    work2->data[i19].re = 0.0;
    work2->data[i19].im = 0.0;
  }

  emxInit_real_T(&rworka, 1);
  SMALL = 2.2250738585072014E-308 * (double)A->size[0] / 2.2204460492503131E-16;
  BIG = 1.0 / SMALL;
  BIGNUM = 1.0 / (2.2250738585072014E-308 * (double)A->size[0]);
  i19 = rworka->size[0];
  rworka->size[0] = A->size[0];
  emxEnsureCapacity((emxArray__common *)rworka, i19, sizeof(double));
  loop_ub = A->size[0];
  for (i19 = 0; i19 < loop_ub; i19++) {
    rworka->data[i19] = 0.0;
  }

  anorm = std::abs(A->data[0].re) + std::abs(A->data[0].im);
  for (j = 1; j - 1 < n; j++) {
    for (loop_ub = 0; loop_ub < j; loop_ub++) {
      rworka->data[j] += std::abs(A->data[loop_ub + A->size[0] * j].re) + std::
        abs(A->data[loop_ub + A->size[0] * j].im);
    }

    y = rworka->data[j] + (std::abs(A->data[j + A->size[0] * j].re) + std::abs
      (A->data[j + A->size[0] * j].im));
    if (y > anorm) {
      anorm = y;
    }
  }

  xmx = anorm;
  if (2.2250738585072014E-308 > anorm) {
    xmx = 2.2250738585072014E-308;
  }

  ascale = 1.0 / xmx;
  i19 = (int)((1.0 + (-1.0 - (double)A->size[0])) / -1.0);
  for (je = 0; je < i19; je++) {
    b_je = n - je;
    xmx = (std::abs(A->data[b_je + A->size[0] * b_je].re) + std::abs(A->
            data[b_je + A->size[0] * b_je].im)) * ascale;
    if (1.0 > xmx) {
      xmx = 1.0;
    }

    temp = 1.0 / xmx;
    xmx = temp * A->data[b_je + A->size[0] * b_je].re;
    scale = temp * A->data[b_je + A->size[0] * b_je].im;
    salpha_re = ascale * xmx;
    salpha_im = ascale * scale;
    acoeff = temp * ascale;
    if ((std::abs(temp) >= 2.2250738585072014E-308) && (std::abs(acoeff) < SMALL))
    {
      lscalea = true;
    } else {
      lscalea = false;
    }

    if ((std::abs(salpha_re) + std::abs(salpha_im) >= 2.2250738585072014E-308) &&
        (std::abs(salpha_re) + std::abs(salpha_im) < SMALL)) {
      lscaleb = true;
    } else {
      lscaleb = false;
    }

    scale = 1.0;
    if (lscalea) {
      xmx = anorm;
      if (BIG < anorm) {
        xmx = BIG;
      }

      scale = SMALL / std::abs(temp) * xmx;
    }

    if (lscaleb) {
      xmx = 1.0;
      if (BIG < 1.0) {
        xmx = BIG;
      }

      y = SMALL / (std::abs(salpha_re) + std::abs(salpha_im)) * xmx;
      if (y > scale) {
        scale = y;
      }
    }

    if (lscalea || lscaleb) {
      xmx = std::abs(acoeff);
      y = std::abs(salpha_re) + std::abs(salpha_im);
      if (1.0 > xmx) {
        xmx = 1.0;
      }

      if (y > xmx) {
        xmx = y;
      }

      y = 1.0 / (2.2250738585072014E-308 * xmx);
      if (y < scale) {
        scale = y;
      }

      if (lscalea) {
        acoeff = ascale * (scale * temp);
      } else {
        acoeff *= scale;
      }

      salpha_re *= scale;
      salpha_im *= scale;
    }

    acoefa = std::abs(acoeff);
    for (jr = 0; jr <= n; jr++) {
      work1->data[jr].re = 0.0;
      work1->data[jr].im = 0.0;
    }

    work1->data[b_je].re = 1.0;
    work1->data[b_je].im = 0.0;
    dmin = 2.2204460492503131E-16 * acoefa * anorm;
    y = 2.2204460492503131E-16 * (std::abs(salpha_re) + std::abs(salpha_im));
    if (y > dmin) {
      dmin = y;
    }

    if (2.2250738585072014E-308 > dmin) {
      dmin = 2.2250738585072014E-308;
    }

    for (jr = 0; jr < b_je; jr++) {
      work1->data[jr].re = acoeff * A->data[jr + A->size[0] * b_je].re;
      work1->data[jr].im = acoeff * A->data[jr + A->size[0] * b_je].im;
    }

    work1->data[b_je].re = 1.0;
    work1->data[b_je].im = 0.0;
    loop_ub = (int)((1.0 + (-1.0 - ((double)(b_je + 1) - 1.0))) / -1.0);
    for (j = 0; j < loop_ub; j++) {
      b_j = (b_je - j) - 1;
      d_re = acoeff * A->data[b_j + A->size[0] * b_j].re - salpha_re;
      d_im = acoeff * A->data[b_j + A->size[0] * b_j].im - salpha_im;
      if (std::abs(d_re) + std::abs(d_im) <= dmin) {
        d_re = dmin;
        d_im = 0.0;
      }

      if ((std::abs(d_re) + std::abs(d_im) < 1.0) && (std::abs(work1->data[b_j].
            re) + std::abs(work1->data[b_j].im) >= BIGNUM * (std::abs(d_re) +
            std::abs(d_im)))) {
        temp = 1.0 / (std::abs(work1->data[b_j].re) + std::abs(work1->data[b_j].
          im));
        for (jr = 0; jr <= b_je; jr++) {
          work1->data[jr].re *= temp;
          work1->data[jr].im *= temp;
        }
      }

      temp = -work1->data[b_j].re;
      work1_im = -work1->data[b_j].im;
      if (d_im == 0.0) {
        if (work1_im == 0.0) {
          work1->data[b_j].re = temp / d_re;
          work1->data[b_j].im = 0.0;
        } else if (temp == 0.0) {
          work1->data[b_j].re = 0.0;
          work1->data[b_j].im = work1_im / d_re;
        } else {
          work1->data[b_j].re = temp / d_re;
          work1->data[b_j].im = work1_im / d_re;
        }
      } else if (d_re == 0.0) {
        if (temp == 0.0) {
          work1->data[b_j].re = work1_im / d_im;
          work1->data[b_j].im = 0.0;
        } else if (work1_im == 0.0) {
          work1->data[b_j].re = 0.0;
          work1->data[b_j].im = -(temp / d_im);
        } else {
          work1->data[b_j].re = work1_im / d_im;
          work1->data[b_j].im = -(temp / d_im);
        }
      } else {
        y = std::abs(d_re);
        xmx = std::abs(d_im);
        if (y > xmx) {
          scale = d_im / d_re;
          xmx = d_re + scale * d_im;
          work1->data[b_j].re = (temp + scale * work1_im) / xmx;
          work1->data[b_j].im = (work1_im - scale * temp) / xmx;
        } else if (xmx == y) {
          if (d_re > 0.0) {
            scale = 0.5;
          } else {
            scale = -0.5;
          }

          if (d_im > 0.0) {
            xmx = 0.5;
          } else {
            xmx = -0.5;
          }

          work1->data[b_j].re = (temp * scale + work1_im * xmx) / y;
          work1->data[b_j].im = (work1_im * scale - temp * xmx) / y;
        } else {
          scale = d_re / d_im;
          xmx = d_im + scale * d_re;
          work1->data[b_j].re = (scale * temp + work1_im) / xmx;
          work1->data[b_j].im = (scale * work1_im - temp) / xmx;
        }
      }

      if (b_j + 1 > 1) {
        if (std::abs(work1->data[b_j].re) + std::abs(work1->data[b_j].im) > 1.0)
        {
          temp = 1.0 / (std::abs(work1->data[b_j].re) + std::abs(work1->data[b_j]
            .im));
          if (acoefa * rworka->data[b_j] >= BIGNUM * temp) {
            for (jr = 0; jr <= b_je; jr++) {
              work1->data[jr].re *= temp;
              work1->data[jr].im *= temp;
            }
          }
        }

        d_re = acoeff * work1->data[b_j].re;
        d_im = acoeff * work1->data[b_j].im;
        for (jr = 0; jr < b_j; jr++) {
          xmx = d_re * A->data[jr + A->size[0] * b_j].re - d_im * A->data[jr +
            A->size[0] * b_j].im;
          scale = d_re * A->data[jr + A->size[0] * b_j].im + d_im * A->data[jr +
            A->size[0] * b_j].re;
          work1->data[jr].re += xmx;
          work1->data[jr].im += scale;
        }
      }
    }

    for (jr = 0; jr <= n; jr++) {
      work2->data[jr].re = 0.0;
      work2->data[jr].im = 0.0;
    }

    for (loop_ub = 0; loop_ub <= b_je; loop_ub++) {
      for (jr = 0; jr <= n; jr++) {
        xmx = V->data[jr + V->size[0] * loop_ub].re * work1->data[loop_ub].re -
          V->data[jr + V->size[0] * loop_ub].im * work1->data[loop_ub].im;
        scale = V->data[jr + V->size[0] * loop_ub].re * work1->data[loop_ub].im
          + V->data[jr + V->size[0] * loop_ub].im * work1->data[loop_ub].re;
        work2->data[jr].re += xmx;
        work2->data[jr].im += scale;
      }
    }

    xmx = std::abs(work2->data[0].re) + std::abs(work2->data[0].im);
    if (n + 1 > 1) {
      for (jr = 1; jr - 1 < n; jr++) {
        y = std::abs(work2->data[jr].re) + std::abs(work2->data[jr].im);
        if (y > xmx) {
          xmx = y;
        }
      }
    }

    if (xmx > 2.2250738585072014E-308) {
      temp = 1.0 / xmx;
      for (jr = 0; jr <= n; jr++) {
        V->data[jr + V->size[0] * b_je].re = temp * work2->data[jr].re;
        V->data[jr + V->size[0] * b_je].im = temp * work2->data[jr].im;
      }
    } else {
      for (jr = 0; jr <= n; jr++) {
        V->data[jr + V->size[0] * b_je].re = 0.0;
        V->data[jr + V->size[0] * b_je].im = 0.0;
      }
    }
  }

  emxFree_real_T(&rworka);
  emxFree_creal_T(&work2);
  emxFree_creal_T(&work1);
}

/* End of code generation (xztgevc.cpp) */
