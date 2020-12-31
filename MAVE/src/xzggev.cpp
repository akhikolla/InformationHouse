/*
 * Student License - for use by students to meet course requirements and
 * perform academic research at degree granting institutions only.  Not
 * for government, commercial, or other organizational use.
 *
 * xzggev.cpp
 *
 * Code generation for function 'xzggev'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "MAVEfast.h"
#include "xzggev.h"
#include "MAVEfast_emxutil.h"
#include "xzlartg.h"
#include "xztgevc.h"
#include "xzhgeqz.h"
#include "xzggbal.h"
#include "xzlarfg.h"
#include "MAVEfast_rtwutil.h"

/* Function Definitions */
void xzggev(emxArray_creal_T *A, int *info, emxArray_creal_T *alpha1,
            emxArray_creal_T *beta1, emxArray_creal_T *V)
{
  int b_info;
  int n;
  int j;
  int loop_ub;
  double anrm;
  boolean_T ilascl;
  boolean_T notdone;
  int jcol;
  boolean_T exitg1;
  double anrmto;
  double absxk;
  emxArray_int32_T *rscale;
  double ctoc;
  emxArray_int8_T *I;
  int ilo;
  int ihi;
  int b_n;
  double cfrom1;
  double cto1;
  double stemp_im;
  int jrow;
  creal_T b_A;
  creal_T c_A;
  double c;
  creal_T tmp;
  double stemp_re;
  b_info = 0;
  n = A->size[0];
  j = alpha1->size[0];
  alpha1->size[0] = A->size[0];
  emxEnsureCapacity((emxArray__common *)alpha1, j, sizeof(creal_T));
  loop_ub = A->size[0];
  for (j = 0; j < loop_ub; j++) {
    alpha1->data[j].re = 0.0;
    alpha1->data[j].im = 0.0;
  }

  j = beta1->size[0];
  beta1->size[0] = A->size[0];
  emxEnsureCapacity((emxArray__common *)beta1, j, sizeof(creal_T));
  loop_ub = A->size[0];
  for (j = 0; j < loop_ub; j++) {
    beta1->data[j].re = 0.0;
    beta1->data[j].im = 0.0;
  }

  j = V->size[0] * V->size[1];
  V->size[0] = A->size[0];
  V->size[1] = A->size[0];
  emxEnsureCapacity((emxArray__common *)V, j, sizeof(creal_T));
  loop_ub = A->size[0] * A->size[0];
  for (j = 0; j < loop_ub; j++) {
    V->data[j].re = 0.0;
    V->data[j].im = 0.0;
  }

  if (!((A->size[0] == 0) || (A->size[1] == 0))) {
    anrm = 0.0;
    ilascl = (A->size[0] == 0);
    notdone = (A->size[1] == 0);
    if (!(ilascl || notdone)) {
      jcol = 0;
      exitg1 = false;
      while ((!exitg1) && (jcol <= A->size[0] * A->size[1] - 1)) {
        absxk = rt_hypotd_snf(A->data[jcol].re, A->data[jcol].im);
        if (rtIsNaN(absxk)) {
          anrm = rtNaN;
          exitg1 = true;
        } else {
          if (absxk > anrm) {
            anrm = absxk;
          }

          jcol++;
        }
      }
    }

    if (!((!rtIsInf(anrm)) && (!rtIsNaN(anrm)))) {
      j = alpha1->size[0];
      alpha1->size[0] = A->size[0];
      emxEnsureCapacity((emxArray__common *)alpha1, j, sizeof(creal_T));
      loop_ub = A->size[0];
      for (j = 0; j < loop_ub; j++) {
        alpha1->data[j].re = rtNaN;
        alpha1->data[j].im = 0.0;
      }

      j = beta1->size[0];
      beta1->size[0] = A->size[0];
      emxEnsureCapacity((emxArray__common *)beta1, j, sizeof(creal_T));
      loop_ub = A->size[0];
      for (j = 0; j < loop_ub; j++) {
        beta1->data[j].re = rtNaN;
        beta1->data[j].im = 0.0;
      }

      j = V->size[0] * V->size[1];
      V->size[0] = A->size[0];
      V->size[1] = A->size[0];
      emxEnsureCapacity((emxArray__common *)V, j, sizeof(creal_T));
      loop_ub = A->size[0] * A->size[0];
      for (j = 0; j < loop_ub; j++) {
        V->data[j].re = rtNaN;
        V->data[j].im = 0.0;
      }
    } else {
      ilascl = false;
      anrmto = anrm;
      if ((anrm > 0.0) && (anrm < 6.7178761075670888E-139)) {
        anrmto = 6.7178761075670888E-139;
        ilascl = true;
      } else {
        if (anrm > 1.4885657073574029E+138) {
          anrmto = 1.4885657073574029E+138;
          ilascl = true;
        }
      }

      if (ilascl) {
        absxk = anrm;
        ctoc = anrmto;
        notdone = true;
        while (notdone) {
          cfrom1 = absxk * 2.0041683600089728E-292;
          cto1 = ctoc / 4.9896007738368E+291;
          if ((cfrom1 > ctoc) && (ctoc != 0.0)) {
            stemp_im = 2.0041683600089728E-292;
            absxk = cfrom1;
          } else if (cto1 > absxk) {
            stemp_im = 4.9896007738368E+291;
            ctoc = cto1;
          } else {
            stemp_im = ctoc / absxk;
            notdone = false;
          }

          j = A->size[0] * A->size[1];
          emxEnsureCapacity((emxArray__common *)A, j, sizeof(creal_T));
          loop_ub = A->size[1];
          for (j = 0; j < loop_ub; j++) {
            jcol = A->size[0];
            for (jrow = 0; jrow < jcol; jrow++) {
              A->data[jrow + A->size[0] * j].re *= stemp_im;
              A->data[jrow + A->size[0] * j].im *= stemp_im;
            }
          }
        }
      }

      emxInit_int32_T(&rscale, 1);
      emxInit_int8_T(&I, 2);
      xzggbal(A, &ilo, &ihi, rscale);
      b_n = A->size[0];
      j = I->size[0] * I->size[1];
      I->size[0] = A->size[0];
      I->size[1] = A->size[0];
      emxEnsureCapacity((emxArray__common *)I, j, sizeof(signed char));
      loop_ub = A->size[0] * A->size[0];
      for (j = 0; j < loop_ub; j++) {
        I->data[j] = 0;
      }

      if (A->size[0] > 0) {
        for (jcol = 0; jcol + 1 <= b_n; jcol++) {
          I->data[jcol + I->size[0] * jcol] = 1;
        }
      }

      j = V->size[0] * V->size[1];
      V->size[0] = I->size[0];
      V->size[1] = I->size[1];
      emxEnsureCapacity((emxArray__common *)V, j, sizeof(creal_T));
      loop_ub = I->size[0] * I->size[1];
      for (j = 0; j < loop_ub; j++) {
        V->data[j].re = I->data[j];
        V->data[j].im = 0.0;
      }

      emxFree_int8_T(&I);
      if ((!(A->size[0] <= 1)) && (!(ihi < ilo + 2))) {
        for (jcol = ilo - 1; jcol + 1 < ihi - 1; jcol++) {
          for (jrow = ihi - 1; jrow + 1 > jcol + 2; jrow--) {
            b_A = A->data[(jrow + A->size[0] * jcol) - 1];
            c_A = A->data[jrow + A->size[0] * jcol];
            xzlartg(b_A, c_A, &c, &tmp, &A->data[(jrow + A->size[0] * jcol) - 1]);
            A->data[jrow + A->size[0] * jcol].re = 0.0;
            A->data[jrow + A->size[0] * jcol].im = 0.0;
            for (j = jcol + 1; j + 1 <= b_n; j++) {
              absxk = tmp.re * A->data[jrow + A->size[0] * j].re - tmp.im *
                A->data[jrow + A->size[0] * j].im;
              ctoc = tmp.re * A->data[jrow + A->size[0] * j].im + tmp.im *
                A->data[jrow + A->size[0] * j].re;
              stemp_re = c * A->data[(jrow + A->size[0] * j) - 1].re + absxk;
              stemp_im = c * A->data[(jrow + A->size[0] * j) - 1].im + ctoc;
              absxk = A->data[(jrow + A->size[0] * j) - 1].re;
              ctoc = A->data[(jrow + A->size[0] * j) - 1].im;
              cfrom1 = A->data[(jrow + A->size[0] * j) - 1].im;
              cto1 = A->data[(jrow + A->size[0] * j) - 1].re;
              A->data[jrow + A->size[0] * j].re = c * A->data[jrow + A->size[0] *
                j].re - (tmp.re * absxk + tmp.im * ctoc);
              A->data[jrow + A->size[0] * j].im = c * A->data[jrow + A->size[0] *
                j].im - (tmp.re * cfrom1 - tmp.im * cto1);
              A->data[(jrow + A->size[0] * j) - 1].re = stemp_re;
              A->data[(jrow + A->size[0] * j) - 1].im = stemp_im;
            }

            tmp.re = -tmp.re;
            tmp.im = -tmp.im;
            for (loop_ub = 0; loop_ub + 1 <= ihi; loop_ub++) {
              absxk = tmp.re * A->data[loop_ub + A->size[0] * (jrow - 1)].re -
                tmp.im * A->data[loop_ub + A->size[0] * (jrow - 1)].im;
              ctoc = tmp.re * A->data[loop_ub + A->size[0] * (jrow - 1)].im +
                tmp.im * A->data[loop_ub + A->size[0] * (jrow - 1)].re;
              stemp_re = c * A->data[loop_ub + A->size[0] * jrow].re + absxk;
              stemp_im = c * A->data[loop_ub + A->size[0] * jrow].im + ctoc;
              absxk = A->data[loop_ub + A->size[0] * jrow].re;
              ctoc = A->data[loop_ub + A->size[0] * jrow].im;
              cfrom1 = A->data[loop_ub + A->size[0] * jrow].im;
              cto1 = A->data[loop_ub + A->size[0] * jrow].re;
              A->data[loop_ub + A->size[0] * (jrow - 1)].re = c * A->
                data[loop_ub + A->size[0] * (jrow - 1)].re - (tmp.re * absxk +
                tmp.im * ctoc);
              A->data[loop_ub + A->size[0] * (jrow - 1)].im = c * A->
                data[loop_ub + A->size[0] * (jrow - 1)].im - (tmp.re * cfrom1 -
                tmp.im * cto1);
              A->data[loop_ub + A->size[0] * jrow].re = stemp_re;
              A->data[loop_ub + A->size[0] * jrow].im = stemp_im;
            }

            for (loop_ub = 0; loop_ub + 1 <= b_n; loop_ub++) {
              absxk = tmp.re * V->data[loop_ub + V->size[0] * (jrow - 1)].re -
                tmp.im * V->data[loop_ub + V->size[0] * (jrow - 1)].im;
              ctoc = tmp.re * V->data[loop_ub + V->size[0] * (jrow - 1)].im +
                tmp.im * V->data[loop_ub + V->size[0] * (jrow - 1)].re;
              stemp_re = c * V->data[loop_ub + V->size[0] * jrow].re + absxk;
              stemp_im = c * V->data[loop_ub + V->size[0] * jrow].im + ctoc;
              absxk = V->data[loop_ub + V->size[0] * jrow].re;
              ctoc = V->data[loop_ub + V->size[0] * jrow].im;
              cfrom1 = V->data[loop_ub + V->size[0] * jrow].im;
              cto1 = V->data[loop_ub + V->size[0] * jrow].re;
              V->data[loop_ub + V->size[0] * (jrow - 1)].re = c * V->
                data[loop_ub + V->size[0] * (jrow - 1)].re - (tmp.re * absxk +
                tmp.im * ctoc);
              V->data[loop_ub + V->size[0] * (jrow - 1)].im = c * V->
                data[loop_ub + V->size[0] * (jrow - 1)].im - (tmp.re * cfrom1 -
                tmp.im * cto1);
              V->data[loop_ub + V->size[0] * jrow].re = stemp_re;
              V->data[loop_ub + V->size[0] * jrow].im = stemp_im;
            }
          }
        }
      }

      xzhgeqz(A, ilo, ihi, V, &b_info, alpha1, beta1);
      if (b_info == 0) {
        xztgevc(A, V);
        b_n = V->size[0];
        jrow = V->size[1];
        if (ilo > 1) {
          for (loop_ub = ilo - 2; loop_ub + 1 >= 1; loop_ub--) {
            jcol = rscale->data[loop_ub] - 1;
            if (rscale->data[loop_ub] != loop_ub + 1) {
              for (j = 0; j + 1 <= jrow; j++) {
                tmp = V->data[loop_ub + V->size[0] * j];
                V->data[loop_ub + V->size[0] * j] = V->data[jcol + V->size[0] *
                  j];
                V->data[jcol + V->size[0] * j] = tmp;
              }
            }
          }
        }

        if (ihi < b_n) {
          while (ihi + 1 <= b_n) {
            jcol = rscale->data[ihi] - 1;
            if (rscale->data[ihi] != ihi + 1) {
              for (j = 0; j + 1 <= jrow; j++) {
                tmp = V->data[ihi + V->size[0] * j];
                V->data[ihi + V->size[0] * j] = V->data[jcol + V->size[0] * j];
                V->data[jcol + V->size[0] * j] = tmp;
              }
            }

            ihi++;
          }
        }

        for (jcol = 0; jcol < n; jcol++) {
          absxk = std::abs(V->data[V->size[0] * jcol].re) + std::abs(V->data
            [V->size[0] * jcol].im);
          if (n > 1) {
            for (jrow = 1; jrow - 1 <= n - 2; jrow++) {
              ctoc = std::abs(V->data[jrow + V->size[0] * jcol].re) + std::abs
                (V->data[jrow + V->size[0] * jcol].im);
              if (ctoc > absxk) {
                absxk = ctoc;
              }
            }
          }

          if (absxk >= 6.7178761075670888E-139) {
            absxk = 1.0 / absxk;
            for (jrow = 0; jrow < n; jrow++) {
              V->data[jrow + V->size[0] * jcol].re *= absxk;
              V->data[jrow + V->size[0] * jcol].im *= absxk;
            }
          }
        }

        if (ilascl) {
          notdone = true;
          while (notdone) {
            cfrom1 = anrmto * 2.0041683600089728E-292;
            cto1 = anrm / 4.9896007738368E+291;
            if ((cfrom1 > anrm) && (anrm != 0.0)) {
              stemp_im = 2.0041683600089728E-292;
              anrmto = cfrom1;
            } else if (cto1 > anrmto) {
              stemp_im = 4.9896007738368E+291;
              anrm = cto1;
            } else {
              stemp_im = anrm / anrmto;
              notdone = false;
            }

            j = alpha1->size[0];
            emxEnsureCapacity((emxArray__common *)alpha1, j, sizeof(creal_T));
            loop_ub = alpha1->size[0];
            for (j = 0; j < loop_ub; j++) {
              alpha1->data[j].re *= stemp_im;
              alpha1->data[j].im *= stemp_im;
            }
          }
        }
      }

      emxFree_int32_T(&rscale);
    }
  }

  *info = b_info;
}

/* End of code generation (xzggev.cpp) */
