/*
 * Student License - for use by students to meet course requirements and
 * perform academic research at degree granting institutions only.  Not
 * for government, commercial, or other organizational use.
 *
 * xzggbal.cpp
 *
 * Code generation for function 'xzggbal'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "MAVEfast.h"
#include "xzggbal.h"
#include "MAVEfast_emxutil.h"

/* Function Definitions */
void xzggbal(emxArray_creal_T *A, int *ilo, int *ihi, emxArray_int32_T *rscale)
{
  int nzcount;
  int ii;
  int b_ilo;
  int b_ihi;
  int exitg2;
  int i;
  int j;
  boolean_T found;
  boolean_T exitg3;
  int jj;
  boolean_T exitg4;
  double atmp_re;
  int exitg1;
  double atmp_im;
  boolean_T b_A;
  nzcount = rscale->size[0];
  rscale->size[0] = A->size[0];
  emxEnsureCapacity((emxArray__common *)rscale, nzcount, sizeof(int));
  ii = A->size[0];
  for (nzcount = 0; nzcount < ii; nzcount++) {
    rscale->data[nzcount] = 1;
  }

  b_ilo = 0;
  b_ihi = A->size[0];
  if (A->size[0] <= 1) {
    b_ihi = 1;
  } else {
    do {
      exitg2 = 0;
      i = 0;
      j = 0;
      found = false;
      ii = b_ihi;
      exitg3 = false;
      while ((!exitg3) && (ii > 0)) {
        nzcount = 0;
        i = ii;
        j = b_ihi;
        jj = 1;
        exitg4 = false;
        while ((!exitg4) && (jj <= b_ihi)) {
          b_A = ((A->data[(ii + A->size[0] * (jj - 1)) - 1].re != 0.0) ||
                 (A->data[(ii + A->size[0] * (jj - 1)) - 1].im != 0.0));
          if (b_A || (ii == jj)) {
            if (nzcount == 0) {
              j = jj;
              nzcount = 1;
              jj++;
            } else {
              nzcount = 2;
              exitg4 = true;
            }
          } else {
            jj++;
          }
        }

        if (nzcount < 2) {
          found = true;
          exitg3 = true;
        } else {
          ii--;
        }
      }

      if (!found) {
        exitg2 = 2;
      } else {
        nzcount = A->size[0];
        if (i != b_ihi) {
          for (ii = 0; ii + 1 <= nzcount; ii++) {
            atmp_re = A->data[(i + A->size[0] * ii) - 1].re;
            atmp_im = A->data[(i + A->size[0] * ii) - 1].im;
            A->data[(i + A->size[0] * ii) - 1] = A->data[(b_ihi + A->size[0] *
              ii) - 1];
            A->data[(b_ihi + A->size[0] * ii) - 1].re = atmp_re;
            A->data[(b_ihi + A->size[0] * ii) - 1].im = atmp_im;
          }
        }

        if (j != b_ihi) {
          for (ii = 0; ii + 1 <= b_ihi; ii++) {
            atmp_re = A->data[ii + A->size[0] * (j - 1)].re;
            atmp_im = A->data[ii + A->size[0] * (j - 1)].im;
            A->data[ii + A->size[0] * (j - 1)] = A->data[ii + A->size[0] *
              (b_ihi - 1)];
            A->data[ii + A->size[0] * (b_ihi - 1)].re = atmp_re;
            A->data[ii + A->size[0] * (b_ihi - 1)].im = atmp_im;
          }
        }

        rscale->data[b_ihi - 1] = j;
        b_ihi--;
        if (b_ihi == 1) {
          rscale->data[0] = 1;
          exitg2 = 1;
        }
      }
    } while (exitg2 == 0);

    if (exitg2 == 1) {
    } else {
      do {
        exitg1 = 0;
        i = 0;
        j = 0;
        found = false;
        jj = b_ilo + 1;
        exitg3 = false;
        while ((!exitg3) && (jj <= b_ihi)) {
          nzcount = 0;
          i = b_ihi;
          j = jj;
          ii = b_ilo + 1;
          exitg4 = false;
          while ((!exitg4) && (ii <= b_ihi)) {
            b_A = ((A->data[(ii + A->size[0] * (jj - 1)) - 1].re != 0.0) ||
                   (A->data[(ii + A->size[0] * (jj - 1)) - 1].im != 0.0));
            if (b_A || (ii == jj)) {
              if (nzcount == 0) {
                i = ii;
                nzcount = 1;
                ii++;
              } else {
                nzcount = 2;
                exitg4 = true;
              }
            } else {
              ii++;
            }
          }

          if (nzcount < 2) {
            found = true;
            exitg3 = true;
          } else {
            jj++;
          }
        }

        if (!found) {
          exitg1 = 1;
        } else {
          nzcount = A->size[0];
          if (i != b_ilo + 1) {
            for (ii = b_ilo; ii + 1 <= nzcount; ii++) {
              atmp_re = A->data[(i + A->size[0] * ii) - 1].re;
              atmp_im = A->data[(i + A->size[0] * ii) - 1].im;
              A->data[(i + A->size[0] * ii) - 1] = A->data[b_ilo + A->size[0] *
                ii];
              A->data[b_ilo + A->size[0] * ii].re = atmp_re;
              A->data[b_ilo + A->size[0] * ii].im = atmp_im;
            }
          }

          if (j != b_ilo + 1) {
            for (ii = 0; ii + 1 <= b_ihi; ii++) {
              atmp_re = A->data[ii + A->size[0] * (j - 1)].re;
              atmp_im = A->data[ii + A->size[0] * (j - 1)].im;
              A->data[ii + A->size[0] * (j - 1)] = A->data[ii + A->size[0] *
                b_ilo];
              A->data[ii + A->size[0] * b_ilo].re = atmp_re;
              A->data[ii + A->size[0] * b_ilo].im = atmp_im;
            }
          }

          rscale->data[b_ilo] = j;
          b_ilo++;
          if (b_ilo + 1 == b_ihi) {
            rscale->data[b_ilo] = b_ilo + 1;
            exitg1 = 1;
          }
        }
      } while (exitg1 == 0);
    }
  }

  *ilo = b_ilo + 1;
  *ihi = b_ihi;
}

/* End of code generation (xzggbal.cpp) */
