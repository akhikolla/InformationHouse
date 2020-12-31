/*
 * Student License - for use by students to meet course requirements and
 * perform academic research at degree granting institutions only.  Not
 * for government, commercial, or other organizational use.
 *
 * xzhgeqz.cpp
 *
 * Code generation for function 'xzhgeqz'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "MAVEfast.h"
#include "xzhgeqz.h"
#include "MAVEfast_emxutil.h"
#include "xzlartg.h"
#include "sqrt.h"

/* Function Definitions */
void xzhgeqz(emxArray_creal_T *A, int ilo, int ihi, emxArray_creal_T *Z, int
             *info, emxArray_creal_T *alpha1, emxArray_creal_T *beta1)
{
  int b_info;
  boolean_T compz;
  int n;
  int ifirst;
  int loop_ub;
  double eshift_re;
  double eshift_im;
  creal_T ctemp;
  double anorm;
  double scale;
  double reAij;
  double sumsq;
  double b_atol;
  boolean_T firstNonZero;
  int j;
  double ascale;
  double bscale;
  int jp1;
  boolean_T failed;
  double imAij;
  boolean_T guard1 = false;
  boolean_T guard2 = false;
  int istart;
  double temp2;
  int ilast;
  int ilastm1;
  int ifrstm;
  int ilastm;
  int iiter;
  boolean_T goto60;
  boolean_T goto70;
  boolean_T goto90;
  int jiter;
  int exitg1;
  int jm1;
  boolean_T exitg2;
  creal_T b_ascale;
  creal_T shift;
  creal_T b_A;
  double ad22_re;
  double ad22_im;
  double t1_im;
  b_info = -1;
  compz = !((Z->size[0] == 0) || (Z->size[1] == 0));
  if ((A->size[0] == 1) && (A->size[1] == 1)) {
    ihi = 1;
  }

  n = A->size[0];
  ifirst = alpha1->size[0];
  alpha1->size[0] = A->size[0];
  emxEnsureCapacity((emxArray__common *)alpha1, ifirst, sizeof(creal_T));
  loop_ub = A->size[0];
  for (ifirst = 0; ifirst < loop_ub; ifirst++) {
    alpha1->data[ifirst].re = 0.0;
    alpha1->data[ifirst].im = 0.0;
  }

  ifirst = beta1->size[0];
  beta1->size[0] = A->size[0];
  emxEnsureCapacity((emxArray__common *)beta1, ifirst, sizeof(creal_T));
  loop_ub = A->size[0];
  for (ifirst = 0; ifirst < loop_ub; ifirst++) {
    beta1->data[ifirst].re = 1.0;
    beta1->data[ifirst].im = 0.0;
  }

  eshift_re = 0.0;
  eshift_im = 0.0;
  ctemp.re = 0.0;
  ctemp.im = 0.0;
  anorm = 0.0;
  if (!(ilo > ihi)) {
    scale = 0.0;
    sumsq = 0.0;
    firstNonZero = true;
    for (j = ilo; j <= ihi; j++) {
      ifirst = j + 1;
      if (ihi < j + 1) {
        ifirst = ihi;
      }

      for (jp1 = ilo; jp1 <= ifirst; jp1++) {
        reAij = A->data[(jp1 + A->size[0] * (j - 1)) - 1].re;
        imAij = A->data[(jp1 + A->size[0] * (j - 1)) - 1].im;
        if (reAij != 0.0) {
          anorm = std::abs(reAij);
          if (firstNonZero) {
            sumsq = 1.0;
            scale = anorm;
            firstNonZero = false;
          } else if (scale < anorm) {
            temp2 = scale / anorm;
            sumsq = 1.0 + sumsq * temp2 * temp2;
            scale = anorm;
          } else {
            temp2 = anorm / scale;
            sumsq += temp2 * temp2;
          }
        }

        if (imAij != 0.0) {
          anorm = std::abs(imAij);
          if (firstNonZero) {
            sumsq = 1.0;
            scale = anorm;
            firstNonZero = false;
          } else if (scale < anorm) {
            temp2 = scale / anorm;
            sumsq = 1.0 + sumsq * temp2 * temp2;
            scale = anorm;
          } else {
            temp2 = anorm / scale;
            sumsq += temp2 * temp2;
          }
        }
      }
    }

    anorm = scale * std::sqrt(sumsq);
  }

  reAij = 2.2204460492503131E-16 * anorm;
  b_atol = 2.2250738585072014E-308;
  if (reAij > 2.2250738585072014E-308) {
    b_atol = reAij;
  }

  reAij = 2.2250738585072014E-308;
  if (anorm > 2.2250738585072014E-308) {
    reAij = anorm;
  }

  ascale = 1.0 / reAij;
  bscale = 1.0 / std::sqrt((double)A->size[0]);
  failed = true;
  for (j = ihi; j + 1 <= n; j++) {
    alpha1->data[j] = A->data[j + A->size[0] * j];
  }

  guard1 = false;
  guard2 = false;
  if (ihi >= ilo) {
    ifirst = ilo;
    istart = ilo;
    ilast = ihi - 1;
    ilastm1 = ihi - 2;
    if (compz) {
      ifrstm = 1;
      ilastm = n;
    } else {
      ifrstm = ilo;
      ilastm = ihi;
    }

    iiter = 0;
    goto60 = false;
    goto70 = false;
    goto90 = false;
    jiter = 1;
    do {
      exitg1 = 0;
      if (jiter <= 30 * ((ihi - ilo) + 1)) {
        if (ilast + 1 == ilo) {
          goto60 = true;
        } else if (std::abs(A->data[ilast + A->size[0] * ilastm1].re) + std::abs
                   (A->data[ilast + A->size[0] * ilastm1].im) <= b_atol) {
          A->data[ilast + A->size[0] * ilastm1].re = 0.0;
          A->data[ilast + A->size[0] * ilastm1].im = 0.0;
          goto60 = true;
        } else {
          j = ilastm1;
          exitg2 = false;
          while ((!exitg2) && (j + 1 >= ilo)) {
            if (j + 1 == ilo) {
              firstNonZero = true;
            } else if (std::abs(A->data[j + A->size[0] * (j - 1)].re) + std::abs
                       (A->data[j + A->size[0] * (j - 1)].im) <= b_atol) {
              A->data[j + A->size[0] * (j - 1)].re = 0.0;
              A->data[j + A->size[0] * (j - 1)].im = 0.0;
              firstNonZero = true;
            } else {
              firstNonZero = false;
            }

            if (firstNonZero) {
              ifirst = j + 1;
              goto70 = true;
              exitg2 = true;
            } else {
              j--;
            }
          }
        }

        if (goto60 || goto70) {
          firstNonZero = true;
        } else {
          firstNonZero = false;
        }

        if (!firstNonZero) {
          jp1 = alpha1->size[0];
          ifirst = alpha1->size[0];
          alpha1->size[0] = jp1;
          emxEnsureCapacity((emxArray__common *)alpha1, ifirst, sizeof(creal_T));
          for (ifirst = 0; ifirst < jp1; ifirst++) {
            alpha1->data[ifirst].re = rtNaN;
            alpha1->data[ifirst].im = 0.0;
          }

          jp1 = beta1->size[0];
          ifirst = beta1->size[0];
          beta1->size[0] = jp1;
          emxEnsureCapacity((emxArray__common *)beta1, ifirst, sizeof(creal_T));
          for (ifirst = 0; ifirst < jp1; ifirst++) {
            beta1->data[ifirst].re = rtNaN;
            beta1->data[ifirst].im = 0.0;
          }

          if (compz) {
            ifirst = Z->size[0] * Z->size[1];
            emxEnsureCapacity((emxArray__common *)Z, ifirst, sizeof(creal_T));
            loop_ub = Z->size[1];
            for (ifirst = 0; ifirst < loop_ub; ifirst++) {
              jp1 = Z->size[0];
              for (jm1 = 0; jm1 < jp1; jm1++) {
                Z->data[jm1 + Z->size[0] * ifirst].re = rtNaN;
                Z->data[jm1 + Z->size[0] * ifirst].im = 0.0;
              }
            }
          }

          b_info = 0;
          exitg1 = 1;
        } else if (goto60) {
          goto60 = false;
          alpha1->data[ilast] = A->data[ilast + A->size[0] * ilast];
          ilast = ilastm1;
          ilastm1--;
          if (ilast + 1 < ilo) {
            failed = false;
            guard2 = true;
            exitg1 = 1;
          } else {
            iiter = 0;
            eshift_re = 0.0;
            eshift_im = 0.0;
            if (!compz) {
              ilastm = ilast + 1;
              if (ifrstm > ilast + 1) {
                ifrstm = ilo;
              }
            }

            jiter++;
          }
        } else {
          if (goto70) {
            goto70 = false;
            iiter++;
            if (!compz) {
              ifrstm = ifirst;
            }

            if (iiter - iiter / 10 * 10 != 0) {
              anorm = ascale * A->data[ilastm1 + A->size[0] * ilastm1].re;
              reAij = ascale * A->data[ilastm1 + A->size[0] * ilastm1].im;
              if (reAij == 0.0) {
                shift.re = anorm / bscale;
                shift.im = 0.0;
              } else if (anorm == 0.0) {
                shift.re = 0.0;
                shift.im = reAij / bscale;
              } else {
                shift.re = anorm / bscale;
                shift.im = reAij / bscale;
              }

              anorm = ascale * A->data[ilast + A->size[0] * ilast].re;
              reAij = ascale * A->data[ilast + A->size[0] * ilast].im;
              if (reAij == 0.0) {
                ad22_re = anorm / bscale;
                ad22_im = 0.0;
              } else if (anorm == 0.0) {
                ad22_re = 0.0;
                ad22_im = reAij / bscale;
              } else {
                ad22_re = anorm / bscale;
                ad22_im = reAij / bscale;
              }

              temp2 = 0.5 * (shift.re + ad22_re);
              t1_im = 0.5 * (shift.im + ad22_im);
              anorm = ascale * A->data[ilastm1 + A->size[0] * ilast].re;
              reAij = ascale * A->data[ilastm1 + A->size[0] * ilast].im;
              if (reAij == 0.0) {
                sumsq = anorm / bscale;
                imAij = 0.0;
              } else if (anorm == 0.0) {
                sumsq = 0.0;
                imAij = reAij / bscale;
              } else {
                sumsq = anorm / bscale;
                imAij = reAij / bscale;
              }

              anorm = ascale * A->data[ilast + A->size[0] * ilastm1].re;
              reAij = ascale * A->data[ilast + A->size[0] * ilastm1].im;
              if (reAij == 0.0) {
                scale = anorm / bscale;
                anorm = 0.0;
              } else if (anorm == 0.0) {
                scale = 0.0;
                anorm = reAij / bscale;
              } else {
                scale = anorm / bscale;
                anorm = reAij / bscale;
              }

              reAij = shift.re * ad22_im + shift.im * ad22_re;
              shift.re = ((temp2 * temp2 - t1_im * t1_im) + (sumsq * scale -
                imAij * anorm)) - (shift.re * ad22_re - shift.im * ad22_im);
              shift.im = ((temp2 * t1_im + t1_im * temp2) + (sumsq * anorm +
                imAij * scale)) - reAij;
              c_sqrt(&shift);
              if ((temp2 - ad22_re) * shift.re + (t1_im - ad22_im) * shift.im <=
                  0.0) {
                shift.re += temp2;
                shift.im += t1_im;
              } else {
                shift.re = temp2 - shift.re;
                shift.im = t1_im - shift.im;
              }
            } else {
              anorm = ascale * A->data[ilast + A->size[0] * ilastm1].re;
              reAij = ascale * A->data[ilast + A->size[0] * ilastm1].im;
              if (reAij == 0.0) {
                sumsq = anorm / bscale;
                imAij = 0.0;
              } else if (anorm == 0.0) {
                sumsq = 0.0;
                imAij = reAij / bscale;
              } else {
                sumsq = anorm / bscale;
                imAij = reAij / bscale;
              }

              eshift_re += sumsq;
              eshift_im += imAij;
              shift.re = eshift_re;
              shift.im = eshift_im;
            }

            j = ilastm1;
            jp1 = ilastm1 + 1;
            exitg2 = false;
            while ((!exitg2) && (j + 1 > ifirst)) {
              istart = j + 1;
              ctemp.re = ascale * A->data[j + A->size[0] * j].re - shift.re *
                bscale;
              ctemp.im = ascale * A->data[j + A->size[0] * j].im - shift.im *
                bscale;
              anorm = std::abs(ctemp.re) + std::abs(ctemp.im);
              temp2 = ascale * (std::abs(A->data[jp1 + A->size[0] * j].re) + std::
                                abs(A->data[jp1 + A->size[0] * j].im));
              reAij = anorm;
              if (temp2 > anorm) {
                reAij = temp2;
              }

              if ((reAij < 1.0) && (reAij != 0.0)) {
                anorm /= reAij;
                temp2 /= reAij;
              }

              if ((std::abs(A->data[j + A->size[0] * (j - 1)].re) + std::abs
                   (A->data[j + A->size[0] * (j - 1)].im)) * temp2 <= anorm *
                  b_atol) {
                goto90 = true;
                exitg2 = true;
              } else {
                jp1 = j;
                j--;
              }
            }

            if (!goto90) {
              istart = ifirst;
              ctemp.re = ascale * A->data[(ifirst + A->size[0] * (ifirst - 1)) -
                1].re - shift.re * bscale;
              ctemp.im = ascale * A->data[(ifirst + A->size[0] * (ifirst - 1)) -
                1].im - shift.im * bscale;
              goto90 = true;
            }
          }

          if (goto90) {
            goto90 = false;
            b_ascale.re = ascale * A->data[istart + A->size[0] * (istart - 1)].
              re;
            b_ascale.im = ascale * A->data[istart + A->size[0] * (istart - 1)].
              im;
            b_xzlartg(ctemp, b_ascale, &imAij, &shift);
            j = istart;
            jm1 = istart - 2;
            while (j < ilast + 1) {
              if (j > istart) {
                b_ascale = A->data[(j + A->size[0] * jm1) - 1];
                b_A = A->data[j + A->size[0] * jm1];
                xzlartg(b_ascale, b_A, &imAij, &shift, &A->data[(j + A->size[0] *
                         jm1) - 1]);
                A->data[j + A->size[0] * jm1].re = 0.0;
                A->data[j + A->size[0] * jm1].im = 0.0;
              }

              for (loop_ub = j - 1; loop_ub + 1 <= ilastm; loop_ub++) {
                anorm = shift.re * A->data[j + A->size[0] * loop_ub].re -
                  shift.im * A->data[j + A->size[0] * loop_ub].im;
                reAij = shift.re * A->data[j + A->size[0] * loop_ub].im +
                  shift.im * A->data[j + A->size[0] * loop_ub].re;
                ad22_re = imAij * A->data[(j + A->size[0] * loop_ub) - 1].re +
                  anorm;
                ad22_im = imAij * A->data[(j + A->size[0] * loop_ub) - 1].im +
                  reAij;
                anorm = A->data[(j + A->size[0] * loop_ub) - 1].re;
                reAij = A->data[(j + A->size[0] * loop_ub) - 1].im;
                scale = A->data[(j + A->size[0] * loop_ub) - 1].im;
                sumsq = A->data[(j + A->size[0] * loop_ub) - 1].re;
                A->data[j + A->size[0] * loop_ub].re = imAij * A->data[j +
                  A->size[0] * loop_ub].re - (shift.re * anorm + shift.im *
                  reAij);
                A->data[j + A->size[0] * loop_ub].im = imAij * A->data[j +
                  A->size[0] * loop_ub].im - (shift.re * scale - shift.im *
                  sumsq);
                A->data[(j + A->size[0] * loop_ub) - 1].re = ad22_re;
                A->data[(j + A->size[0] * loop_ub) - 1].im = ad22_im;
              }

              shift.re = -shift.re;
              shift.im = -shift.im;
              loop_ub = j;
              if (ilast + 1 < j + 2) {
                loop_ub = ilast - 1;
              }

              for (jp1 = ifrstm - 1; jp1 + 1 <= loop_ub + 2; jp1++) {
                anorm = shift.re * A->data[jp1 + A->size[0] * (j - 1)].re -
                  shift.im * A->data[jp1 + A->size[0] * (j - 1)].im;
                reAij = shift.re * A->data[jp1 + A->size[0] * (j - 1)].im +
                  shift.im * A->data[jp1 + A->size[0] * (j - 1)].re;
                ad22_re = imAij * A->data[jp1 + A->size[0] * j].re + anorm;
                ad22_im = imAij * A->data[jp1 + A->size[0] * j].im + reAij;
                anorm = A->data[jp1 + A->size[0] * j].re;
                reAij = A->data[jp1 + A->size[0] * j].im;
                scale = A->data[jp1 + A->size[0] * j].im;
                sumsq = A->data[jp1 + A->size[0] * j].re;
                A->data[jp1 + A->size[0] * (j - 1)].re = imAij * A->data[jp1 +
                  A->size[0] * (j - 1)].re - (shift.re * anorm + shift.im *
                  reAij);
                A->data[jp1 + A->size[0] * (j - 1)].im = imAij * A->data[jp1 +
                  A->size[0] * (j - 1)].im - (shift.re * scale - shift.im *
                  sumsq);
                A->data[jp1 + A->size[0] * j].re = ad22_re;
                A->data[jp1 + A->size[0] * j].im = ad22_im;
              }

              if (compz) {
                for (jp1 = 0; jp1 + 1 <= n; jp1++) {
                  anorm = shift.re * Z->data[jp1 + Z->size[0] * (j - 1)].re -
                    shift.im * Z->data[jp1 + Z->size[0] * (j - 1)].im;
                  reAij = shift.re * Z->data[jp1 + Z->size[0] * (j - 1)].im +
                    shift.im * Z->data[jp1 + Z->size[0] * (j - 1)].re;
                  ad22_re = imAij * Z->data[jp1 + Z->size[0] * j].re + anorm;
                  ad22_im = imAij * Z->data[jp1 + Z->size[0] * j].im + reAij;
                  anorm = Z->data[jp1 + Z->size[0] * j].re;
                  reAij = Z->data[jp1 + Z->size[0] * j].im;
                  scale = Z->data[jp1 + Z->size[0] * j].im;
                  sumsq = Z->data[jp1 + Z->size[0] * j].re;
                  Z->data[jp1 + Z->size[0] * (j - 1)].re = imAij * Z->data[jp1 +
                    Z->size[0] * (j - 1)].re - (shift.re * anorm + shift.im *
                    reAij);
                  Z->data[jp1 + Z->size[0] * (j - 1)].im = imAij * Z->data[jp1 +
                    Z->size[0] * (j - 1)].im - (shift.re * scale - shift.im *
                    sumsq);
                  Z->data[jp1 + Z->size[0] * j].re = ad22_re;
                  Z->data[jp1 + Z->size[0] * j].im = ad22_im;
                }
              }

              jm1 = j - 1;
              j++;
            }
          }

          jiter++;
        }
      } else {
        guard2 = true;
        exitg1 = 1;
      }
    } while (exitg1 == 0);
  } else {
    guard1 = true;
  }

  if (guard2) {
    if (failed) {
      b_info = ilast;
      for (jp1 = 0; jp1 + 1 <= ilast + 1; jp1++) {
        alpha1->data[jp1].re = rtNaN;
        alpha1->data[jp1].im = 0.0;
        beta1->data[jp1].re = rtNaN;
        beta1->data[jp1].im = 0.0;
      }

      if (compz) {
        ifirst = Z->size[0] * Z->size[1];
        emxEnsureCapacity((emxArray__common *)Z, ifirst, sizeof(creal_T));
        loop_ub = Z->size[1];
        for (ifirst = 0; ifirst < loop_ub; ifirst++) {
          jp1 = Z->size[0];
          for (jm1 = 0; jm1 < jp1; jm1++) {
            Z->data[jm1 + Z->size[0] * ifirst].re = rtNaN;
            Z->data[jm1 + Z->size[0] * ifirst].im = 0.0;
          }
        }
      }
    } else {
      guard1 = true;
    }
  }

  if (guard1) {
    for (j = 0; j + 1 < ilo; j++) {
      alpha1->data[j] = A->data[j + A->size[0] * j];
    }
  }

  *info = b_info + 1;
}

/* End of code generation (xzhgeqz.cpp) */
