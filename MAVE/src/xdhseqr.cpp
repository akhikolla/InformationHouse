/*
 * Student License - for use by students to meet course requirements and
 * perform academic research at degree granting institutions only.  Not
 * for government, commercial, or other organizational use.
 *
 * xdhseqr.cpp
 *
 * Code generation for function 'xdhseqr'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "MAVEfast.h"
#include "xdhseqr.h"
#include "xrot.h"
#include "xdlanv2.h"
#include "xzlarfg.h"

/* Function Definitions */
int eml_dlahqr(emxArray_real_T *h, emxArray_real_T *z)
{
  int info;
  int n;
  int ldh;
  int ldz;
  int j;
  double SMLNUM;
  int i;
  boolean_T exitg1;
  int L;
  boolean_T goto150;
  int its;
  boolean_T exitg2;
  int k;
  boolean_T exitg3;
  double tst;
  double htmp1;
  double aa;
  double ab;
  double ba;
  double rt2r;
  double rt1r;
  double s;
  double sn;
  int m;
  double b_SMLNUM;
  int b_k;
  double v[3];
  int nr;
  int hoffset;
  n = h->size[0];
  ldh = h->size[0];
  ldz = z->size[0];
  info = 0;
  if ((h->size[0] != 0) && (1 != h->size[0])) {
    for (j = 0; j + 1 <= n - 3; j++) {
      h->data[(j + h->size[0] * j) + 2] = 0.0;
      h->data[(j + h->size[0] * j) + 3] = 0.0;
    }

    if (1 <= n - 2) {
      h->data[(n + h->size[0] * (n - 3)) - 1] = 0.0;
    }

    SMLNUM = 2.2250738585072014E-308 * ((double)n / 2.2204460492503131E-16);
    i = n - 1;
    exitg1 = false;
    while ((!exitg1) && (i + 1 >= 1)) {
      L = 1;
      goto150 = false;
      its = 0;
      exitg2 = false;
      while ((!exitg2) && (its < 31)) {
        k = i;
        exitg3 = false;
        while ((!exitg3) && ((k + 1 > L) && (!(std::abs(h->data[k + h->size[0] *
                   (k - 1)]) <= SMLNUM)))) {
          tst = std::abs(h->data[(k + h->size[0] * (k - 1)) - 1]) + std::abs
            (h->data[k + h->size[0] * k]);
          if (tst == 0.0) {
            if (k - 1 >= 1) {
              tst = std::abs(h->data[(k + h->size[0] * (k - 2)) - 1]);
            }

            if (k + 2 <= n) {
              tst += std::abs(h->data[(k + h->size[0] * k) + 1]);
            }
          }

          if (std::abs(h->data[k + h->size[0] * (k - 1)]) <=
              2.2204460492503131E-16 * tst) {
            htmp1 = std::abs(h->data[k + h->size[0] * (k - 1)]);
            tst = std::abs(h->data[(k + h->size[0] * k) - 1]);
            if (htmp1 > tst) {
              ab = htmp1;
              ba = tst;
            } else {
              ab = tst;
              ba = htmp1;
            }

            htmp1 = std::abs(h->data[k + h->size[0] * k]);
            tst = std::abs(h->data[(k + h->size[0] * (k - 1)) - 1] - h->data[k +
                           h->size[0] * k]);
            if (htmp1 > tst) {
              aa = htmp1;
              htmp1 = tst;
            } else {
              aa = tst;
            }

            s = aa + ab;
            tst = 2.2204460492503131E-16 * (htmp1 * (aa / s));
            if ((SMLNUM > tst) || rtIsNaN(tst)) {
              b_SMLNUM = SMLNUM;
            } else {
              b_SMLNUM = tst;
            }

            if (ba * (ab / s) <= b_SMLNUM) {
              exitg3 = true;
            } else {
              k--;
            }
          } else {
            k--;
          }
        }

        L = k + 1;
        if (k + 1 > 1) {
          h->data[k + h->size[0] * (k - 1)] = 0.0;
        }

        if (k + 1 >= i) {
          goto150 = true;
          exitg2 = true;
        } else {
          if (its == 10) {
            s = std::abs(h->data[(k + h->size[0] * k) + 1]) + std::abs(h->data
              [(k + h->size[0] * (k + 1)) + 2]);
            tst = 0.75 * s + h->data[k + h->size[0] * k];
            aa = -0.4375 * s;
            htmp1 = s;
            ba = tst;
          } else if (its == 20) {
            s = std::abs(h->data[i + h->size[0] * (i - 1)]) + std::abs(h->data
              [(i + h->size[0] * (i - 2)) - 1]);
            tst = 0.75 * s + h->data[i + h->size[0] * i];
            aa = -0.4375 * s;
            htmp1 = s;
            ba = tst;
          } else {
            tst = h->data[(i + h->size[0] * (i - 1)) - 1];
            htmp1 = h->data[i + h->size[0] * (i - 1)];
            aa = h->data[(i + h->size[0] * i) - 1];
            ba = h->data[i + h->size[0] * i];
          }

          s = ((std::abs(tst) + std::abs(aa)) + std::abs(htmp1)) + std::abs(ba);
          if (s == 0.0) {
            rt1r = 0.0;
            ab = 0.0;
            rt2r = 0.0;
            ba = 0.0;
          } else {
            tst /= s;
            htmp1 /= s;
            aa /= s;
            ba /= s;
            ab = (tst + ba) / 2.0;
            tst = (tst - ab) * (ba - ab) - aa * htmp1;
            htmp1 = std::sqrt(std::abs(tst));
            if (tst >= 0.0) {
              rt1r = ab * s;
              rt2r = rt1r;
              ab = htmp1 * s;
              ba = -ab;
            } else {
              rt1r = ab + htmp1;
              rt2r = ab - htmp1;
              if (std::abs(rt1r - ba) <= std::abs(rt2r - ba)) {
                rt1r *= s;
                rt2r = rt1r;
              } else {
                rt2r *= s;
                rt1r = rt2r;
              }

              ab = 0.0;
              ba = 0.0;
            }
          }

          m = i - 1;
          exitg3 = false;
          while ((!exitg3) && (m >= k + 1)) {
            s = (std::abs(h->data[(m + h->size[0] * (m - 1)) - 1] - rt2r) + std::
                 abs(ba)) + std::abs(h->data[m + h->size[0] * (m - 1)]);
            tst = h->data[m + h->size[0] * (m - 1)] / s;
            v[0] = (tst * h->data[(m + h->size[0] * m) - 1] + (h->data[(m +
                      h->size[0] * (m - 1)) - 1] - rt1r) * ((h->data[(m +
                       h->size[0] * (m - 1)) - 1] - rt2r) / s)) - ab * (ba / s);
            v[1] = tst * (((h->data[(m + h->size[0] * (m - 1)) - 1] + h->data[m
                            + h->size[0] * m]) - rt1r) - rt2r);
            v[2] = tst * h->data[(m + h->size[0] * m) + 1];
            s = (std::abs(v[0]) + std::abs(v[1])) + std::abs(v[2]);
            tst = v[0] / s;
            v[0] /= s;
            htmp1 = v[1] / s;
            v[1] /= s;
            aa = v[2] / s;
            v[2] /= s;
            if ((m == k + 1) || (std::abs(h->data[(m + h->size[0] * (m - 2)) - 1])
                                 * (std::abs(htmp1) + std::abs(aa)) <=
                                 2.2204460492503131E-16 * std::abs(tst) * ((std::
                   abs(h->data[(m + h->size[0] * (m - 2)) - 2]) + std::abs
                   (h->data[(m + h->size[0] * (m - 1)) - 1])) + std::abs(h->
                   data[m + h->size[0] * m])))) {
              exitg3 = true;
            } else {
              m--;
            }
          }

          for (b_k = m; b_k <= i; b_k++) {
            nr = (i - b_k) + 2;
            if (3 < nr) {
              nr = 3;
            }

            if (b_k > m) {
              hoffset = b_k + ldh * (b_k - 2);
              for (j = -1; j + 2 <= nr; j++) {
                v[j + 1] = h->data[j + hoffset];
              }
            }

            tst = v[0];
            rt2r = b_xzlarfg(nr, &tst, v);
            v[0] = tst;
            if (b_k > m) {
              h->data[(b_k + h->size[0] * (b_k - 2)) - 1] = tst;
              h->data[b_k + h->size[0] * (b_k - 2)] = 0.0;
              if (b_k < i) {
                h->data[(b_k + h->size[0] * (b_k - 2)) + 1] = 0.0;
              }
            } else {
              if (m > k + 1) {
                h->data[(b_k + h->size[0] * (b_k - 2)) - 1] *= 1.0 - rt2r;
              }
            }

            tst = v[1];
            htmp1 = rt2r * v[1];
            if (nr == 3) {
              ab = v[2];
              ba = rt2r * v[2];
              for (j = b_k - 1; j + 1 <= n; j++) {
                aa = (h->data[(b_k + h->size[0] * j) - 1] + tst * h->data[b_k +
                      h->size[0] * j]) + ab * h->data[(b_k + h->size[0] * j) + 1];
                h->data[(b_k + h->size[0] * j) - 1] -= aa * rt2r;
                h->data[b_k + h->size[0] * j] -= aa * htmp1;
                h->data[(b_k + h->size[0] * j) + 1] -= aa * ba;
              }

              if (b_k + 3 < i + 1) {
                nr = b_k;
              } else {
                nr = i - 2;
              }

              for (j = 0; j + 1 <= nr + 3; j++) {
                aa = (h->data[j + h->size[0] * (b_k - 1)] + tst * h->data[j +
                      h->size[0] * b_k]) + ab * h->data[j + h->size[0] * (b_k +
                  1)];
                h->data[j + h->size[0] * (b_k - 1)] -= aa * rt2r;
                h->data[j + h->size[0] * b_k] -= aa * htmp1;
                h->data[j + h->size[0] * (b_k + 1)] -= aa * ba;
              }

              for (j = 0; j + 1 <= n; j++) {
                aa = (z->data[j + z->size[0] * (b_k - 1)] + tst * z->data[j +
                      z->size[0] * b_k]) + ab * z->data[j + z->size[0] * (b_k +
                  1)];
                z->data[j + z->size[0] * (b_k - 1)] -= aa * rt2r;
                z->data[j + z->size[0] * b_k] -= aa * htmp1;
                z->data[j + z->size[0] * (b_k + 1)] -= aa * ba;
              }
            } else {
              if (nr == 2) {
                for (j = b_k - 1; j + 1 <= n; j++) {
                  aa = h->data[(b_k + h->size[0] * j) - 1] + tst * h->data[b_k +
                    h->size[0] * j];
                  h->data[(b_k + h->size[0] * j) - 1] -= aa * rt2r;
                  h->data[b_k + h->size[0] * j] -= aa * htmp1;
                }

                for (j = 0; j + 1 <= i + 1; j++) {
                  aa = h->data[j + h->size[0] * (b_k - 1)] + tst * h->data[j +
                    h->size[0] * b_k];
                  h->data[j + h->size[0] * (b_k - 1)] -= aa * rt2r;
                  h->data[j + h->size[0] * b_k] -= aa * htmp1;
                }

                for (j = 0; j + 1 <= n; j++) {
                  aa = z->data[j + z->size[0] * (b_k - 1)] + tst * z->data[j +
                    z->size[0] * b_k];
                  z->data[j + z->size[0] * (b_k - 1)] -= aa * rt2r;
                  z->data[j + z->size[0] * b_k] -= aa * htmp1;
                }
              }
            }
          }

          its++;
        }
      }

      if (!goto150) {
        info = i + 1;
        exitg1 = true;
      } else {
        if ((L != i + 1) && (L == i)) {
          tst = h->data[(i + h->size[0] * i) - 1];
          htmp1 = h->data[i + h->size[0] * (i - 1)];
          aa = h->data[i + h->size[0] * i];
          xdlanv2(&h->data[(i + h->size[0] * (i - 1)) - 1], &tst, &htmp1, &aa,
                  &ab, &ba, &rt2r, &rt1r, &s, &sn);
          h->data[(i + h->size[0] * i) - 1] = tst;
          h->data[i + h->size[0] * (i - 1)] = htmp1;
          h->data[i + h->size[0] * i] = aa;
          if (n > i + 1) {
            xrot((n - i) - 1, h, i + (i + 1) * ldh, ldh, (i + (i + 1) * ldh) + 1,
                 ldh, s, sn);
          }

          b_xrot(i - 1, h, 1 + (i - 1) * ldh, 1 + i * ldh, s, sn);
          b_xrot(n, z, 1 + (i - 1) * ldz, 1 + i * ldz, s, sn);
        }

        i = L - 2;
      }
    }
  }

  return info;
}

/* End of code generation (xdhseqr.cpp) */
