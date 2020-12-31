/*
 * Student License - for use by students to meet course requirements and
 * perform academic research at degree granting institutions only.  Not
 * for government, commercial, or other organizational use.
 *
 * sortIdx.cpp
 *
 * Code generation for function 'sortIdx'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "MAVEfast.h"
#include "sortIdx.h"
#include "MAVEfast_emxutil.h"
#include "mergesort.h"

/* Function Declarations */
static void b_merge(emxArray_int32_T *idx, emxArray_real_T *x, int offset, int
                    np, int nq, emxArray_int32_T *iwork, emxArray_real_T *xwork);
static void b_merge_block(emxArray_int32_T *idx, emxArray_real_T *x, int offset,
  int n, int preSortLevel, emxArray_int32_T *iwork, emxArray_real_T *xwork);
static void merge(emxArray_int32_T *idx, emxArray_real_T *x, int offset, int np,
                  int nq, emxArray_int32_T *iwork, emxArray_real_T *xwork);
static void merge_block(emxArray_int32_T *idx, emxArray_real_T *x, int offset,
  int n, int preSortLevel, emxArray_int32_T *iwork, emxArray_real_T *xwork);

/* Function Definitions */
static void b_merge(emxArray_int32_T *idx, emxArray_real_T *x, int offset, int
                    np, int nq, emxArray_int32_T *iwork, emxArray_real_T *xwork)
{
  int n;
  int qend;
  int p;
  int iout;
  int exitg1;
  if (nq != 0) {
    n = np + nq;
    for (qend = 0; qend + 1 <= n; qend++) {
      iwork->data[qend] = idx->data[offset + qend];
      xwork->data[qend] = x->data[offset + qend];
    }

    p = 0;
    n = np;
    qend = np + nq;
    iout = offset - 1;
    do {
      exitg1 = 0;
      iout++;
      if (xwork->data[p] <= xwork->data[n]) {
        idx->data[iout] = iwork->data[p];
        x->data[iout] = xwork->data[p];
        if (p + 1 < np) {
          p++;
        } else {
          exitg1 = 1;
        }
      } else {
        idx->data[iout] = iwork->data[n];
        x->data[iout] = xwork->data[n];
        if (n + 1 < qend) {
          n++;
        } else {
          n = (iout - p) + 1;
          while (p + 1 <= np) {
            idx->data[n + p] = iwork->data[p];
            x->data[n + p] = xwork->data[p];
            p++;
          }

          exitg1 = 1;
        }
      }
    } while (exitg1 == 0);
  }
}

static void b_merge_block(emxArray_int32_T *idx, emxArray_real_T *x, int offset,
  int n, int preSortLevel, emxArray_int32_T *iwork, emxArray_real_T *xwork)
{
  int nPairs;
  int bLen;
  int tailOffset;
  int nTail;
  nPairs = n >> preSortLevel;
  bLen = 1 << preSortLevel;
  while (nPairs > 1) {
    if ((nPairs & 1) != 0) {
      nPairs--;
      tailOffset = bLen * nPairs;
      nTail = n - tailOffset;
      if (nTail > bLen) {
        b_merge(idx, x, offset + tailOffset, bLen, nTail - bLen, iwork, xwork);
      }
    }

    tailOffset = bLen << 1;
    nPairs >>= 1;
    for (nTail = 1; nTail <= nPairs; nTail++) {
      b_merge(idx, x, offset + (nTail - 1) * tailOffset, bLen, bLen, iwork,
              xwork);
    }

    bLen = tailOffset;
  }

  if (n > bLen) {
    b_merge(idx, x, offset, bLen, n - bLen, iwork, xwork);
  }
}

static void merge(emxArray_int32_T *idx, emxArray_real_T *x, int offset, int np,
                  int nq, emxArray_int32_T *iwork, emxArray_real_T *xwork)
{
  int n;
  int qend;
  int p;
  int iout;
  int exitg1;
  if (nq != 0) {
    n = np + nq;
    for (qend = 0; qend + 1 <= n; qend++) {
      iwork->data[qend] = idx->data[offset + qend];
      xwork->data[qend] = x->data[offset + qend];
    }

    p = 0;
    n = np;
    qend = np + nq;
    iout = offset - 1;
    do {
      exitg1 = 0;
      iout++;
      if (xwork->data[p] >= xwork->data[n]) {
        idx->data[iout] = iwork->data[p];
        x->data[iout] = xwork->data[p];
        if (p + 1 < np) {
          p++;
        } else {
          exitg1 = 1;
        }
      } else {
        idx->data[iout] = iwork->data[n];
        x->data[iout] = xwork->data[n];
        if (n + 1 < qend) {
          n++;
        } else {
          n = (iout - p) + 1;
          while (p + 1 <= np) {
            idx->data[n + p] = iwork->data[p];
            x->data[n + p] = xwork->data[p];
            p++;
          }

          exitg1 = 1;
        }
      }
    } while (exitg1 == 0);
  }
}

static void merge_block(emxArray_int32_T *idx, emxArray_real_T *x, int offset,
  int n, int preSortLevel, emxArray_int32_T *iwork, emxArray_real_T *xwork)
{
  int nPairs;
  int bLen;
  int tailOffset;
  int nTail;
  nPairs = n >> preSortLevel;
  bLen = 1 << preSortLevel;
  while (nPairs > 1) {
    if ((nPairs & 1) != 0) {
      nPairs--;
      tailOffset = bLen * nPairs;
      nTail = n - tailOffset;
      if (nTail > bLen) {
        merge(idx, x, offset + tailOffset, bLen, nTail - bLen, iwork, xwork);
      }
    }

    tailOffset = bLen << 1;
    nPairs >>= 1;
    for (nTail = 1; nTail <= nPairs; nTail++) {
      merge(idx, x, offset + (nTail - 1) * tailOffset, bLen, bLen, iwork, xwork);
    }

    bLen = tailOffset;
  }

  if (n > bLen) {
    merge(idx, x, offset, bLen, n - bLen, iwork, xwork);
  }
}

void b_sortIdx(emxArray_real_T *x, emxArray_int32_T *idx)
{
  int ib;
  int i;
  int n;
  int b_n;
  emxArray_int32_T *iwork;
  double x4[4];
  int idx4[4];
  emxArray_real_T *xwork;
  int nNaNs;
  int k;
  int wOffset;
  signed char perm[4];
  int nNonNaN;
  int i3;
  int i4;
  int nBlocks;
  int bLen2;
  int nPairs;
  int b_iwork[256];
  double b_xwork[256];
  int exitg1;
  ib = x->size[0];
  i = idx->size[0];
  idx->size[0] = ib;
  emxEnsureCapacity((emxArray__common *)idx, i, sizeof(int));
  for (i = 0; i < ib; i++) {
    idx->data[i] = 0;
  }

  n = x->size[0];
  b_n = x->size[0];
  for (i = 0; i < 4; i++) {
    x4[i] = 0.0;
    idx4[i] = 0;
  }

  emxInit_int32_T(&iwork, 1);
  i = iwork->size[0];
  iwork->size[0] = ib;
  emxEnsureCapacity((emxArray__common *)iwork, i, sizeof(int));
  ib = iwork->size[0];
  i = iwork->size[0];
  iwork->size[0] = ib;
  emxEnsureCapacity((emxArray__common *)iwork, i, sizeof(int));
  for (i = 0; i < ib; i++) {
    iwork->data[i] = 0;
  }

  emxInit_real_T(&xwork, 1);
  ib = x->size[0];
  i = xwork->size[0];
  xwork->size[0] = ib;
  emxEnsureCapacity((emxArray__common *)xwork, i, sizeof(double));
  ib = xwork->size[0];
  i = xwork->size[0];
  xwork->size[0] = ib;
  emxEnsureCapacity((emxArray__common *)xwork, i, sizeof(double));
  for (i = 0; i < ib; i++) {
    xwork->data[i] = 0.0;
  }

  nNaNs = 0;
  ib = 0;
  for (k = 0; k + 1 <= b_n; k++) {
    if (rtIsNaN(x->data[k])) {
      idx->data[(b_n - nNaNs) - 1] = k + 1;
      xwork->data[(b_n - nNaNs) - 1] = x->data[k];
      nNaNs++;
    } else {
      ib++;
      idx4[ib - 1] = k + 1;
      x4[ib - 1] = x->data[k];
      if (ib == 4) {
        i = k - nNaNs;
        if (x4[0] >= x4[1]) {
          ib = 1;
          wOffset = 2;
        } else {
          ib = 2;
          wOffset = 1;
        }

        if (x4[2] >= x4[3]) {
          i3 = 3;
          i4 = 4;
        } else {
          i3 = 4;
          i4 = 3;
        }

        if (x4[ib - 1] >= x4[i3 - 1]) {
          if (x4[wOffset - 1] >= x4[i3 - 1]) {
            perm[0] = (signed char)ib;
            perm[1] = (signed char)wOffset;
            perm[2] = (signed char)i3;
            perm[3] = (signed char)i4;
          } else if (x4[wOffset - 1] >= x4[i4 - 1]) {
            perm[0] = (signed char)ib;
            perm[1] = (signed char)i3;
            perm[2] = (signed char)wOffset;
            perm[3] = (signed char)i4;
          } else {
            perm[0] = (signed char)ib;
            perm[1] = (signed char)i3;
            perm[2] = (signed char)i4;
            perm[3] = (signed char)wOffset;
          }
        } else if (x4[ib - 1] >= x4[i4 - 1]) {
          if (x4[wOffset - 1] >= x4[i4 - 1]) {
            perm[0] = (signed char)i3;
            perm[1] = (signed char)ib;
            perm[2] = (signed char)wOffset;
            perm[3] = (signed char)i4;
          } else {
            perm[0] = (signed char)i3;
            perm[1] = (signed char)ib;
            perm[2] = (signed char)i4;
            perm[3] = (signed char)wOffset;
          }
        } else {
          perm[0] = (signed char)i3;
          perm[1] = (signed char)i4;
          perm[2] = (signed char)ib;
          perm[3] = (signed char)wOffset;
        }

        idx->data[i - 3] = idx4[perm[0] - 1];
        idx->data[i - 2] = idx4[perm[1] - 1];
        idx->data[i - 1] = idx4[perm[2] - 1];
        idx->data[i] = idx4[perm[3] - 1];
        x->data[i - 3] = x4[perm[0] - 1];
        x->data[i - 2] = x4[perm[1] - 1];
        x->data[i - 1] = x4[perm[2] - 1];
        x->data[i] = x4[perm[3] - 1];
        ib = 0;
      }
    }
  }

  wOffset = (b_n - nNaNs) - 1;
  if (ib > 0) {
    for (i = 0; i < 4; i++) {
      perm[i] = 0;
    }

    if (ib == 1) {
      perm[0] = 1;
    } else if (ib == 2) {
      if (x4[0] >= x4[1]) {
        perm[0] = 1;
        perm[1] = 2;
      } else {
        perm[0] = 2;
        perm[1] = 1;
      }
    } else if (x4[0] >= x4[1]) {
      if (x4[1] >= x4[2]) {
        perm[0] = 1;
        perm[1] = 2;
        perm[2] = 3;
      } else if (x4[0] >= x4[2]) {
        perm[0] = 1;
        perm[1] = 3;
        perm[2] = 2;
      } else {
        perm[0] = 3;
        perm[1] = 1;
        perm[2] = 2;
      }
    } else if (x4[0] >= x4[2]) {
      perm[0] = 2;
      perm[1] = 1;
      perm[2] = 3;
    } else if (x4[1] >= x4[2]) {
      perm[0] = 2;
      perm[1] = 3;
      perm[2] = 1;
    } else {
      perm[0] = 3;
      perm[1] = 2;
      perm[2] = 1;
    }

    for (k = 1; k <= ib; k++) {
      idx->data[(wOffset - ib) + k] = idx4[perm[k - 1] - 1];
      x->data[(wOffset - ib) + k] = x4[perm[k - 1] - 1];
    }
  }

  i = nNaNs >> 1;
  for (k = 1; k <= i; k++) {
    ib = idx->data[wOffset + k];
    idx->data[wOffset + k] = idx->data[b_n - k];
    idx->data[b_n - k] = ib;
    x->data[wOffset + k] = xwork->data[b_n - k];
    x->data[b_n - k] = xwork->data[wOffset + k];
  }

  if ((nNaNs & 1) != 0) {
    x->data[(wOffset + i) + 1] = xwork->data[(wOffset + i) + 1];
  }

  nNonNaN = n - nNaNs;
  ib = 2;
  if (nNonNaN > 1) {
    if (n >= 256) {
      nBlocks = nNonNaN >> 8;
      if (nBlocks > 0) {
        for (i3 = 1; i3 <= nBlocks; i3++) {
          i4 = ((i3 - 1) << 8) - 1;
          for (b_n = 0; b_n < 6; b_n++) {
            n = 1 << (b_n + 2);
            bLen2 = n << 1;
            nPairs = 256 >> (b_n + 3);
            for (k = 1; k <= nPairs; k++) {
              ib = i4 + (k - 1) * bLen2;
              for (i = 1; i <= bLen2; i++) {
                b_iwork[i - 1] = idx->data[ib + i];
                b_xwork[i - 1] = x->data[ib + i];
              }

              wOffset = 0;
              i = n;
              do {
                exitg1 = 0;
                ib++;
                if (b_xwork[wOffset] >= b_xwork[i]) {
                  idx->data[ib] = b_iwork[wOffset];
                  x->data[ib] = b_xwork[wOffset];
                  if (wOffset + 1 < n) {
                    wOffset++;
                  } else {
                    exitg1 = 1;
                  }
                } else {
                  idx->data[ib] = b_iwork[i];
                  x->data[ib] = b_xwork[i];
                  if (i + 1 < bLen2) {
                    i++;
                  } else {
                    i = ib - wOffset;
                    while (wOffset + 1 <= n) {
                      idx->data[(i + wOffset) + 1] = b_iwork[wOffset];
                      x->data[(i + wOffset) + 1] = b_xwork[wOffset];
                      wOffset++;
                    }

                    exitg1 = 1;
                  }
                }
              } while (exitg1 == 0);
            }
          }
        }

        ib = nBlocks << 8;
        i = nNonNaN - ib;
        if (i > 0) {
          merge_block(idx, x, ib, i, 2, iwork, xwork);
        }

        ib = 8;
      }
    }

    merge_block(idx, x, 0, nNonNaN, ib, iwork, xwork);
  }

  if ((nNaNs > 0) && (nNonNaN > 0)) {
    for (k = 0; k + 1 <= nNaNs; k++) {
      xwork->data[k] = x->data[nNonNaN + k];
      iwork->data[k] = idx->data[nNonNaN + k];
    }

    for (k = nNonNaN - 1; k + 1 > 0; k--) {
      x->data[nNaNs + k] = x->data[k];
      idx->data[nNaNs + k] = idx->data[k];
    }

    for (k = 0; k + 1 <= nNaNs; k++) {
      x->data[k] = xwork->data[k];
      idx->data[k] = iwork->data[k];
    }
  }

  emxFree_real_T(&xwork);
  emxFree_int32_T(&iwork);
}

void c_sortIdx(emxArray_real_T *x, emxArray_int32_T *idx)
{
  int ib;
  int m;
  int n;
  int b_n;
  int i;
  emxArray_int32_T *iwork;
  double x4[4];
  int idx4[4];
  emxArray_real_T *xwork;
  int nNaNs;
  int k;
  int wOffset;
  signed char perm[4];
  int nNonNaN;
  int i4;
  int nBlocks;
  int nPairs;
  int b_iwork[256];
  double b_xwork[256];
  int exitg1;
  ib = x->size[0];
  m = idx->size[0];
  idx->size[0] = ib;
  emxEnsureCapacity((emxArray__common *)idx, m, sizeof(int));
  for (m = 0; m < ib; m++) {
    idx->data[m] = 0;
  }

  n = x->size[0];
  b_n = x->size[0];
  for (i = 0; i < 4; i++) {
    x4[i] = 0.0;
    idx4[i] = 0;
  }

  emxInit_int32_T(&iwork, 1);
  m = iwork->size[0];
  iwork->size[0] = ib;
  emxEnsureCapacity((emxArray__common *)iwork, m, sizeof(int));
  i = iwork->size[0];
  m = iwork->size[0];
  iwork->size[0] = i;
  emxEnsureCapacity((emxArray__common *)iwork, m, sizeof(int));
  for (m = 0; m < i; m++) {
    iwork->data[m] = 0;
  }

  emxInit_real_T(&xwork, 1);
  i = x->size[0];
  m = xwork->size[0];
  xwork->size[0] = i;
  emxEnsureCapacity((emxArray__common *)xwork, m, sizeof(double));
  i = xwork->size[0];
  m = xwork->size[0];
  xwork->size[0] = i;
  emxEnsureCapacity((emxArray__common *)xwork, m, sizeof(double));
  for (m = 0; m < i; m++) {
    xwork->data[m] = 0.0;
  }

  nNaNs = 1;
  ib = 0;
  for (k = 0; k + 1 <= b_n; k++) {
    if (rtIsNaN(x->data[k])) {
      idx->data[b_n - nNaNs] = k + 1;
      xwork->data[b_n - nNaNs] = x->data[k];
      nNaNs++;
    } else {
      ib++;
      idx4[ib - 1] = k + 1;
      x4[ib - 1] = x->data[k];
      if (ib == 4) {
        i = k - nNaNs;
        if (x4[0] <= x4[1]) {
          ib = 1;
          m = 2;
        } else {
          ib = 2;
          m = 1;
        }

        if (x4[2] <= x4[3]) {
          wOffset = 3;
          i4 = 4;
        } else {
          wOffset = 4;
          i4 = 3;
        }

        if (x4[ib - 1] <= x4[wOffset - 1]) {
          if (x4[m - 1] <= x4[wOffset - 1]) {
            perm[0] = (signed char)ib;
            perm[1] = (signed char)m;
            perm[2] = (signed char)wOffset;
            perm[3] = (signed char)i4;
          } else if (x4[m - 1] <= x4[i4 - 1]) {
            perm[0] = (signed char)ib;
            perm[1] = (signed char)wOffset;
            perm[2] = (signed char)m;
            perm[3] = (signed char)i4;
          } else {
            perm[0] = (signed char)ib;
            perm[1] = (signed char)wOffset;
            perm[2] = (signed char)i4;
            perm[3] = (signed char)m;
          }
        } else if (x4[ib - 1] <= x4[i4 - 1]) {
          if (x4[m - 1] <= x4[i4 - 1]) {
            perm[0] = (signed char)wOffset;
            perm[1] = (signed char)ib;
            perm[2] = (signed char)m;
            perm[3] = (signed char)i4;
          } else {
            perm[0] = (signed char)wOffset;
            perm[1] = (signed char)ib;
            perm[2] = (signed char)i4;
            perm[3] = (signed char)m;
          }
        } else {
          perm[0] = (signed char)wOffset;
          perm[1] = (signed char)i4;
          perm[2] = (signed char)ib;
          perm[3] = (signed char)m;
        }

        idx->data[i - 2] = idx4[perm[0] - 1];
        idx->data[i - 1] = idx4[perm[1] - 1];
        idx->data[i] = idx4[perm[2] - 1];
        idx->data[i + 1] = idx4[perm[3] - 1];
        x->data[i - 2] = x4[perm[0] - 1];
        x->data[i - 1] = x4[perm[1] - 1];
        x->data[i] = x4[perm[2] - 1];
        x->data[i + 1] = x4[perm[3] - 1];
        ib = 0;
      }
    }
  }

  wOffset = b_n - nNaNs;
  if (ib > 0) {
    for (i = 0; i < 4; i++) {
      perm[i] = 0;
    }

    if (ib == 1) {
      perm[0] = 1;
    } else if (ib == 2) {
      if (x4[0] <= x4[1]) {
        perm[0] = 1;
        perm[1] = 2;
      } else {
        perm[0] = 2;
        perm[1] = 1;
      }
    } else if (x4[0] <= x4[1]) {
      if (x4[1] <= x4[2]) {
        perm[0] = 1;
        perm[1] = 2;
        perm[2] = 3;
      } else if (x4[0] <= x4[2]) {
        perm[0] = 1;
        perm[1] = 3;
        perm[2] = 2;
      } else {
        perm[0] = 3;
        perm[1] = 1;
        perm[2] = 2;
      }
    } else if (x4[0] <= x4[2]) {
      perm[0] = 2;
      perm[1] = 1;
      perm[2] = 3;
    } else if (x4[1] <= x4[2]) {
      perm[0] = 2;
      perm[1] = 3;
      perm[2] = 1;
    } else {
      perm[0] = 3;
      perm[1] = 2;
      perm[2] = 1;
    }

    for (k = 1; k <= ib; k++) {
      idx->data[(wOffset - ib) + k] = idx4[perm[k - 1] - 1];
      x->data[(wOffset - ib) + k] = x4[perm[k - 1] - 1];
    }
  }

  m = (nNaNs - 1) >> 1;
  for (k = 1; k <= m; k++) {
    i = idx->data[wOffset + k];
    idx->data[wOffset + k] = idx->data[b_n - k];
    idx->data[b_n - k] = i;
    x->data[wOffset + k] = xwork->data[b_n - k];
    x->data[b_n - k] = xwork->data[wOffset + k];
  }

  if (((nNaNs - 1) & 1) != 0) {
    x->data[(wOffset + m) + 1] = xwork->data[(wOffset + m) + 1];
  }

  nNonNaN = (n - nNaNs) + 1;
  i = 2;
  if (nNonNaN > 1) {
    if (n >= 256) {
      nBlocks = nNonNaN >> 8;
      if (nBlocks > 0) {
        for (wOffset = 1; wOffset <= nBlocks; wOffset++) {
          i4 = ((wOffset - 1) << 8) - 1;
          for (nNaNs = 0; nNaNs < 6; nNaNs++) {
            b_n = 1 << (nNaNs + 2);
            n = b_n << 1;
            nPairs = 256 >> (nNaNs + 3);
            for (k = 1; k <= nPairs; k++) {
              ib = i4 + (k - 1) * n;
              for (i = 1; i <= n; i++) {
                b_iwork[i - 1] = idx->data[ib + i];
                b_xwork[i - 1] = x->data[ib + i];
              }

              m = 0;
              i = b_n;
              do {
                exitg1 = 0;
                ib++;
                if (b_xwork[m] <= b_xwork[i]) {
                  idx->data[ib] = b_iwork[m];
                  x->data[ib] = b_xwork[m];
                  if (m + 1 < b_n) {
                    m++;
                  } else {
                    exitg1 = 1;
                  }
                } else {
                  idx->data[ib] = b_iwork[i];
                  x->data[ib] = b_xwork[i];
                  if (i + 1 < n) {
                    i++;
                  } else {
                    i = ib - m;
                    while (m + 1 <= b_n) {
                      idx->data[(i + m) + 1] = b_iwork[m];
                      x->data[(i + m) + 1] = b_xwork[m];
                      m++;
                    }

                    exitg1 = 1;
                  }
                }
              } while (exitg1 == 0);
            }
          }
        }

        i = nBlocks << 8;
        ib = nNonNaN - i;
        if (ib > 0) {
          b_merge_block(idx, x, i, ib, 2, iwork, xwork);
        }

        i = 8;
      }
    }

    b_merge_block(idx, x, 0, nNonNaN, i, iwork, xwork);
  }

  emxFree_real_T(&xwork);
  emxFree_int32_T(&iwork);
}

void sortIdx(const emxArray_real_T *x, emxArray_int32_T *idx)
{
  int i11;
  int loop_ub;
  i11 = idx->size[0] * idx->size[1];
  idx->size[0] = 1;
  idx->size[1] = x->size[1];
  emxEnsureCapacity((emxArray__common *)idx, i11, sizeof(int));
  loop_ub = x->size[1];
  for (i11 = 0; i11 < loop_ub; i11++) {
    idx->data[i11] = 0;
  }

  if (x->size[1] != 0) {
    b_mergesort(idx, x, x->size[1]);
  }
}

/* End of code generation (sortIdx.cpp) */
