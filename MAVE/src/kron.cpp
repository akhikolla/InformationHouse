/*
 * Student License - for use by students to meet course requirements and
 * perform academic research at degree granting institutions only.  Not
 * for government, commercial, or other organizational use.
 *
 * kron.cpp
 *
 * Code generation for function 'kron'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "MAVEfast.h"
#include "kron.h"
#include "MAVEfast_emxutil.h"

/* Function Definitions */
void kron(const emxArray_real_T *A, const emxArray_real_T *B, emxArray_real_T *K)
{
  int kidx;
  int unnamed_idx_1;
  int j2;
  int i1;
  int i2;
  kidx = A->size[0] * B->size[0];
  unnamed_idx_1 = A->size[1] * B->size[1];
  j2 = K->size[0] * K->size[1];
  K->size[0] = kidx;
  K->size[1] = unnamed_idx_1;
  emxEnsureCapacity((emxArray__common *)K, j2, sizeof(double));
  kidx = -1;
  for (unnamed_idx_1 = 1; unnamed_idx_1 <= A->size[1]; unnamed_idx_1++) {
    for (j2 = 1; j2 <= B->size[1]; j2++) {
      for (i1 = 1; i1 <= A->size[0]; i1++) {
        for (i2 = 1; i2 <= B->size[0]; i2++) {
          kidx++;
          K->data[kidx] = A->data[(i1 + A->size[0] * (unnamed_idx_1 - 1)) - 1] *
            B->data[(i2 + B->size[0] * (j2 - 1)) - 1];
        }
      }
    }
  }
}

/* End of code generation (kron.cpp) */
