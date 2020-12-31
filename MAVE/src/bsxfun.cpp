/*
 * Student License - for use by students to meet course requirements and
 * perform academic research at degree granting institutions only.  Not
 * for government, commercial, or other organizational use.
 *
 * bsxfun.cpp
 *
 * Code generation for function 'bsxfun'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "MAVEfast.h"
#include "bsxfun.h"
#include "MAVEfast_emxutil.h"

/* Function Definitions */
void b_bsxfun(const emxArray_real_T *a, const emxArray_real_T *b,
              emxArray_real_T *c)
{
  int csz_idx_0;
  int csz_idx_1;
  int i4;
  emxArray_real_T *bv;
  unsigned int b_idx_0;
  int asub;
  int ak;
  int nc1;
  int ck;
  emxArray_real_T *cv;
  double b_a;
  csz_idx_0 = b->size[0];
  csz_idx_1 = a->size[1];
  i4 = c->size[0] * c->size[1];
  c->size[0] = csz_idx_0;
  c->size[1] = csz_idx_1;
  emxEnsureCapacity((emxArray__common *)c, i4, sizeof(double));
  if (!((c->size[0] == 0) || (c->size[1] == 0))) {
    emxInit_real_T(&bv, 1);
    b_idx_0 = (unsigned int)b->size[0];
    i4 = bv->size[0];
    bv->size[0] = (int)b_idx_0;
    emxEnsureCapacity((emxArray__common *)bv, i4, sizeof(double));
    asub = 1;
    ak = 0;
    nc1 = c->size[0];
    i4 = c->size[0] * c->size[1] - c->size[0];
    ck = 0;
    emxInit_real_T(&cv, 1);
    while (ck <= i4) {
      for (csz_idx_0 = 0; csz_idx_0 + 1 <= b->size[0]; csz_idx_0++) {
        bv->data[csz_idx_0] = b->data[csz_idx_0];
      }

      b_a = a->data[ak];
      csz_idx_0 = cv->size[0];
      cv->size[0] = bv->size[0];
      emxEnsureCapacity((emxArray__common *)cv, csz_idx_0, sizeof(double));
      csz_idx_1 = bv->size[0];
      for (csz_idx_0 = 0; csz_idx_0 < csz_idx_1; csz_idx_0++) {
        cv->data[csz_idx_0] = b_a + bv->data[csz_idx_0];
      }

      for (csz_idx_0 = 0; csz_idx_0 + 1 <= nc1; csz_idx_0++) {
        c->data[ck + csz_idx_0] = cv->data[csz_idx_0];
      }

      if (asub < a->size[1]) {
        ak++;
        asub++;
      } else {
        asub = 1;
      }

      ck += nc1;
    }

    emxFree_real_T(&cv);
    emxFree_real_T(&bv);
  }
}

void bsxfun(const emxArray_real_T *a, const emxArray_real_T *b, emxArray_real_T *
            c)
{
  int csz_idx_0;
  int csz_idx_1;
  int i2;
  emxArray_real_T *av;
  unsigned int a_idx_0;
  int bsub;
  int bk;
  int nc1;
  int ck;
  emxArray_real_T *cv;
  double b_b;
  csz_idx_0 = a->size[0];
  csz_idx_1 = b->size[1];
  i2 = c->size[0] * c->size[1];
  c->size[0] = csz_idx_0;
  c->size[1] = csz_idx_1;
  emxEnsureCapacity((emxArray__common *)c, i2, sizeof(double));
  if (!((c->size[0] == 0) || (c->size[1] == 0))) {
    emxInit_real_T(&av, 1);
    a_idx_0 = (unsigned int)a->size[0];
    i2 = av->size[0];
    av->size[0] = (int)a_idx_0;
    emxEnsureCapacity((emxArray__common *)av, i2, sizeof(double));
    bsub = 1;
    bk = 0;
    nc1 = c->size[0];
    i2 = c->size[0] * c->size[1] - c->size[0];
    ck = 0;
    emxInit_real_T(&cv, 1);
    while (ck <= i2) {
      for (csz_idx_0 = 0; csz_idx_0 + 1 <= a->size[0]; csz_idx_0++) {
        av->data[csz_idx_0] = a->data[csz_idx_0];
      }

      b_b = b->data[bk];
      csz_idx_0 = cv->size[0];
      cv->size[0] = av->size[0];
      emxEnsureCapacity((emxArray__common *)cv, csz_idx_0, sizeof(double));
      csz_idx_1 = av->size[0];
      for (csz_idx_0 = 0; csz_idx_0 < csz_idx_1; csz_idx_0++) {
        cv->data[csz_idx_0] = av->data[csz_idx_0] + b_b;
      }

      for (csz_idx_0 = 0; csz_idx_0 + 1 <= nc1; csz_idx_0++) {
        c->data[ck + csz_idx_0] = cv->data[csz_idx_0];
      }

      if (bsub < b->size[1]) {
        bk++;
        bsub++;
      } else {
        bsub = 1;
      }

      ck += nc1;
    }

    emxFree_real_T(&cv);
    emxFree_real_T(&av);
  }
}

/* End of code generation (bsxfun.cpp) */
