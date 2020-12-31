/*
 * Student License - for use by students to meet course requirements and
 * perform academic research at degree granting institutions only.  Not
 * for government, commercial, or other organizational use.
 *
 * MAVEfast_emxAPI.cpp
 *
 * Code generation for function 'MAVEfast_emxAPI'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "MAVEfast.h"
#include "MAVEfast_emxAPI.h"
#include "MAVEfast_emxutil.h"

/* Function Definitions */
emxArray_char_T *emxCreateND_char_T(int numDimensions, int *size)
{
  emxArray_char_T *emx;
  int numEl;
  int i;
  emxInit_char_T(&emx, numDimensions);
  numEl = 1;
  for (i = 0; i < numDimensions; i++) {
    numEl *= size[i];
    emx->size[i] = size[i];
  }

  emx->data = (char *)calloc((unsigned int)numEl, sizeof(char));
  emx->numDimensions = numDimensions;
  emx->allocatedSize = numEl;
  return emx;
}

emxArray_real_T *emxCreateND_real_T(int numDimensions, int *size)
{
  emxArray_real_T *emx;
  int numEl;
  int i;
  emxInit_real_T2(&emx, numDimensions);
  numEl = 1;
  for (i = 0; i < numDimensions; i++) {
    numEl *= size[i];
    emx->size[i] = size[i];
  }

  emx->data = (double *)calloc((unsigned int)numEl, sizeof(double));
  emx->numDimensions = numDimensions;
  emx->allocatedSize = numEl;
  return emx;
}

emxArray_char_T *emxCreateWrapperND_char_T(char *data, int numDimensions, int
  *size)
{
  emxArray_char_T *emx;
  int numEl;
  int i;
  emxInit_char_T(&emx, numDimensions);
  numEl = 1;
  for (i = 0; i < numDimensions; i++) {
    numEl *= size[i];
    emx->size[i] = size[i];
  }

  emx->data = data;
  emx->numDimensions = numDimensions;
  emx->allocatedSize = numEl;
  emx->canFreeData = false;
  return emx;
}

emxArray_real_T *emxCreateWrapperND_real_T(double *data, int numDimensions, int *
  size)
{
  emxArray_real_T *emx;
  int numEl;
  int i;
  emxInit_real_T2(&emx, numDimensions);
  numEl = 1;
  for (i = 0; i < numDimensions; i++) {
    numEl *= size[i];
    emx->size[i] = size[i];
  }

  emx->data = data;
  emx->numDimensions = numDimensions;
  emx->allocatedSize = numEl;
  emx->canFreeData = false;
  return emx;
}

emxArray_char_T *emxCreateWrapper_char_T(char *data, int rows, int cols)
{
  emxArray_char_T *emx;
  int size[2];
  int numEl;
  int i;
  size[0] = rows;
  size[1] = cols;
  emxInit_char_T(&emx, 2);
  numEl = 1;
  for (i = 0; i < 2; i++) {
    numEl *= size[i];
    emx->size[i] = size[i];
  }

  emx->data = data;
  emx->numDimensions = 2;
  emx->allocatedSize = numEl;
  emx->canFreeData = false;
  return emx;
}

emxArray_real_T *emxCreateWrapper_real_T(double *data, int rows, int cols)
{
  emxArray_real_T *emx;
  int size[2];
  int numEl;
  int i;
  size[0] = rows;
  size[1] = cols;
  emxInit_real_T2(&emx, 2);
  numEl = 1;
  for (i = 0; i < 2; i++) {
    numEl *= size[i];
    emx->size[i] = size[i];
  }

  emx->data = data;
  emx->numDimensions = 2;
  emx->allocatedSize = numEl;
  emx->canFreeData = false;
  return emx;
}

emxArray_char_T *emxCreate_char_T(int rows, int cols)
{
  emxArray_char_T *emx;
  int size[2];
  int numEl;
  int i;
  size[0] = rows;
  size[1] = cols;
  emxInit_char_T(&emx, 2);
  numEl = 1;
  for (i = 0; i < 2; i++) {
    numEl *= size[i];
    emx->size[i] = size[i];
  }

  emx->data = (char *)calloc((unsigned int)numEl, sizeof(char));
  emx->numDimensions = 2;
  emx->allocatedSize = numEl;
  return emx;
}

emxArray_real_T *emxCreate_real_T(int rows, int cols)
{
  emxArray_real_T *emx;
  int size[2];
  int numEl;
  int i;
  size[0] = rows;
  size[1] = cols;
  emxInit_real_T2(&emx, 2);
  numEl = 1;
  for (i = 0; i < 2; i++) {
    numEl *= size[i];
    emx->size[i] = size[i];
  }

  emx->data = (double *)calloc((unsigned int)numEl, sizeof(double));
  emx->numDimensions = 2;
  emx->allocatedSize = numEl;
  return emx;
}

void emxDestroyArray_char_T(emxArray_char_T *emxArray)
{
  emxFree_char_T(&emxArray);
}

void emxDestroyArray_real_T(emxArray_real_T *emxArray)
{
  emxFree_real_T(&emxArray);
}

void emxInitArray_char_T(emxArray_char_T **pEmxArray, int numDimensions)
{
  emxInit_char_T(pEmxArray, numDimensions);
}

void emxInitArray_real_T(emxArray_real_T **pEmxArray, int numDimensions)
{
  emxInit_real_T2(pEmxArray, numDimensions);
}

/* End of code generation (MAVEfast_emxAPI.cpp) */
