/*
 * Student License - for use by students to meet course requirements and
 * perform academic research at degree granting institutions only.  Not
 * for government, commercial, or other organizational use.
 *
 * strcmp.cpp
 *
 * Code generation for function 'strcmp'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "MAVEfast.h"
#include "strcmp.h"

/* Function Definitions */
boolean_T b_strcmp(const emxArray_char_T *a)
{
  boolean_T b_bool;
  int kstr;
  int exitg1;
  static const char cv1[8] = { 'M', 'E', 'A', 'N', 'M', 'A', 'V', 'E' };

  b_bool = false;
  if ((a->size[0] != 1) || (a->size[1] != 8)) {
  } else {
    kstr = 0;
    do {
      exitg1 = 0;
      if (kstr + 1 < 9) {
        if (a->data[kstr] != cv1[kstr]) {
          exitg1 = 1;
        } else {
          kstr++;
        }
      } else {
        b_bool = true;
        exitg1 = 1;
      }
    } while (exitg1 == 0);
  }

  return b_bool;
}

boolean_T c_strcmp(const emxArray_char_T *a)
{
  boolean_T b_bool;
  int kstr;
  int exitg1;
  static const char cv2[7] = { 'M', 'E', 'A', 'N', 'O', 'P', 'G' };

  b_bool = false;
  if ((a->size[0] != 1) || (a->size[1] != 7)) {
  } else {
    kstr = 0;
    do {
      exitg1 = 0;
      if (kstr + 1 < 8) {
        if (a->data[kstr] != cv2[kstr]) {
          exitg1 = 1;
        } else {
          kstr++;
        }
      } else {
        b_bool = true;
        exitg1 = 1;
      }
    } while (exitg1 == 0);
  }

  return b_bool;
}

boolean_T d_strcmp(const emxArray_char_T *a)
{
  boolean_T b_bool;
  int kstr;
  int exitg1;
  static const char cv3[4] = { 'K', 'S', 'I', 'R' };

  b_bool = false;
  if ((a->size[0] != 1) || (a->size[1] != 4)) {
  } else {
    kstr = 0;
    do {
      exitg1 = 0;
      if (kstr + 1 < 5) {
        if (a->data[kstr] != cv3[kstr]) {
          exitg1 = 1;
        } else {
          kstr++;
        }
      } else {
        b_bool = true;
        exitg1 = 1;
      }
    } while (exitg1 == 0);
  }

  return b_bool;
}

boolean_T e_strcmp(const emxArray_char_T *a)
{
  boolean_T b_bool;
  int kstr;
  int exitg1;
  static const char cv4[5] = { 'C', 'S', 'O', 'P', 'G' };

  b_bool = false;
  if ((a->size[0] != 1) || (a->size[1] != 5)) {
  } else {
    kstr = 0;
    do {
      exitg1 = 0;
      if (kstr + 1 < 6) {
        if (a->data[kstr] != cv4[kstr]) {
          exitg1 = 1;
        } else {
          kstr++;
        }
      } else {
        b_bool = true;
        exitg1 = 1;
      }
    } while (exitg1 == 0);
  }

  return b_bool;
}

/* End of code generation (strcmp.cpp) */
