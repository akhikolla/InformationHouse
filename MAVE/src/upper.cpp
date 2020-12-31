/*
 * Student License - for use by students to meet course requirements and
 * perform academic research at degree granting institutions only.  Not
 * for government, commercial, or other organizational use.
 *
 * upper.cpp
 *
 * Code generation for function 'upper'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "MAVEfast.h"
#include "upper.h"
#include "MAVEfast_emxutil.h"

/* Function Definitions */
void upper(const emxArray_char_T *x, emxArray_char_T *y)
{
  int i8;
  int k;
  static const char cv0[128] = { '\x00', '\x01', '\x02', '\x03', '\x04', '\x05',
    '\x06', '\x07', '\x08', '\x09', '\x0a', '\x0b', '\x0c', '\x0d', '\x0e',
    '\x0f', '\x10', '\x11', '\x12', '\x13', '\x14', '\x15', '\x16', '\x17',
    '\x18', '\x19', '\x1a', '\x1b', '\x1c', '\x1d', '\x1e', '\x1f', ' ', '!',
    '\"', '#', '$', '%', '&', '\'', '(', ')', '*', '+', ',', '-', '.', '/', '0',
    '1', '2', '3', '4', '5', '6', '7', '8', '9', ':', ';', '<', '=', '>', '?',
    '@', 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N',
    'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z', '[', '\\', ']',
    '^', '_', '`', 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L',
    'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z', '{',
    '|', '}', '~', '\x7f' };

  i8 = y->size[0] * y->size[1];
  y->size[0] = x->size[0];
  y->size[1] = x->size[1];
  emxEnsureCapacity((emxArray__common *)y, i8, sizeof(char));
  i8 = x->size[0] * x->size[1];
  for (k = 0; k < i8; k++) {
    y->data[k] = cv0[(unsigned char)x->data[k] & 127];
  }
}

/* End of code generation (upper.cpp) */
