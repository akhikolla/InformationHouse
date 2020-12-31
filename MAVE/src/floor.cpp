/*
 * Student License - for use by students to meet course requirements and
 * perform academic research at degree granting institutions only.  Not
 * for government, commercial, or other organizational use.
 *
 * floor.cpp
 *
 * Code generation for function 'floor'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "MAVEfast.h"
#include "floor.h"

/* Function Definitions */
void b_floor(double x_data[], int x_size[2])
{
  int nx;
  int k;
  nx = x_size[1];
  for (k = 0; k + 1 <= nx; k++) {
    x_data[k] = std::floor(x_data[k]);
  }
}

/* End of code generation (floor.cpp) */
