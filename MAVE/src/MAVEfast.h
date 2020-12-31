/*
 * Student License - for use by students to meet course requirements and
 * perform academic research at degree granting institutions only.  Not
 * for government, commercial, or other organizational use.
 *
 * MAVEfast.h
 *
 * Code generation for function 'MAVEfast'
 *
 */

#ifndef MAVEFAST_H
#define MAVEFAST_H

/* Include files */
#include <cmath>
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include "rt_nonfinite.h"
#include "rtwtypes.h"
#include "MAVEfast_types.h"

/* Function Declarations */
extern void MAVEfast(emxArray_real_T *x, const emxArray_real_T *y,
                     emxArray_char_T *method, double max_dim, double screen,
                     emxArray_real_T *BB, emxArray_real_T *ky, emxArray_real_T
                     *BBvs, emxArray_real_T *idx, emxArray_real_T *C);

#endif

/* End of code generation (MAVEfast.h) */
