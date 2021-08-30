/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * f.c
 *
 * Code generation for function 'f'
 *
 */

/* Include files */
#include "f.h"
#include "rt_nonfinite.h"

/* Function Definitions */
void f(const real_T uj[2], real_T result[2])
{
  if (uj[0] == 0.0) {
    result[0] = 0.0;
    result[1] = 0.0;
  } else {
    result[0] = uj[1];
    result[1] = 4.905 * (uj[0] * uj[0]) + uj[1] * uj[1] / uj[0];
  }
}

/* End of code generation (f.c) */
