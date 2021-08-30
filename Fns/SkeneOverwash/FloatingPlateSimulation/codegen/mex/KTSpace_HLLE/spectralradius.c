/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * spectralradius.c
 *
 * Code generation for function 'spectralradius'
 *
 */

/* Include files */
#include "spectralradius.h"
#include "rt_nonfinite.h"
#include "mwmathutil.h"

/* Function Definitions */
void spectralradius(real_T u[2], real_T lambda[2])
{
  real_T d;
  real_T lambda_tmp;
  if (u[0] < 0.0) {
    u[0] = 0.0;
  }

  d = u[1] / u[0];
  if (muDoubleScalarIsNaN(d) || muDoubleScalarIsInf(d)) {
    lambda[0] = 0.0;
    lambda[1] = 0.0;
  } else {
    lambda_tmp = muDoubleScalarSqrt(9.81 * u[0]);
    lambda[0] = d + lambda_tmp;
    lambda[1] = d - lambda_tmp;
  }
}

/* End of code generation (spectralradius.c) */
