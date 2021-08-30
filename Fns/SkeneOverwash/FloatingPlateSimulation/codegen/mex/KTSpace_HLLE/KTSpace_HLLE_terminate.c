/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * KTSpace_HLLE_terminate.c
 *
 * Code generation for function 'KTSpace_HLLE_terminate'
 *
 */

/* Include files */
#include "KTSpace_HLLE_terminate.h"
#include "KTSpace_HLLE_data.h"
#include "_coder_KTSpace_HLLE_mex.h"
#include "rt_nonfinite.h"

/* Function Definitions */
void KTSpace_HLLE_atexit(void)
{
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  mexFunctionCreateRootTLS();
  st.tls = emlrtRootTLSGlobal;
  emlrtEnterRtStackR2012b(&st);
  emlrtLeaveRtStackR2012b(&st);
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
  emlrtExitTimeCleanup(&emlrtContextGlobal);
}

void KTSpace_HLLE_terminate(void)
{
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  st.tls = emlrtRootTLSGlobal;
  emlrtLeaveRtStackR2012b(&st);
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
}

/* End of code generation (KTSpace_HLLE_terminate.c) */
