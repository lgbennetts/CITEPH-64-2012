/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * _coder_KTSpace_HLLE_mex.c
 *
 * Code generation for function '_coder_KTSpace_HLLE_mex'
 *
 */

/* Include files */
#include "_coder_KTSpace_HLLE_mex.h"
#include "KTSpace_HLLE_data.h"
#include "KTSpace_HLLE_initialize.h"
#include "KTSpace_HLLE_terminate.h"
#include "_coder_KTSpace_HLLE_api.h"
#include "rt_nonfinite.h"

/* Function Definitions */
void KTSpace_HLLE_mexFunction(int32_T nlhs, mxArray *plhs[2], int32_T nrhs,
  const mxArray *prhs[5])
{
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  const mxArray *outputs[2];
  int32_T b_nlhs;
  st.tls = emlrtRootTLSGlobal;

  /* Check for proper number of arguments. */
  if (nrhs != 5) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:WrongNumberOfInputs", 5, 12, 5, 4,
                        12, "KTSpace_HLLE");
  }

  if (nlhs > 2) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:TooManyOutputArguments", 3, 4, 12,
                        "KTSpace_HLLE");
  }

  /* Call the function. */
  KTSpace_HLLE_api(prhs, nlhs, outputs);

  /* Copy over outputs to the caller. */
  if (nlhs < 1) {
    b_nlhs = 1;
  } else {
    b_nlhs = nlhs;
  }

  emlrtReturnArrays(b_nlhs, plhs, outputs);
}

void mexFunction(int32_T nlhs, mxArray *plhs[], int32_T nrhs, const mxArray
                 *prhs[])
{
  mexAtExit(&KTSpace_HLLE_atexit);

  /* Module initialization. */
  KTSpace_HLLE_initialize();

  /* Dispatch the entry-point. */
  KTSpace_HLLE_mexFunction(nlhs, plhs, nrhs, prhs);

  /* Module termination. */
  KTSpace_HLLE_terminate();
}

emlrtCTX mexFunctionCreateRootTLS(void)
{
  emlrtCreateRootTLS(&emlrtRootTLSGlobal, &emlrtContextGlobal, NULL, 1);
  return emlrtRootTLSGlobal;
}

/* End of code generation (_coder_KTSpace_HLLE_mex.c) */
