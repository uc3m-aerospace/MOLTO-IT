/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: _coder_propagate_spirals_mex.cpp
 *
 * MATLAB Coder version            : 4.1
 * C/C++ source code generated on  : 26-Jun-2019 16:58:51
 */

/* Include Files */
#include "_coder_propagate_spirals_api.h"
#include "_coder_propagate_spirals_mex.h"

/* Function Declarations */
static void propagate_spirals_mexFunction(int32_T nlhs, mxArray *plhs[6],
  int32_T nrhs, const mxArray *prhs[6]);

/* Function Definitions */

/*
 * Arguments    : int32_T nlhs
 *                mxArray *plhs[6]
 *                int32_T nrhs
 *                const mxArray *prhs[6]
 * Return Type  : void
 */
static void propagate_spirals_mexFunction(int32_T nlhs, mxArray *plhs[6],
  int32_T nrhs, const mxArray *prhs[6])
{
  const mxArray *outputs[6];
  int32_T b_nlhs;
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  st.tls = emlrtRootTLSGlobal;

  /* Check for proper number of arguments. */
  if (nrhs != 6) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:WrongNumberOfInputs", 5, 12, 6, 4,
                        17, "propagate_spirals");
  }

  if (nlhs > 6) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:TooManyOutputArguments", 3, 4, 17,
                        "propagate_spirals");
  }

  /* Call the function. */
  propagate_spirals_api(prhs, nlhs, outputs);

  /* Copy over outputs to the caller. */
  if (nlhs < 1) {
    b_nlhs = 1;
  } else {
    b_nlhs = nlhs;
  }

  emlrtReturnArrays(b_nlhs, plhs, outputs);
}

/*
 * Arguments    : int32_T nlhs
 *                mxArray * const plhs[]
 *                int32_T nrhs
 *                const mxArray * const prhs[]
 * Return Type  : void
 */
void mexFunction(int32_T nlhs, mxArray *plhs[], int32_T nrhs, const mxArray
                 *prhs[])
{
  mexAtExit(propagate_spirals_atexit);

  /* Module initialization. */
  propagate_spirals_initialize();

  /* Dispatch the entry-point. */
  propagate_spirals_mexFunction(nlhs, plhs, nrhs, prhs);

  /* Module termination. */
  propagate_spirals_terminate();
}

/*
 * Arguments    : void
 * Return Type  : emlrtCTX
 */
emlrtCTX mexFunctionCreateRootTLS(void)
{
  emlrtCreateRootTLS(&emlrtRootTLSGlobal, &emlrtContextGlobal, NULL, 1);
  return emlrtRootTLSGlobal;
}

/*
 * File trailer for _coder_propagate_spirals_mex.cpp
 *
 * [EOF]
 */
