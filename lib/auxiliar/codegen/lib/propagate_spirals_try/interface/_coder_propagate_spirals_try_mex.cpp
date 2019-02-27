/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * _coder_propagate_spirals_try_mex.cpp
 *
 * Code generation for function '_coder_propagate_spirals_try_mex'
 *
 */

/* Include files */
#include "_coder_propagate_spirals_try_api.h"
#include "_coder_propagate_spirals_try_mex.h"

/* Function Declarations */
static void c_propagate_spirals_try_mexFunc(int32_T nlhs, mxArray *plhs[6],
  int32_T nrhs, const mxArray *prhs[6]);

/* Function Definitions */
static void c_propagate_spirals_try_mexFunc(int32_T nlhs, mxArray *plhs[6],
  int32_T nrhs, const mxArray *prhs[6])
{
  const mxArray *inputs[6];
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
                        21, "propagate_spirals_try");
  }

  if (nlhs > 6) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:TooManyOutputArguments", 3, 4, 21,
                        "propagate_spirals_try");
  }

  /* Temporary copy for mex inputs. */
  if (0 <= nrhs - 1) {
    memcpy((void *)&inputs[0], (void *)&prhs[0], (uint32_T)(nrhs * (int32_T)
            sizeof(const mxArray *)));
  }

  /* Call the function. */
  propagate_spirals_try_api(inputs, outputs);

  /* Copy over outputs to the caller. */
  if (nlhs < 1) {
    b_nlhs = 1;
  } else {
    b_nlhs = nlhs;
  }

  emlrtReturnArrays(b_nlhs, plhs, outputs);

  /* Module termination. */
  propagate_spirals_try_terminate();
}

void mexFunction(int32_T nlhs, mxArray *plhs[], int32_T nrhs, const mxArray
                 *prhs[])
{
  mexAtExit(propagate_spirals_try_atexit);

  /* Initialize the memory manager. */
  /* Module initialization. */
  propagate_spirals_try_initialize();

  /* Dispatch the entry-point. */
  c_propagate_spirals_try_mexFunc(nlhs, plhs, nrhs, prhs);
}

emlrtCTX mexFunctionCreateRootTLS(void)
{
  emlrtCreateRootTLS(&emlrtRootTLSGlobal, &emlrtContextGlobal, NULL, 1);
  return emlrtRootTLSGlobal;
}

/* End of code generation (_coder_propagate_spirals_try_mex.cpp) */
