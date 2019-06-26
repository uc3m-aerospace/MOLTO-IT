/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: _coder_propagate_spirals_api.c
 *
 * MATLAB Coder version            : 4.1
 * C/C++ source code generated on  : 26-Jun-2019 16:58:51
 */

/* Include Files */
#include "tmwtypes.h"
#include "_coder_propagate_spirals_api.h"
#include "_coder_propagate_spirals_mex.h"

/* Variable Definitions */
emlrtCTX emlrtRootTLSGlobal = NULL;
emlrtContext emlrtContextGlobal = { true,/* bFirstTime */
  false,                               /* bInitialized */
  131467U,                             /* fVersionInfo */
  NULL,                                /* fErrorFunction */
  "propagate_spirals",                 /* fFunctionName */
  NULL,                                /* fRTCallStack */
  false,                               /* bDebugMode */
  { 2045744189U, 2170104910U, 2743257031U, 4284093946U },/* fSigWrd */
  NULL                                 /* fSigMem */
};

/* Function Declarations */
static real_T b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId);
static real_T c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId);
static real_T emlrt_marshallIn(const emlrtStack *sp, const mxArray *v0, const
  char_T *identifier);
static const mxArray *emlrt_marshallOut(const real_T u);

/* Function Definitions */

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *u
 *                const emlrtMsgIdentifier *parentId
 * Return Type  : real_T
 */
static real_T b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId)
{
  real_T y;
  y = c_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *src
 *                const emlrtMsgIdentifier *msgId
 * Return Type  : real_T
 */
static real_T c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId)
{
  real_T ret;
  static const int32_T dims = 0;
  emlrtCheckBuiltInR2012b(sp, (const emlrtMsgIdentifier *)msgId, src, "double",
    false, 0U, (int32_T *)&dims);
  ret = *(real_T *)emlrtMxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *v0
 *                const char_T *identifier
 * Return Type  : real_T
 */
static real_T emlrt_marshallIn(const emlrtStack *sp, const mxArray *v0, const
  char_T *identifier)
{
  real_T y;
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = (const char *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = b_emlrt_marshallIn(sp, emlrtAlias(v0), &thisId);
  emlrtDestroyArray(&v0);
  return y;
}

/*
 * Arguments    : const real_T u
 * Return Type  : const mxArray *
 */
static const mxArray *emlrt_marshallOut(const real_T u)
{
  const mxArray *y;
  const mxArray *m0;
  y = NULL;
  m0 = emlrtCreateDoubleScalar(u);
  emlrtAssign(&y, m0);
  return y;
}

/*
 * Arguments    : const mxArray * const prhs[6]
 *                int32_T nlhs
 *                const mxArray *plhs[6]
 * Return Type  : void
 */
void propagate_spirals_api(const mxArray * const prhs[6], int32_T nlhs, const
  mxArray *plhs[6])
{
  real_T v0;
  real_T r0;
  real_T theta0;
  real_T psi0;
  real_T thetaf;
  real_T ee;
  real_T Time;
  real_T v;
  real_T r;
  real_T theta;
  real_T psi;
  real_T flag;
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  st.tls = emlrtRootTLSGlobal;

  /* Marshall function inputs */
  v0 = emlrt_marshallIn(&st, emlrtAliasP(prhs[0]), "v0");
  r0 = emlrt_marshallIn(&st, emlrtAliasP(prhs[1]), "r0");
  theta0 = emlrt_marshallIn(&st, emlrtAliasP(prhs[2]), "theta0");
  psi0 = emlrt_marshallIn(&st, emlrtAliasP(prhs[3]), "psi0");
  thetaf = emlrt_marshallIn(&st, emlrtAliasP(prhs[4]), "thetaf");
  ee = emlrt_marshallIn(&st, emlrtAliasP(prhs[5]), "ee");

  /* Invoke the target function */
  propagate_spirals(v0, r0, theta0, psi0, thetaf, ee, &Time, &v, &r, &theta,
                    &psi, &flag);

  /* Marshall function outputs */
  plhs[0] = emlrt_marshallOut(Time);
  if (nlhs > 1) {
    plhs[1] = emlrt_marshallOut(v);
  }

  if (nlhs > 2) {
    plhs[2] = emlrt_marshallOut(r);
  }

  if (nlhs > 3) {
    plhs[3] = emlrt_marshallOut(theta);
  }

  if (nlhs > 4) {
    plhs[4] = emlrt_marshallOut(psi);
  }

  if (nlhs > 5) {
    plhs[5] = emlrt_marshallOut(flag);
  }
}

/*
 * Arguments    : void
 * Return Type  : void
 */
void propagate_spirals_atexit(void)
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
  propagate_spirals_xil_terminate();
}

/*
 * Arguments    : void
 * Return Type  : void
 */
void propagate_spirals_initialize(void)
{
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  mexFunctionCreateRootTLS();
  st.tls = emlrtRootTLSGlobal;
  emlrtClearAllocCountR2012b(&st, false, 0U, 0);
  emlrtEnterRtStackR2012b(&st);
  emlrtFirstTimeR2012b(emlrtRootTLSGlobal);
}

/*
 * Arguments    : void
 * Return Type  : void
 */
void propagate_spirals_terminate(void)
{
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  st.tls = emlrtRootTLSGlobal;
  emlrtLeaveRtStackR2012b(&st);
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
}

/*
 * File trailer for _coder_propagate_spirals_api.c
 *
 * [EOF]
 */
