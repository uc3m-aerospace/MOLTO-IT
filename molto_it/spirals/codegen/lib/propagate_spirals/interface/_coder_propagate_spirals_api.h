/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: _coder_propagate_spirals_api.h
 *
 * MATLAB Coder version            : 4.1
 * C/C++ source code generated on  : 26-Jun-2019 16:58:51
 */

#ifndef _CODER_PROPAGATE_SPIRALS_API_H
#define _CODER_PROPAGATE_SPIRALS_API_H

/* Include Files */
#include "tmwtypes.h"
#include "mex.h"
#include "emlrt.h"
#include <stddef.h>
#include <stdlib.h>
#include "_coder_propagate_spirals_api.h"

/* Variable Declarations */
extern emlrtCTX emlrtRootTLSGlobal;
extern emlrtContext emlrtContextGlobal;

/* Function Declarations */
extern void propagate_spirals(real_T v0, real_T r0, real_T theta0, real_T psi0,
  real_T thetaf, real_T ee, real_T *Time, real_T *v, real_T *r, real_T *theta,
  real_T *psi, real_T *flag);
extern void propagate_spirals_api(const mxArray * const prhs[6], int32_T nlhs,
  const mxArray *plhs[6]);
extern void propagate_spirals_atexit(void);
extern void propagate_spirals_initialize(void);
extern void propagate_spirals_terminate(void);
extern void propagate_spirals_xil_terminate(void);

#endif

/*
 * File trailer for _coder_propagate_spirals_api.h
 *
 * [EOF]
 */
