/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * propagate_spirals_try.h
 *
 * Code generation for function 'propagate_spirals_try'
 *
 */

#ifndef PROPAGATE_SPIRALS_TRY_H
#define PROPAGATE_SPIRALS_TRY_H

/* Include files */
#include <cmath>
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include "rt_defines.h"
#include "rt_nonfinite.h"
#include "rtwtypes.h"
#include "propagate_spirals_try_types.h"

/* Function Declarations */
extern void propagate_spirals_try(double v0, double r0, double theta0, double
  psi0, double thetaf, double ee, double *Time, double *v, double *r, double
  *theta, double *psi, double *flag);

#endif

/* End of code generation (propagate_spirals_try.h) */
