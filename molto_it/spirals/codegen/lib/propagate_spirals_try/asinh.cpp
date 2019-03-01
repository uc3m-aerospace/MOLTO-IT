/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * asinh.cpp
 *
 * Code generation for function 'asinh'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "propagate_spirals_try.h"
#include "asinh.h"

/* Function Definitions */
void b_asinh(double *x)
{
  boolean_T xneg;
  double t;
  double absz;
  xneg = (*x < 0.0);
  if (xneg) {
    *x = -*x;
  }

  if (*x >= 2.68435456E+8) {
    *x = std::log(*x) + 0.69314718055994529;
  } else if (*x > 2.0) {
    *x = std::log(2.0 * *x + 1.0 / (std::sqrt(*x * *x + 1.0) + *x));
  } else {
    t = *x * *x;
    t = *x + t / (1.0 + std::sqrt(1.0 + t));
    *x = t;
    absz = std::abs(t);
    if ((absz > 4.503599627370496E+15) || (!((!rtIsInf(t)) && (!rtIsNaN(t))))) {
      *x = std::log(1.0 + t);
    } else {
      if (!(absz < 2.2204460492503131E-16)) {
        *x = std::log(1.0 + t) * (t / ((1.0 + t) - 1.0));
      }
    }
  }

  if (xneg) {
    *x = -*x;
  }
}

/* End of code generation (asinh.cpp) */
