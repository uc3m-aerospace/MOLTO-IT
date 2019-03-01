/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * acosh.cpp
 *
 * Code generation for function 'acosh'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "propagate_spirals_try.h"
#include "acosh.h"

/* Function Definitions */
void b_acosh(double *x)
{
  double z;
  double absz;
  if (*x < 1.0) {
    *x = rtNaN;
  } else if (*x >= 2.68435456E+8) {
    *x = std::log(*x) + 0.69314718055994529;
  } else if (*x > 2.0) {
    *x = std::log(*x + std::sqrt(*x * *x - 1.0));
  } else {
    (*x)--;
    z = *x + std::sqrt(2.0 * *x + *x * *x);
    *x = z;
    absz = std::abs(z);
    if ((absz > 4.503599627370496E+15) || (!((!rtIsInf(z)) && (!rtIsNaN(z))))) {
      *x = std::log(1.0 + z);
    } else {
      if (!(absz < 2.2204460492503131E-16)) {
        *x = std::log(1.0 + z) * (z / ((1.0 + z) - 1.0));
      }
    }
  }
}

/* End of code generation (acosh.cpp) */
