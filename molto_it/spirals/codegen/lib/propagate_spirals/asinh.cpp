//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: asinh.cpp
//
// MATLAB Coder version            : 4.1
// C/C++ source code generated on  : 26-Jun-2019 16:58:51
//

// Include Files
#include <cmath>
#include "rt_nonfinite.h"
#include "propagate_spirals.h"
#include "asinh.h"

// Function Definitions

//
// Arguments    : double *x
// Return Type  : void
//
void b_asinh(double *x)
{
  boolean_T xneg;
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
    absz = *x * *x;
    *x += absz / (1.0 + std::sqrt(1.0 + absz));
    absz = std::abs(*x);
    if ((absz > 4.503599627370496E+15) || (rtIsInf(*x) || rtIsNaN(*x))) {
      (*x)++;
      *x = std::log(*x);
    } else {
      if (!(absz < 2.2204460492503131E-16)) {
        *x = std::log(1.0 + *x) * (*x / ((1.0 + *x) - 1.0));
      }
    }
  }

  if (xneg) {
    *x = -*x;
  }
}

//
// File trailer for asinh.cpp
//
// [EOF]
//
