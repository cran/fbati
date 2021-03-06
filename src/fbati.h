#ifndef _fbati_h_
#define _fbati_h_

#include "defines.h"
#include "random.h"
//#include <math.h>

#include <R.h>
#include <Rmath.h>

double sumProd( double *a, double *b, int N );

extern "C" {
  void fbati_cpp( double *pRes, // the result

              int *pN,   // the length of the upcoming vectors
              double *xmxbar,  // all of these vectors are sorted by group
              double *zmzbar,
              int *group,

              int *iter );
}

#endif
