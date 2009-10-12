#ifndef _fbatmeev_h_
#define _fbatmeev_h_

#include "rvector.h"
#include "rmatrix.h"
#include "datamatrix.h"
#include "fbatdist.h"

void REXP_fbatmeeev( RMatrix &data, RVector &marker,  double trait, double model,  RVector &RET_stat );

#endif
