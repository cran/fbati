#ifndef _joint_h_
#define _joint_h_

#include "rvector.h"
#include "rmatrix.h"
#include "datamatrix.h"
#include "fbatdist.h"

void REXP_joint( RMatrix &ddata,  RVector &marker, double trait, double env,  double model,  RVector &RET_a, RMatrix &RET_b,  RVector &RET_numInf );

#endif

