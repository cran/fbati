#ifndef _rmatrix_h_
#define _rmatrix_h_

#include "defines.h"

class RMatrix
{
  public:
    //setting it up
    void set( double *p_data, unsigned int p_R, unsigned int p_C );
    void set( double *p_data, int *p_dimData );
    //accessor function
    double& operator() ( unsigned r, unsigned c );
    double& elt( unsigned r, unsigned c );
  public:
    double *data;
    unsigned int R, C;
};

#endif
