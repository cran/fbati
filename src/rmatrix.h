#ifndef _rmatrix_h_
#define _rmatrix_h_

#include "defines.h"

class RMatrix
{
  public:
    //setting it up
    void set( double *data, unsigned int R, unsigned int C );
    void set( double *data, int *dimData );
    //accessor function
    double& operator() ( unsigned r, unsigned c );
    double& elt( unsigned r, unsigned c );
  public:
    double *data;
    unsigned int R, C;
};

#endif
