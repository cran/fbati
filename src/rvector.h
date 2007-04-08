#ifndef _rvector_h_
#define _rvector_h_

#include "defines.h"

class RVector
{
  public:
    // setting it up
    void set( double *data, int *dimData );
    void set( double *data, int dimData );
    // accessor function
    double& operator () ( unsigned i );
    double& elt( unsigned i );
  public:
    double *data;
    unsigned I;
};

#endif
