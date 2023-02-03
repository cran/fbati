#ifndef _rvector_h_
#define _rvector_h_

#include "defines.h"

class RVector
{
  public:
    // setting it up
    void set( double *p_data, int *p_dimData );
    void set( double *p_data, int p_dimData );
    // accessor function
    double& operator () ( unsigned i );
    double& elt( unsigned i );
  public:
    double *ddata;
    unsigned I;
};

#endif
