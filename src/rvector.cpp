#include <R.h>

#include "rvector.h"
#include <iostream>
using namespace std;

// setting it up
void RVector::set( double *p_data, int *p_dimData )
{
  this->data = p_data;
  this->I = *p_dimData;
}

void RVector::set( double *p_data, int p_dimData )
{
  this->data = p_data;
  this->I = p_dimData;
}
// accessor function
double& RVector::operator () (unsigned i)
{
#ifdef DEBUG_RVECTOR
  //if( i<0 || i>=I ) { // can't be <0, by unsigned int casting
  if( i>=I ) {
    Rprintf("Index (%d) is out of range (%d).\n", i, I);
  }
#endif
  return( data[i] );
}

// repeated code from above
double& RVector::elt( unsigned i )
{
#ifdef DEBUG_RVECTOR
  //if( i<0 || i>=I ) {
  if( i>=I ) {
    Rprintf("Index (%d) is out of range (%d).\n", i, I);
  }
#endif
  return( data[i] );
}
