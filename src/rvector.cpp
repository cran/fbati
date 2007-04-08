#include "rvector.h"
#include <iostream>
using namespace std;

// setting it up
void RVector::set( double *data, int *dimData )
{
  this->data = data;
  this->I = *dimData;
}

void RVector::set( double *data, int dimData )
{
  this->data = data;
  this->I = dimData;
}
// accessor function
double& RVector::operator () (unsigned i)
{
#ifdef DEBUG_RVECTOR
  if( i<0 || i>=I ) {
    cout << "Index (" << i << ") is out of range " << I << "." << endl;
  }
#endif
  return( data[i] );
}

// repeated code from above
double& RVector::elt( unsigned i )
{
#ifdef DEBUG_RVECTOR
  if( i<0 || i>=I ) {
    cout << "Index (" << i << ") is out of range " << I << "." << endl;
  }
#endif
  return( data[i] );
}
