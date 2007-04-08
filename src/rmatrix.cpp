#include "rmatrix.h"
#include <iostream>
using namespace std;

// Call to initialize
void RMatrix::set( double *data, unsigned int R, unsigned int C )
{
  this->data = data;
  this->R = R;
  this->C = C;
}
void RMatrix::set( double *data, int *dimData )
{
  set( data, dimData[0], dimData[1] );
}

// accessor function
double& RMatrix::operator() ( unsigned int r, unsigned int c )
{
#ifdef DEBUG_RMATRIX
  if( r<0 || r>=R || c<0 || c>=C ) {
    cout << "Index (" << r << "," << c << ") is out of range (" << R << "," << C << "). Seg-fault, R choose you!" << endl;
  }
#endif
  return( data[ r + c*R ] );
}

double& RMatrix::elt( unsigned r, unsigned c )
{
#ifdef DEBUG_RMATRIX
  if( r<0 || r>=R || c<0 || c>=C ) {
    cout << "Index (" << r << "," << c << ") is out of range (" << R << "," << C << "). Seg-fault, R choose you!" << endl;
  }
#endif
  return( data[ r + c*R ] );
}


/*
// R CMD SHLIB rmatrix.cpp
extern "C" {
  void pmat(double *data, int *R, int *C)
  {
    RMatrix m;
    m.set( data, *R, *C );

    cout << *R << " " << *C << endl;

    for( int r=0; r<m.R; r++ ) {
      for( int c=0; c<m.C; c++ )
        cout << m(r,c) << " ";
      cout << endl;
    }
  }
}
*/
