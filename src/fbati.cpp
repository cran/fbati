#include "fbati.h"

#include <R.h>

#include <iostream>
using namespace std;

/*
// not numerically stable!
double sumProd( double *a, double *b, int N ) {
  double sum=0.0;
  for( int n=0; n<N; n++ )
    sum += a[n]*b[n];
  return(sum);
}*/


// WARNING: v will be overwritten
double sum( double* v, int LEN )
{
  // some special cases
  if( LEN == 0 ) {
    return(R_NaN); //( NAN );
  }else if( LEN == 1 ){
    return( v[0] );
  }else if( LEN == 2 ){
    return( v[0] + v[1] );
  }

  // now on to the general algorithm
  int stepsize = 2;
  while( (stepsize/2) < LEN ){
    for( int i=0; i<LEN; i+=stepsize ){
      if( i+(stepsize/2) < LEN )
        v[i] += v[i+(stepsize/2)];
    }
    stepsize *= 2;
  }
  return( v[0] );
}

double sumProd( double *a, double *b, int N ) {
  //double c[N];
  //for( int n=0; n<N; n++ )
  //  c[n] = a[n] * b[n];
  //return( sum( c, N ) );

  double *c = new double[N];
  for( int n=0; n<N; n++ )
    c[n] = a[n] * b[n];
  double ssum = sum( c, N );
  delete [] c;
  return( ssum );
}

extern "C" {
  void fbati( double *pRes, // the result

              int *pN,   // the length of the upcoming vectors
              double *xmxbar,  // all of these vectors are sorted by group
              double *zmzbar,
              int *group,

              int *iter ) {
    int N = *pN;

    rndAttach();  // NEW NEW NEW!!!

    // get the value as it's currently
    double stat = sumProd( xmxbar, zmzbar, N );

    //cout << "STAT: " << stat << endl;

    // set up the permutation cutoffs
    int ngroups=0;
    int groupStart[100];
    int groupEnd[100];
    for( int n=0; n<N; n++ ) {
      if( n==0 ) {
        ngroups=1;
        groupStart[ngroups-1] = n;
      }else if( group[n-1] != group[n] ) {
        groupEnd[ngroups-1] = n-1;
        groupStart[ngroups] = n;
        ngroups++;
      }
    }
    groupEnd[ngroups-1] = N-1;

    int moreExtreme = 0; // actually AS or more extreme
    // now get it for all the permutations
    int niter = *iter;
    for( int i=0; i<niter; i++ ) {
      for( int g=0; g<ngroups; g++ ) {
        int LEN = groupEnd[g]-groupStart[g]+1;

        // permute z within that group
        for( int k=0; k<LEN; k++ ) {
          int randN = RandInt( k, LEN-1 );
          double temp = zmzbar[ groupStart[g] + k ];
          zmzbar[ groupStart[g] + k ] = zmzbar[ groupStart[g] + randN ];
          zmzbar[ groupStart[g] + randN ] = temp;
        }
      }

      // and get the new sum
      double newstat = sumProd( xmxbar, zmzbar, N );
      //cout << newstat << " ";
      //cout << "";

      if( fabs(newstat) >= fabs(stat) )
        moreExtreme++;
      //cout << (newstat-stat) << " ";
    }
    //cout << endl;
    //cout << "moreExtreme " << moreExtreme << endl;

    // lastly, compute the pvalue
    pRes[0] = ((double)moreExtreme)/((double)niter);

    // NEW NEW NEW!!!
    rndDetach();
  }
}
