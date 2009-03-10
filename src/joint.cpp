#include "joint.h"
#include "fbatdist.h"

// hmm... maybe...
//#define isnan(x) ((x) != (x))

const int MAX_CHILD = 20;

// then with R we need to compute a^T %*% b^- %*% a ~ chi-sq(rank(b))
void REXP_joint( RMatrix &data,  RVector &marker, double trait, double env,  double model,  RVector &RET_a, RMatrix &RET_b )
{
  DataMatrix d;
  d.setGen( data.data, data.R, data.C );

  /*
  // print out the data -- debug
  for( int rr=0; rr<d.R; rr++ ) {
    for( int cc=0; cc<d.C; cc++ ) {
      cout << d(rr,cc) << " ";
    }
    cout << endl;
  }*/

  // zero out the a and b RMatrix objects
  unsigned int r=0, c=0, m=0;
  for( r=0; r<RET_a.I; r++ )
    RET_a(r) = 0.0;
  for( r=0; r<RET_b.R; r++ )
    for( c=0; c<RET_b.C; c++ )
      RET_b(r,c) = 0.0;

  /*
  cout << "dim(data) " << d.R << " " << d.C << endl;
  cout << "zeroed out" << endl;
  cout << "trait " << (int)trait << endl;
  cout << "env " << (int)env << endl;
  cout << "model " << (int)model << endl;
  */

  int start=-1, end=-1;
  while( d.getNextFamily( start, end ) ) {
    // Stores \sum_j U_ij (since effectively looping over 'i' above)
    // double tempVec[marker.I]; // sunCC hell
    double *tempVec = new double[marker.I];

    // go across each marker
    for( m=0; m<marker.I; m+=2 ) {
      int m0 = (int)marker(m);
      int m1 = (int)marker(m+1);

      // first find all the children and parental alleles
      int p1[2]; p1[0] = 0; p1[1] = 0;
      int p2[2]; p2[0] = 0; p2[1] = 0;
      int curP=0;
      int ca[MAX_CHILD];
      int cb[MAX_CHILD];
      double cEnv[MAX_CHILD];
      int curC=0; // n
      double y[MAX_CHILD];

      //int crow[MAX_CHILD]; // DEBUG ONLY

      for( int i=start; i<=end; i++ ) {
        //if( d(i, m0)!=0 && d(i,m1)!=0 ) {
        if( d(i, C_FATH)==0 && d(i, C_MOTH)==0 ) {
          // it's a parent
          if( curP>1 ) {
            cout << "Too many parents in family " << d(i,C_PID) << ", ignoring parents 3 and higher." << endl;
          }else{
            if( curP==0 ) {
              p1[0] = (int)d(i,m0);
              p1[1] = (int)d(i,m1);
            }else if( curP==1 ) {
              p2[0] = (int)d(i,m0);
              p2[1] = (int)d(i,m1);
            }
            curP++;
          }
        }else{
          // assume a child then
          ca[curC] = (int)d(i,m0);
          cb[curC] = (int)d(i,m1);
          cEnv[curC] = d(i,(int)env);
          y[curC] = d(i,(int)trait);
          //crow[curC] = i; // DEBUG ONLY
          curC++;
        }
      }//i

      // now compute the temp vector
      double ex = fbat_EXS( curC,  p1, p2,  ca, cb,  y,  (int)model );
      tempVec[m]   = 0.0;
      tempVec[m+1] = 0.0;
      for( int j=0; j<curC; j++ ) {
        double txmex = y[j] * ( xCode( ca[j], cb[j], (int)model ) - ex ); // POTENTIALLY WANT TO DO (y[j]-OFFSET)
        //cout << crow[j] << " " << (int)d(crow[j],C_PID) << " " << (int)d(crow[j],C_ID) << ": " << y[j] << "; " << xCode( ca[j], cb[j], (int)model ) - ex << "  " << txmex << endl;
        //cout << "ex=" << ex << " trait=" << y[j] << " txmex=" << txmex << " env=" << cEnv[j] << endl;
        if( !isnan(txmex) && !isnan(cEnv[j]) && ca[j]!=0 && cb[j]!=0 ) {
          tempVec[m] += txmex;
          tempVec[m+1] += txmex * cEnv[j];
        }
      }//j
    }//m

    /*
    cout << "tempVec ";
    for( m=0; m<marker.I; m++ ) cout << tempVec[m] << " ";
    cout << endl;
    */

    // now use the tempVec information to update a and b
    for( m=0; m<marker.I; m++ ) {
      RET_a(m) += tempVec[m];
      for( unsigned int mp=0; mp<marker.I; mp++ )
        RET_b(m,mp) += tempVec[m] * tempVec[mp];
    }//m

    delete [] tempVec;
  }

  /*
  cout << "RET_a= "
  for( m=0; m<marker.I; m++ )
    cout << RET_a(m) << " ";
  cout << endl;
  */
}
