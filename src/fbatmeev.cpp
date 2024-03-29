#include <R.h>
#include "fbatmeev.h"

const int MAX_CHILD_MEEV = 20;

void REXP_fbatmeev(RMatrix &ddata, RVector &marker,  double trait, double model,  RVector &RET_stat, RVector &RET_numInf)
{
  DataMatrix d;
  d.setGen( ddata.ddata, ddata.R, ddata.C );

  // zero out the returns
  unsigned int m;
  for( m=0; m<marker.I/2; m++ ) {
    RET_stat(m) = 0.0;
    RET_numInf(m) = 0.0;
  }

  double *num = new double[ marker.I / 2 ]; // ALLOC
  double *den = new double[ marker.I / 2 ]; // ALLOC

  int start=-1, end=-1;
  while( d.getNextFamily( start, end ) ) {
    // Stores \sum_j U_ij (since effectively looping over 'i' above)
    // double tempVec[marker.I]; // sunCC hell
    double *tempVec = new double[marker.I/2];

    // go across each marker
    for( m=0; m<marker.I; m+=2 ) {
      int m0 = (int)marker(m);
      int m1 = (int)marker(m+1);

      // first find all the children and parental alleles
      int p1[2]; p1[0] = 0; p1[1] = 0;
      int p2[2]; p2[0] = 0; p2[1] = 0;
      int curP=0;
      int ca[MAX_CHILD_MEEV];
      int cb[MAX_CHILD_MEEV];
      double cEnv[MAX_CHILD_MEEV];
      int curC=0; // n
      double y[MAX_CHILD_MEEV];
      
      // 2020-07-07 -- prevent some compiler target_mem_ref unfounded warnings...
      memset(ca, 0, MAX_CHILD_MEEV*sizeof(int));
      memset(cb, 0, MAX_CHILD_MEEV*sizeof(int));
      memset(cEnv, 0, MAX_CHILD_MEEV*sizeof(double));
      memset(y, 0, MAX_CHILD_MEEV*sizeof(double));
      // 2020-07-07 end

      //int crow[MAX_CHILD]; // DEBUG ONLY

      for( int i=start; i<=end; i++ ) {
        //if( d(i, m0)!=0 && d(i,m1)!=0 ) {
        if( d(i, C_FATH)==0 && d(i, C_MOTH)==0 ) {
          // it's a parent
          if( curP>1 ) {
            Rprintf("Too many parents in family %d, ignoring parents 3 and higher.\n", (int)d(i,C_PID));
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
          y[curC] = d(i,(int)trait);
          //crow[curC] = i; // DEBUG ONLY
          curC++;
        }
      }//i

      // now compute the temp vector
      double ex = fbat_EXS( curC,  p1, p2,  ca, cb,  y,  (int)model );
      tempVec[m/2] = 0.0;
      for( int j=0; j<curC; j++ ) {
        double txmex = y[j] * ( xCode( ca[j], cb[j], (int)model ) - ex );
        if( !ISNAN(txmex) && !ISNAN(cEnv[j]) && ca[j]!=0 && cb[j]!=0 ) {
          tempVec[m/2] += txmex;
        }
      }//j
    }//m

    /*
    cout << "tempVec ";
    for( m=0; m<marker.I; m++ ) cout << tempVec[m] << " ";
    cout << endl;
    */

    // now use the tempVec information to update a and b
    for( m=0; m<marker.I/2; m++ ) {
      num[m] += tempVec[m];
      den[m] += tempVec[m] * tempVec[m];

      if( tempVec[m] > 0.0005 || tempVec[m] < -0.0005 )
        RET_numInf(m)++;
    }//m

    delete [] tempVec;
  }

  // finally compute the test statistics...
  for( m=0; m<marker.I/2; m++ ) {
    RET_stat(m) = num[m]*num[m] / den[m]; // chi-squared distribution, 1 d.f.
    //cout << "num = " << num[m] << " den = " << den[m] << " num*num/den = " << RET_stat(m) << endl;
  }

  // delete
  delete [] num; // DE-ALLOC
  delete [] den; // DE-ALLOC
}

extern "C" {
// AUTO-GENERATED by RExport.R
  void eREXP_fbatmeev(  double *ddata, int *dimdata, double *marker, int *dimmarker, double *trait, double *model, double *RET_stat, int *dimRET_stat, double *RET_numInf, int *dimRET_numInf ) {
    RMatrix _data;
    _data.set(ddata,dimdata);
    RVector _marker;
    _marker.set(marker,dimmarker);
    RVector _RET_stat;
    _RET_stat.set(RET_stat,dimRET_stat);
    RVector _RET_numInf;
    _RET_numInf.set(RET_numInf,dimRET_numInf);

    REXP_fbatmeev(_data, _marker, *trait, *model, _RET_stat, _RET_numInf );
  }
}





void REXP_fbatme(RMatrix &ddata, RVector &marker,  double trait, double model,  RVector &RET_stat, RVector &RET_numInf)
{
  DataMatrix d;
  d.setGen( ddata.ddata, ddata.R, ddata.C );

  // zero out the returns
  unsigned int m;
  for( m=0; m<marker.I/2; m++ ) {
    RET_stat(m) = 0.0;
    RET_numInf(m) = 0.0;
  }

  double *num = new double[ marker.I / 2 ]; // ALLOC
  double *den = new double[ marker.I / 2 ]; // ALLOC

  int start=-1, end=-1;
  while( d.getNextFamily( start, end ) ) {
    // go across each marker
    for( m=0; m<marker.I; m+=2 ) {
      int m0 = (int)marker(m);
      int m1 = (int)marker(m+1);

      // first find all the children and parental alleles
      int p1[2]; p1[0] = 0; p1[1] = 0;
      int p2[2]; p2[0] = 0; p2[1] = 0;
      int curP=0;
      int ca[MAX_CHILD_MEEV];
      int cb[MAX_CHILD_MEEV];
      //////double cEnv[MAX_CHILD_MEEV];
      int curC=0; // n
      double y[MAX_CHILD_MEEV];

      //int crow[MAX_CHILD]; // DEBUG ONLY

      for( int i=start; i<=end; i++ ) {
        //if( d(i, m0)!=0 && d(i,m1)!=0 ) {
        if( d(i, C_FATH)==0 && d(i, C_MOTH)==0 ) {
          // it's a parent
          if( curP>1 ) {
            Rprintf("Too many parents in family %d, ignoring parents 3 and higher.\n", (int)d(i,C_PID));
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
          y[curC] = d(i,(int)trait);
          //crow[curC] = i; // DEBUG ONLY
          curC++;
        }
      }//i

      double Vi = 0.0;
      //cout << curC << endl;
      double Si = fbat_Si( curC, p1, p2, ca, cb, y, (int)model, Vi, 0, curC );
      //cout << "Vi=" << Vi << " Si=" << Si << endl;

      num[m/2] += Si;
      den[m/2] += Vi;

      if( Si > 0.0005 || Si < -0.0005 )
        RET_numInf(m/2)++;
    }//m
  }

  // finally compute the test statistics...
  for( m=0; m<marker.I/2; m++ ) {
    RET_stat(m) = num[m]*num[m] / den[m]; // chi-squared distribution, 1 d.f.
    //cout << "num = " << num[m] << " den = " << den[m] << " num*num/den = " << RET_stat(m) << endl;
  }

  // delete
  delete [] num; // DE-ALLOC
  delete [] den; // DE-ALLOC
}

extern "C" {
// AUTO-GENERATED by RExport.R
  void eREXP_fbatme(  double *ddata, int *dimdata, double *marker, int *dimmarker, double *trait, double *model, double *RET_stat, int *dimRET_stat, double *RET_numInf, int *dimRET_numInf ) {
    RMatrix _data;
    _data.set(ddata,dimdata);
    RVector _marker;
    _marker.set(marker,dimmarker);
    RVector _RET_stat;
    _RET_stat.set(RET_stat,dimRET_stat);
    RVector _RET_numInf;
    _RET_numInf.set(RET_numInf,dimRET_numInf);

    REXP_fbatme(_data, _marker, *trait, *model, _RET_stat, _RET_numInf );
  }
}
