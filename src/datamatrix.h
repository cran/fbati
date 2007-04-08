/* Provides a few extra constants and vars
 *  to the RMatrix class to deal with the data
 *  in the gene by environment case.
 *
 * The format is essentially a merged pedigree/phenotype file
 */

#ifndef _datamatrix_h_
#define _datamatrix_h_

#include <cstring>
#include <string>
#include <stdlib.h>
using namespace std;

#include "rmatrix.h"
#include "fbatdist.h"


// Constants based on the data
const int C_PID  = 0;
const int C_ID   = 1;
const int C_FATH = 2;
const int C_MOTH = 3;
const int C_SEX  = 4;
const int C_AFF  = 5;

class DataMatrix : public RMatrix
{
  public:
    // 'constants' that need to be initialized
    // - these need to be set for computeGroupG
    int c_m0;
    int c_m1;
    // - this is for use elsewhere...
    int c_env;

  public:
    // fills in the constants above as well,
    //  this is for the generation routines...
    void setGen( double *data, unsigned int R, unsigned int C=9 );

    bool getNextFamily( int &start, int &end );
    void computeGroupG( int *groups,
                        double *g0, double *g1, double *g2,
                        int *affected_index,
                        int &affected_index_size,
                        int &data_num_families );

    void genPush( int pid, int id, int idfath, int idmoth,
                  int sex, int affection,
                  int m0, int m1,
                  double env,
                  int curRow
                );
};

extern "C" {
  void dataComputeGroupG( double *data, int *dataDim,
                          int *m0pos, int *m1pos,
                          int *groups,
                          double *g0, double *g1, double *g2,
                          int *affected_index,
                          int *affected_index_size,
                          int *data_num_families );
}

#endif
