/* Thomas Hoffmann
 * nuclify.h
 *
 * Replaces the nuclifyMerged(...) R function
 *  with some more efficient c++ code
 */

#ifndef _nuclify_h_
#define _nuclify_h_

#include "datamatrix.h"

extern "C" {
  void nuclify_cpp( double *data, int *dimData,
                double *dataOut, int *dimDataOut,
                int *failure );
  void strataReduce_cpp( double *data, int *dimData,
                     double *dataOut, int *dimDataOut,
                     int *pEnvCol, int *pm0, int *pm1,
                     int *pMaxSib ); // new 04/08/2008
}
#endif
