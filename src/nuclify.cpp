// R CMD SHLIB nuclify.cpp datamatrix.cpp rmatrix.cpp fbatdist.cpp random.cpp
#include <R.h>

#include <limits>  // for quiet_NaN

#include "nuclify.h"
#define NUCLIFY_MAX_SIBS 256
#define SEX_MALE 1
#define SEX_FEMALE 2


// needed for strataReduce
// -- could be written faster, written for safety!
void strataReduceRemove( int *array, int &arraySize, int elt ) {
  // find element, and replace it with the last element
  for( int i=0; i<arraySize; i++ ) {
    if( array[i]==elt ) {
      array[i] = array[arraySize-1];
      arraySize--;
      return;
    }
  }

  // it wasn't found! uh-oh...
  Rprintf("strataReduceRemove ERROR -- elt %d was not found in the array, and so could not be removed!\n", elt);
}

// clearParents just removes the fathid/mothid references (for nuclify)
// fullyClearParents removes the fathid/mothid alleles as well
int pushDataRow( DataMatrix &din, int dinI, DataMatrix &dout, int doutI, int newPid, bool clearParents=false, bool fullyClearParents=false, bool clearAffection=false, int envCol=-1 ) {
  // copy it in
  for( int c=0; c<(int)din.C; c++ ) {
    dout(doutI,c) = din(dinI,c);
    //cout << " " << c << "  " << din(dinI,c) << endl;
  }
  // and reset the pid
  dout(doutI,C_PID) = newPid;

  // should the parents of this one be erased?
  if( clearParents ) {
    dout(doutI,C_FATH) = 0;
    dout(doutI,C_MOTH) = 0;
  }
  // fully clear the parents? (alleles as well?)
  if( fullyClearParents ) {
    for( int c=C_AFF+1; c<(int)din.C; c++ )
      dout(doutI,c) = 0;
  }
  // and clear the affection (for the other sib)
  if( clearAffection ) {
    dout(doutI,C_AFF) = 0; // missing
    // make the env missing is really the only safe way to do this,
    //  as we might be using different sorts of affection statuses
    //  -- this way it should work with them all
    if( envCol!=-1 )
      dout(doutI,envCol) = std::numeric_limits<double>::quiet_NaN();
  }

  return( doutI + 1 );
}
int pushEmptyRow( DataMatrix &dout, int doutI, int newPid, int id, int sex ) {
  // zero it out
  for( int c=0; c<(int)dout.C; c++ )
    dout(doutI,c) = 0;
  // then setup the other stuff
  dout(doutI,C_PID) = newPid;
  dout(doutI,C_ID) = id;
  //dout(doutI,C_FATH) = 0; // already zeroed out...
  //dout(doutI,C_MOTH) = 0;
  dout(doutI,C_SEX) = sex;
  //dout(doutI,C_AFFECTION) = 0;

  return( doutI + 1 );
}

// prerequisite: dimData[1] = dimDataOut[1]
// up to 100 sibs (or newPid fails)
extern "C" {
  // Chops pedigrees up into nuclear families
  // failure=0 == success
  // failure=1 == indicates the output isn't large enough (suggest doubling it, which is what the R code does in a recursive fashion)
  void nuclify_cpp( double *data, int *dimData,
                double *dataOut, int *dimDataOut,
                int *failure ) {
    *failure = 0;

    // set up the DataMatrix objects
    DataMatrix din, dout;
    din.set( data, dimData );
    dout.set( dataOut, dimDataOut );

    // do the nuclification
    int start=-1, end=-1;
    //int cpid = 1; // apparently unused?
    int curRow = 0;
    // start by going across all families
    while( din.getNextFamily(start,end) ) {
      int newPid = (int)din(start,C_PID) * 100;
      // pull out each unique mother father pair
      for( int i=start; i<=end; i++ ) {
        int idfath = (int)din(i,C_FATH);
        int idmoth = (int)din(i,C_MOTH);

        int idfathRow = -1;
        int idmothRow = -1;

        int numSibs = 0;
        int sibs[NUCLIFY_MAX_SIBS];
        bool validFamily = true;
        for( int ii=start; ii<=end && validFamily; ii++ ) {
          // We actually _want_ when i==ii to fall through...

          if( idfath==din(ii,C_FATH) && idmoth==din(ii,C_MOTH) ) {
            // we know it's a new/unique sub-family if it hasn't been found before...
            if( ii < i ) {
              validFamily = false;
            }else{
              // we've got a new family!
              sibs[numSibs] = ii;
              numSibs++;
            }
          }else if( din(ii,C_ID) == idfath ) {
            idfathRow = ii;
          }else if( din(ii,C_ID) == idmoth ) {
            idmothRow = ii;
          }
        }
        if( idfath==0 && idmoth==0 ) validFamily = false;
        if( validFamily ) {
          // Then need to push onto the output

          // the father, inserting one if not in the pedigree
          if( idfathRow != -1 ) {
            curRow = pushDataRow( din, idfathRow, dout, curRow, newPid, true );
            if(curRow == (int)dout.R) {*failure = 1; return;}
          }else{
            curRow = pushEmptyRow( dout, curRow, newPid, idfath, SEX_MALE );
            if(curRow == (int)dout.R) {*failure = 1; return;}
          }

          // the mother, ""
          if( idmothRow != -1 ) {
            curRow = pushDataRow( din, idmothRow, dout, curRow, newPid, true );
            if(curRow == (int)dout.R) {*failure = 1; return;}
          }else{
            curRow = pushEmptyRow( dout, curRow, newPid, idmoth, SEX_FEMALE );
            if(curRow == (int)dout.R) {*failure = 1; return;}
          }

          // and all of their children
          for( int ch=0; ch<numSibs; ch++ ) {
            curRow = pushDataRow( din, sibs[ch], dout, curRow, newPid );
            if(curRow == (int)dout.R) {*failure = 1; return;}
          }

          newPid++;
        }
      }
    }

    // set the dims on the output
    dimDataOut[0] = curRow;
  }

  // returns a more strata friendly family configuration
  // i.e. a random affected (with env if possible),
  //  their parents ONLY if both present,
  //  otherwise it returns a random sibpair
  //
  // Assumes data has been through the 'nuclify' routine,
  //  so data is in nuclear families and the children have
  //  parents, but their parents' parents are marked
  //  as missing
  void strataReduce_cpp( double *data, int *dimData,
                     double *dataOut, int *dimDataOut,
                     int *pEnvCol, int *pm0, int *pm1,
                     int *pMaxSib ) { // new 04/08/2008

    rndAttach(); // RANDOM NUMBER GENERATOR SETUP

    int envCol = *pEnvCol;
    int m0 = *pm0, m1 = *pm1;

    int maxSib = *pMaxSib;

    //cout << "c++ dimDataOut " << dimDataOut[0] << " " << dimDataOut[1] << endl;

    // set up the DataMatrix objects
    DataMatrix din, dout;
    din.set( data, dimData );
    dout.set( dataOut, dimDataOut );

    // do the processing
    int start=-1, end=-1;
    //int cpid = 1; // and still unused!!
    int curRow = 0;
    // go across all the families; note some will not be used
    while( din.getNextFamily(start,end) ) {
      int pid = (int)din(start,C_PID);

      // find the mother/father id's
      int idmoth=-1;
      int idfath=-1;
      int i;
      for( i=start; i<=end; i++ ) {
        if( din(i,C_MOTH)!=0 )
          idmoth = (int)din(i,C_MOTH);
        if( din(i,C_FATH)!=0 )
          idfath = (int)din(i,C_FATH);
      }

      // find all children and parents locations for current family
      // - basic children/parents loc
      //int childRow[100];   // not actually
      int childRowSize=0;  //  needed?
      int parentRow[20];
      int parentRowSize=0;
      // - affected, environmental info, genotyped
      int childRowAffected[100];
      int childRowAffectedSize=0;
      // - genotyped (for sibpair)
      int childRowGeno[100];
      int childRowGenoSize=0;
      // -- so the childRowAffected set is contained in childRowGeno set
      for( i=start; i<=end; i++ ) {
        if( din(i,C_ID)==idfath || din(i,C_ID)==idmoth ) {
          if( parentRowSize==2 ) {
          }else{
            parentRow[parentRowSize] = i;
            parentRowSize++;
          }
        }else{
          //childRow[childRowSize] = i;
          childRowSize++;

          if( din(i,C_AFF)==2                  // affected
              && !ISNAN(din(i,envCol))         // environment
              && din(i,m0)!=0 && din(i,m1)!=0 // genotyped
          ){
            //cout << "affected, env present, genotyped: " << i << endl;
            childRowAffected[childRowAffectedSize] = i;
            childRowAffectedSize++;
          }

          if( din(i,m0)!=0 && din(i,m1)!=0 ) { // genotyped
            childRowGeno[childRowGenoSize] = i;
            childRowGenoSize++;
          }
        }
      }

      // Now push on the friendlier strata data
      // - are the parents there, and allele's not missing?
      if( parentRowSize==2
          && din(parentRow[0],m0)!=0 && din(parentRow[0],m1)!=0
          && din(parentRow[1],m0)!=0 && din(parentRow[1],m1)!=0 ) {
        // yes, we've got both parents!

        // have at least a usable affected?
        if( childRowAffectedSize>0 ) {
          // push the parents on
          for( int p=0; p<2; p++ )
            curRow = pushDataRow( din, parentRow[p], dout, curRow, pid, true );

          // choose a random affected
          curRow = pushDataRow( din, childRowAffected[RandInt(0,childRowAffectedSize-1)], dout, curRow, pid );
        }
      }else{
        // no, we don't

        // have at least a usable affected, and another genotyped?
        if( childRowAffectedSize>0 && childRowGenoSize>1 ) {

          /* GOOD CODE THAT WORKS FOR SIBPAIRS
          // get a random affected (with env, genoed)
          int raff = RandInt(0,childRowAffectedSize-1);
          // get a random other genotyped individual
          int rother = RandInt(0,childRowGenoSize-1);
          //while( rother==raff ) rother = RandInt(0,childRowGenoSize-1);
          while( childRowAffected[raff] == childRowGeno[rother] ) rother = RandInt(0,childRowGenoSize-1);
          // ABOVE WAS A NASTY LITTLE BUG??

          //cout << "raff " << raff << ", rother " << rother << ", childRowAffectedSize " << childRowAffectedSize << ", childRowGenoSize " << childRowGenoSize << endl;
          //cout << childRowAffected[raff] << " -- " << childRowGeno[rother] << endl;

          // push their parents on (as missing)
          for( int p=0; p<2; p++ )
            curRow = pushDataRow( din, parentRow[p], dout, curRow, pid, true, true );

          // and push them both on
          curRow = pushDataRow( din, childRowAffected[raff], dout, curRow, pid );
          curRow = pushDataRow( din, childRowGeno[rother],   dout, curRow, pid, false, false, true );
          */

          // New code that works for more generalized strata
          int affectedRow = childRowAffected[ RandInt( 0, childRowAffectedSize-1 ) ];
          curRow = pushDataRow( din, affectedRow, dout, curRow, pid );
          strataReduceRemove( childRowGeno, childRowGenoSize, affectedRow );

          if( childRowGenoSize+1 <= maxSib ) { // +1 for affected
            // push all offspring on
            for( int c=0; c<childRowGenoSize; c++ )
              curRow = pushDataRow( din, childRowGeno[c], dout, curRow, pid, false, false, true, envCol );
          }else{
            // random # of offspring to push on
            for( int c=1; c<maxSib && childRowGenoSize>0; c++ ) { // +1 for affected
              int unaffectedRow = childRowGeno[ RandInt( 0, childRowGenoSize-1 ) ]; // not really unaffected, but will be marked that way
              curRow = pushDataRow( din, unaffectedRow, dout, curRow, pid, false, false, true );
              strataReduceRemove( childRowGeno, childRowGenoSize, unaffectedRow );
            }
          }
        }
      }
    }

    // set the dims on the output
    dimDataOut[0] = curRow;
    //cout << "dimDataOut[0] " << dimDataOut[0] << " curRow " << curRow << endl;

    rndDetach(); // RANDOM NUMBER GENERATOR CLEANUP
  }
}
