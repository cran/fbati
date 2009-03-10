#include "datamatrix.h"

void DataMatrix::setGen( double *p_data, unsigned int p_R, unsigned int p_C )
{
  set( p_data, p_R, p_C );
  c_m0 = 6;
  c_m1 = 7;
  c_env = 8;
}


bool DataMatrix::getNextFamily( int &start, int &end )
{
  // set up the new start position
  if( start == -1 )
    start = 0;
  else
    start = end + 1;

  // make sure we haven't gotten past the bounds
  if( start >= (int)R )
    return( false );

  // now fill in the end
  int pid = (int)elt(start,C_PID);
  for( int i=start; i<(int)R; i++ ) {
    if( elt(i,C_PID) == pid )
      end = i;
    else
      return( true ); // already reached the end
  }
  return( true ); // for the last pedigree
}

// Only uses the first affected in the family to return in the affected_index
void DataMatrix::computeGroupG( int *groups,
                                double *g0, double *g1, double *g2,
                                int *affected_index,
                                int &affected_index_size,
                                int &data_num_families )
{
  int start=-1, end=-1;

  int prevInformativePid = -1; // for when firstAffectedOnly=true, which has been altered to _always_ be the case, since this is the analysis portion!
  data_num_families = 0;

  // new precaution
  memset( g0, 255, R*sizeof(double) ); // this is sort of okay I guess (it has desired effect)
  memset( g1, 255, R*sizeof(double) );
  memset( g2, 255, R*sizeof(double) );
  memset( groups, 0, R*sizeof(int) );

  affected_index_size = 0;

  while( getNextFamily( start, end ) ){
    data_num_families++;

    // get the alleles that correspond to this family
    int numChild = 0;
    int ca[100], cb[100]; // the childs alleles
    int childi[100]; // indexing into the children
    int curParent = 0;
    int p[4]; // the parents alleles
    memset( p, 0, sizeof(int)*4 );
    for( int i=start; i<=end; i++ ){
      if( elt(i,C_FATH)==0 && elt(i,C_MOTH)==0 ) {

        // then you know it's a parent -- oy, a requirement!
        if( curParent>1 ) {
          cout << "More than two parents in a pedigree! Current code can only handle nuclear pedigrees where the parents have no parents." << endl;
          cout << "pid=" << elt(i,C_PID) << endl;
          exit(1);
        }

        p[curParent*2]     = (int)elt(i,c_m0);
        p[curParent*2 + 1] = (int)elt(i,c_m1);
        curParent++;
      }else{
        // it's a child
        ca[numChild] = (int)elt(i,c_m0);
        cb[numChild] = (int)elt(i,c_m1);
        childi[numChild] = i;
        numChild++;
      }
    }

    if( numChild==0 ) {
      cout << "No children in pedigree." << endl;
      continue; // don't do any more with it...
    }

    double g[3];
    int group = pG_group( numChild, p, &p[2], ca, cb, g );
    //cout << group << endl; // DEBUG DEBUG DEBUG!!!

    if( group != -1 ) { // this is going to always happen now...
      // then it was informative, so worthwhile to put info in it

      // loop through all the children
      for( int c=0; c<numChild; c++ ){
        // set the data stuff
        groups[ childi[c] ] = group;
        // NOTE: It's sort of coded backwards in fbatDist.cpp than
        //  I've been thinking of it in terms of g0 g1 g2
        //  which, stupid, is why you use constants and things that are
        //  a bit more agreeable. Shame on you.
        g0[ childi[c] ] = g[gBB];
        g1[ childi[c] ] = g[gAB];
        g2[ childi[c] ] = g[gAA];

        // and add it to the affected list if they are affected
        if( elt(childi[c],C_AFF) == 2 ) {
          // but only use the first affected in a family!
          if( prevInformativePid!=elt(childi[c],C_PID) ) {
            affected_index[ affected_index_size ] = childi[c];
            affected_index_size++;

            prevInformativePid = (int)elt(childi[c],C_PID); // addition for sp1a here...
          }
        }
      }
    }
  }

}

void DataMatrix::genPush( int pid, int id, int idfath, int idmoth,
                          int sex, int affection,
                          int m0, int m1,
                          double env,
                          int curRow
                        )
{
  //cout << curRow << ": " << pid << " " << id << " " << idfath << " " << idmoth << " " << sex << " " << affection << " " << m0 << " " << m1 << " " << env << endl;

  elt( curRow, C_PID )  = pid;
  elt( curRow, C_ID )   = id;
  elt( curRow, C_FATH ) = idfath;
  elt( curRow, C_MOTH ) = idmoth;
  elt( curRow, C_SEX )  = sex;
  elt( curRow, C_AFF )  = affection;
  elt( curRow, c_m0 )   = m0;
  elt( curRow, c_m1 )   = m1;
  elt( curRow, c_env )  = env;
}


// R routine
extern "C" {
  void dataComputeGroupG( double *data, int *dataDim,
                          int *m0pos, int *m1pos,
                          int *groups,
                          double *g0, double *g1, double *g2,
                          int *affected_index,
                          int *affected_index_size,
                          int *data_num_families )
  {
    DataMatrix m;
    m.set( data, dataDim );
    m.c_m0 = *m0pos;
    m.c_m1 = *m1pos;

    int p_affected_index_size;
    int p_data_num_families;
    m.computeGroupG( groups,  g0, g1, g2,  affected_index,
                     p_affected_index_size, p_data_num_families );
    *affected_index_size = p_affected_index_size;
    *data_num_families = p_data_num_families;
  }
}
