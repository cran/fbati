// DOES NOT HANDLE MISSING GENOTYPE OR PHENOTYPE INFORMATION PROPERLY YET!!!

// g++ -I/usr/share/R/include      -fpic  -O3 -pipe  -g -c fbatge.cpp -o fbatge.o -pedantic


// ? Does this code enforce that the extra covariates actually exist?

// The usual includes
#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
//#include <math.h>
#include <R.h>
using namespace std;

const double ZEROTOL = 1.490116e-08;


/////////////////////////////
// MISCELLANEOUS FUNCTIONS //
/////////////////////////////

// a more numerically accurate summing routine -- does this really matter, or not?
/*
double sum( double* v, int LEN )
{
  // DEBUG ONLY
  //double standardSum = 0.0;
  //for( int i=0; i<LEN; i++ )
  //  standardSum += v[i];
  // YLNO GUBED

  // some special cases
  if( LEN == 0 ) {
    return( 0.0 );
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

  //cout << "Stable sum = " << v[0] << ", Standard sum = " << standardSum << endl; // DEBUG ONLY

  return( v[0] );
}*/

// converts a double (or integer, really) to a string
string d2s( double d ) {
  ostringstream oss;
  oss << d;
  return( oss.str() );
}//d2s -- DEBUGGED

// Returns all possible permutations of the data
void perms( vector<int> &data, vector< vector<int> > &ret ) {
  if( data.size() == 1 ) {
    ret.push_back( data );
    return;
  }

  for( unsigned int i=0; i<data.size(); i++ ) {
    vector<int> subdata = data;
    subdata.erase( subdata.begin() + i );
    vector< vector<int> > subret;
    perms( subdata, subret );
    for( unsigned int s=0; s<subret.size(); s++ ) {
      subret[s].push_back(data[i]);
      ret.push_back(subret[s]);
    }//s
  }//i
}//perms -- DEBUGGED

// takes current permutations, and adds all possibilities of next...
void fanperms( vector<int> &next, vector< vector<int> > &perm ) {
  unsigned int P = perm.size();

  // special initial case
  if( perm.size() == 0 ) {
    perm.resize( next.size() );
    for( unsigned int n=0; n<next.size(); n++ )
      perm[n].push_back( next[n] );
    return;
  }

  // duplicate the rows in perm
  for( unsigned int n=1; n<next.size(); n++ ) { // note starting at 1, this is correct, as first is already there
    for( unsigned int p=0; p<P; p++ )
      perm.push_back( perm[p] );
  }

  // now push on the next;
  for( unsigned int n=0; n<next.size(); n++ )
    for( unsigned int p=0; p<P; p++ )
      perm[ n*P + p ].push_back( next[n] );
}//fanperms -- DEBUGGED

// like fanperms, but can be weighted now
void fanpermsw( vector<int> &nextPerm, vector<double> &nextWeight, vector< vector<int> > &perm, vector<double> &weight ) {
  unsigned int P = perm.size();

  if( nextPerm.size() != nextWeight.size() ) {
    Rprintf("fanpermsw Error, nextPerm.size() != nextWeight.size()\n");
    return;
  }

  // special initial case
  if( perm.size() == 0 ) {
    perm.resize( nextPerm.size() );
    weight.resize( nextWeight.size() );
    for( unsigned int n=0; n<nextPerm.size(); n++ ) {
      perm[n].push_back( nextPerm[n] );
      weight[n] = nextWeight[n];
    }
    return;
  }

  // duplicate the rows in perm and weight
  for( unsigned int n=1; n<nextPerm.size(); n++ ) {
    for( unsigned int p=0; p<P; p++ ) {
      perm.push_back( perm[p] );
      weight.push_back( weight[p] );
    }
  }

  // now push on the next, and multiply the weights
  for( unsigned int n=0; n<nextPerm.size(); n++ ) {
    for( unsigned int p=0; p<P; p++ ) {
      perm[ n*P + p ].push_back( nextPerm[n] );
      weight[ n*P + p ] *= nextWeight[n];
    }
  }
}//fanpermsw -- debugged

// initialSize -- when nonzero, the size of curPerm, and curPerm will be initialized to catFill
void perm2categories( vector< vector<int> > &perm, vector<int> &curPerm, int curLoc, int nPlace, int catFill, int catBase, int initialSize=0 ) {
  if( initialSize!=0 ) {
    curPerm.resize(initialSize);
    for( int p=0; p<initialSize; p++ )
      curPerm[p] = catBase;

    if( nPlace < 1 ) {
      perm.push_back( curPerm );
      return;
    }
  }

  if( nPlace < 1 ) {
    Rprintf("perm2categories error, nPlace<1 (%d) when it should not be.\n", nPlace);
    return;
  }

  for( unsigned int l=curLoc; l<curPerm.size()-nPlace+1; l++ ) {
    curPerm[l] = catFill;
    if( nPlace == 1 ) {
      // then push it on
      perm.push_back( curPerm );
    }else{
      // recurse
      perm2categories( perm, curPerm, l+1, nPlace-1, catFill, catBase );
    }
    curPerm[l] = catBase;
  }
}//perm2categories -- debugged

// For printing the permutation matrices (when debugging)
void printperms( vector< vector<int> > &perm ) {
  for( unsigned int i=0; i<perm.size(); i++ ) {
    for( unsigned int j=0; j<perm[i].size(); j++ )
      Rprintf("%d ", perm[i][j]);
    Rprintf("\n");
  }
}//printperms -- DEBUGGED
void printpermsw( vector< vector<int> > &perm, vector<double> &w ) {
  if( perm.size() != w.size() ) {
    Rprintf("printpermsw perm.size()=%d, but w.size()=%d\n", perm.size(), w.size());
  }

  for( unsigned int i=0; i<perm.size(); i++ ) {
    for( unsigned int j=0; j<perm[i].size(); j++ )
      Rprintf("%d ", perm[i][j]);
    Rprintf("%d\n", w[i]);
  }
}//printpermsw -- DEBUGGED

////////////////////////
// Basic matrix class //
////////////////////////
class MMatrix{
public:
  vector< vector<double> > m;

  void resize( int nrows, int ncols ) {
    m.resize( nrows );
    for( int r=0; r<nrows; r++ )
      m[r].resize( ncols );
  }//MMatrix::resize -- SIMPLE

  // creates from a vector V V^T
  void createVVt( double *v, int vlength ) {
    resize( vlength, vlength );
    for( int r=0; r<nrows(); r++ )
      for( int c=0; c<ncols(); c++ )
        m[r][c] = v[r]*v[c];
  }//MMatrix::createVVt -- DEBUGGED
  // creates from two vectors V1 V2^T
  void createV1V2t( double *v1, int v1length, double *v2, int v2length ) {
    resize( v1length, v2length );
    for( int r=0; r<v1length; r++ )
      for( int c=0; c<v2length; c++ )
        m[r][c] = v1[r] * v2[c];
  }//MMatrix::createV1V2t -- DEBUGGED
  // makes a vector into a matrix, so we can use these routines...
  void createV( double *v, int vlength, bool transpose=false ) {
    if( !transpose ) {
      resize( vlength, 1 );
      for( int i=0; i<vlength; i++ )
        m[i][0] = v[i];
    }else{
      resize( 1, vlength );
      for( int i=0; i<vlength; i++ )
        m[0][i] = v[i];
    }
  }//MMatrix::creatV -- DEBUGGED

  int nrows() {
    return( m.size() );
  }//MMatrix::nrows -- SIMPLE

  int ncols() {
    if( m.size() > 0 )
      return( m[0].size() );
    return( 0 );
  }//MMatrix::ncols -- SIMPLE

  string pad_d2s( double d, int pad=6 ) {
    string s = d2s(d);
    while( (int)s.length() < pad )
      s+=' '; // kind of inefficient, but only for debug, doesn't matter
    return( s );
  }//MMatrix::pad_d2s -- SIMPLE

  string toString() {
    string s;
    for( int r=0; r<nrows(); r++ ) {
      for( int c=0; c<ncols(); c++ )
        s += pad_d2s(m[r][c]) + " ";
      s += "\n";
    }//r
    return( s );
  }//MMatrix::toString -- DEBUGGED

  void addSelf( MMatrix &rhs ) {
    // verify some preconditions for adding two matrices
    if( nrows()!=rhs.nrows() || ncols()!=rhs.ncols() ) {
      Rprintf("MMatrix::add() -- LHS rows=%d != RHS rows=%d OR LHS cols=%d != RHS cols=%d\n", nrows(), rhs.nrows(), ncols(), rhs.ncols());
      //exit(1);
      return;
    }//fi

    // then add the rhs to this matrix
    for( int r=0; r<nrows(); r++ )
      for( int c=0; c<ncols(); c++ )
        m[r][c] += rhs.m[r][c];
  }//MMatrix::add - DEBUGGED
  void addSelfC( MMatrix rhs ) { addSelf( rhs ); }
  MMatrix add( MMatrix &rhs ) {
    MMatrix out = *this;
    out.addSelf( rhs );
    return( out );
  }//MMatrix::addSelfC -- SIMPLE

  void subtractSelf( MMatrix &rhs ) {
    // verify some preconditions
    if( nrows()!=rhs.nrows() || ncols()!=rhs.ncols() ) {
      Rprintf("MMatrix::substractSelf() -- LHS rows=%d != RHS rows=%d OR LHS cols=%d != RHS cols=%d\n", nrows(), rhs.nrows(), ncols(), rhs.ncols());
      //exit(1);
      return;
    }//fi

    // then subtract the rhs from this matrix
    for( int r=0; r<nrows(); r++ )
      for( int c=0; c<ncols(); c++ )
        m[r][c] -= rhs.m[r][c];
  }//MMatrix::subtractSelf -- SIMPLE
  void subtractSelfC( MMatrix rhs ) { subtractSelf( rhs ); }
  MMatrix subtractC( MMatrix rhs ) {
    MMatrix out = *this;
    out.subtractSelf( rhs );
    return( out );
  }//MMatrix::subtractSelfC -- SIMPLE

  void multiply( MMatrix &rhs, MMatrix &out ) {
    // verify preconditions
    if( ncols() != rhs.nrows() ) {
      Rprintf("MMatrix::multiply -- LHS ncols=%d != RHS nrows=%d\n", ncols(), rhs.nrows());
      //exit(1);
      return;
    }//fi

    // then do the multiplication
    out.resize( nrows(), rhs.ncols() );
    for( int r=0; r<nrows(); r++ ) {
      for( int c=0; c<rhs.ncols(); c++ ) {
        out.m[r][c] = 0;
        for( int k=0; k<ncols(); k++ )
          out.m[r][c] += m[r][k] * rhs.m[k][c];
      }//c
    }//r
  }//MMatrix::multiply -- DEBUGGED
  MMatrix multiply( MMatrix &rhs ) {
    MMatrix out;
    multiply( rhs, out );
    return( out );
  }//MMatrix::multiply overload -- SIMPLE
  MMatrix multiplyC( MMatrix rhs ) {
    return( multiply( rhs ) );
  }//MMatrix::multiplyC -- SIMPLE

  void multiplySelf( double scalar ) {
    for( int r=0; r<nrows(); r++ )
      for( int c=0; c<ncols(); c++ )
        m[r][c] *= scalar;
  }//MMatrix::multiflySelf -- SIMPLE (scalar)
  MMatrix multiply( double scalar ) {
    MMatrix out = *this;
    out.multiplySelf( scalar );
    return( out );
  }//MMatrix::multiply -- SIMPLE (scalar)

  void transpose( MMatrix &out ) {
    out.resize( ncols(), nrows() );
    for( int r=0; r<nrows(); r++ )
      for( int c=0; c<ncols(); c++ )
        out.m[c][r] = m[r][c];
  }//MMatrix::transpose -- DEBUGGED
  MMatrix transpose() {
    MMatrix out;
    transpose( out );
    return( out );
  }//MMatrix::transpose overload -- SIMPLE

  // fills in a matrix (i.e. to zero it)
  void fill( double value ) {
    for( int r=0; r<nrows(); r++ )
      for( int c=0; c<ncols(); c++ )
        m[r][c] = value;
  }//MMatrix::fill -- SIMPLE

  MMatrix inv2x2() {
    MMatrix out;
    if( nrows()!=2 || ncols()!= 2 ) {
      Rprintf("MMatrix::inv2x2, not a 2x2 matrix! Dimensions: %d, %d\n", nrows(), ncols());
      //exit(1);
      return(out);
    }

    out.resize( 2, 2 );

    double a=m[0][0], b=m[0][1], c=m[1][0], d=m[1][1];
    double div = a*d - b*c;
    out.m[0][0] = d / div;
    out.m[0][1] = - b / div;
    out.m[1][0] = - c / div;
    out.m[1][1] = a / div;

    return( out );
  }//MMatrix::inv2x2 - DEBUGGED

  MMatrix subMatrix( int rowStart, int rowEnd, int colStart, int colEnd ) {
    MMatrix out;
    if( rowStart<0 || rowEnd>nrows()-1 || colStart<0 || colEnd>ncols()-1 || rowStart>rowEnd || colStart>colEnd ) {
      Rprintf("MMatrix::subMatrix invalid dimensions supplied (rowStart=%d, rowEnd=%d, colStart=%d, colEnd=%d), the dimensions of the matrix are %dx%d\n", rowStart, rowEnd, colStart, colEnd, nrows(), ncols());
      //exit(1);
      return(out);
    }

    out.resize( rowEnd-rowStart+1, colEnd-colStart+1 );
    for( int r=0; r<out.nrows(); r++ )
      for( int c=0; c<out.ncols(); c++ )
        out.m[r][c] = m[rowStart+r][colStart+c];

    return( out );
  }//MMatrix subMatrix -- SEMI-SIMPLE

  // fills in the matrix easily for use with the debug routine
  void debugFill( int nrows, int ncols ) {
    resize( nrows, ncols );
    for( int r=0; r<nrows; r++ )
      for( int c=0; c<ncols; c++ )
        m[r][c] = r*100 + c;
  }//MMatrix::debugFill -- DEBUGGED

  // Tests all major operations
  static void debug() {
    MMatrix m1, m2, m3, out; // we'll add m1 and m2, mult m2 and m3, and transpose m2
    m1.debugFill( 3, 2 );
    m2.debugFill( 3, 2 );
    m3.debugFill( 2, 3 );

    Rprintf("ADDING\n%s\nAND\n%s\nYIELDS\n", m1.toString().c_str(), m2.toString().c_str());
    m1.addSelf( m2 );
    Rprintf("%s\n\n", m1.toString().c_str());

    Rprintf("MULTIPLYING\n%s\nAND\n%s\nYIELDS\n", m2.toString().c_str(), m3.toString().c_str());
    m2.multiply( m3, out );
    Rprintf("%s\n\n", out.toString().c_str());

    Rprintf("TRANSPOSE OF\n%s\nYIELDS\n", m2.toString().c_str());
    m2.transpose( out );
    Rprintf("%s\n\n", out.toString().c_str());

    Rprintf("SECOND SET DEBUGGING\n\n");
    double v1[] = {1,2,3,4};
    double v2[] = {101,102,103,104,105,106};
    MMatrix v1m, v2mt, v1v1t, v1v2t;
    v1m.createV( v1, 4 );
    v2mt.createV(v2, 6, true );
    v1v1t.createVVt( v1, 4 );
    v1v2t.createV1V2t( v1, 4, v2, 6 );
    Rprintf("V1\n%s\n", v1m.toString().c_str());
    Rprintf("V2\n%s\n", v2mt.toString().c_str());
    Rprintf("V1V1t\n%s\n", v1v1t.toString().c_str());
    Rprintf("v1v2t\n%s\n", v1v2t.toString().c_str());

    MMatrix tt;
    tt.resize( 2, 2 );
    tt.m[0][0] = 4; tt.m[1][0] = 12; tt.m[1][0] = 13; tt.m[1][1] = 44;
    Rprintf("tt\n%s\n", tt.toString().c_str());
    Rprintf("inv(tt)\n%s\n", tt.inv2x2().toString().c_str());
    Rprintf("tt inv(tt)\n%s\n", tt.multiplyC( tt.inv2x2() ).toString().c_str());
  }//MMatrix::debug -- PASSES
};//MMatrix -- DEBUGGED!

extern "C" {
  void cpp_mmatrix_debug() {
    MMatrix::debug();
  }
}


////////////////////////////////////
// Random number generating class //
////////////////////////////////////
class Random{
public:
  static const int MAX_TRY = 10000;
  vector< vector<double> > cholesky; // for normal

  /** attaching and detaching **/
  void attach() { GetRNGstate(); }
  void detach() { PutRNGstate(); }

  /** Normal routines **/
  // chol is the cholesky decomposition of sigma (from R)
  // p is the dimensions of sigma
  void setNormalSigma( double *chol, int p ) {
    cholesky.resize(p);
    for( int r=0; r<p; r++ )
      cholesky[r].resize(p);

    // order starts by going down all the rows in the first column
    int rc=0;
    for( int c=0; c<p; c++ ) {
      for( int r=0; r<p; r++ ) {
        cholesky[r][c] = chol[rc];
        //cout << chol[rc] << endl;
        rc++;
      }//r
    }//c
  }//Random::setNormalSigma -- DEBUGGED

  int P() {
    return cholesky.size();
  }//Random::P -- SIMPLE

  double rnorm() {
    return( norm_rand() );
  }//Random::rnorm -- SIMPLE

  double rnormTrunc( double trunc ) {
    double z;
    for( int t=0; t<MAX_TRY; t++ ) {
      z = norm_rand();
      if( z > -trunc && z < trunc )
        return( z );
    }
    //cout << "Random::rnormTrunc could not ascertain a truncated normal (passed maximal number of tries=" << MAX_TRY << "), continuing anyway..." << endl;
    Rprintf("Random::rnormTrunc count not ascertain a truncated normal (passed maximal number of tries=%d), continuing anyway...\n", MAX_TRY);
    return( z );
  }

  double rnormExt( double mean, double sigma ) {
    return( mean + norm_rand()*sigma );
  }//Random::rnorm -- SIMPLE

  void mvrnorm( vector<double> &v ) {
    // resize if necessary
    if( (int)v.size() != P() )
      v.resize( P() );

    // draw up the univariate random normals
    vector<double> z;
    z.resize( P() );
    for( int i=0; i<P(); i++ )
      z[i] = rnorm();

    // and then multiply by the cholesky decomposition
    for( int i=0; i<P(); i++ ) {
      v[i] = 0.0;
      for( int j=0; j<P(); j++ )
        v[i] += z[j] * cholesky[j][i];
    }//i
  }//Random::mvrnorm -- DEBUGGED

  void mvrnormTrunc( vector<double> &v, double trunc ) {
    for( int t=0; t<MAX_TRY; t++ ) {
      mvrnorm( v );
      double ok=true;
      for( unsigned int i=0; i<v.size(); i++ )
        ok = ok && v[i] > -trunc && v[i] < trunc;
      if( ok ) return; // all set!
    }//t
    Rprintf("Could not ascertain a truncated normal distribution! Proceeding anyway!!!\n");
  }//Random::mvrnormTrunc -- SEMI-SIMPLE

  /** Other random numbers **/
  double runif() {
    return( unif_rand() );
  }//Random::runif -- SIMPLE

  double runifExt( double min, double max ) {
    return( unif_rand() * (max-min) + min );
  }//Random::runifExt -- SIMPLE

  int randInt( int min, int max ) {
    int range = (max-min) + 1;
    return( (int)(runif() * range) + min );
  }//Random::randInt -- SIMPLE

  // Random::debug
  // Debugging -- really just the mvrnorm
  void debug() {
    if( cholesky.size() == 0 ) {
      Rprintf("You need to call 'setNormalSigma' before calling any multivariate random normal routine.\n");
      return;
    }

    // first, can we generate a random number
    vector<double> v;
    mvrnorm( v );
    for( unsigned int i=0; i<v.size(); i++ )
      Rprintf("%d ", v[i]);
    Rprintf("\n");

    // now generate a whole bunch of random numbers, and then print them in such a way that it can be loaded into R, and the correlation matrix checked...
    int N = 500;
    vector< vector<double> > norms;
    for( int n=0; n<N; n++ ) {
      mvrnorm( v );
      norms.push_back( v );
    }//n
    // create a .R file that can be sourced (or vanilla-ed) into R!
    ofstream of;
    of.open( "KILLME_rn_debug.R", ios::out );
    of << "x <- cbind( ";
    for( int p=0; p<P(); p++ ) {
      if( p != 0 ) of << ", ";
      of << "c(";
      for( int n=0; n<N; n++ ) {
        if( n != 0 ) of << ", ";
        of << norms[n][p];
      }//n
      of << ")";
    }//p
    of << ")" << endl;
    of << "cor(x)" << endl;
    of.close();
  }//Random::debug -- DEBUGGED (random multivariate normal only)
};

Random rn;

// RANDOM EXPORTED FUNCTIONS
// - attach before using any routines that need a random number in c++,
//    but then detach before using any random numbers back in R
extern "C" {
  void cpp_rn_attach() { rn.attach(); }
  void cpp_rn_detach() { rn.detach(); }
  void cpp_rn_setNormalSigma( double *chol, int *p ) { rn.setNormalSigma( chol, *p ); }
  void cpp_rn_debug() { rn.debug(); }
}

///////////////////
// GFamily class //
///////////////////
class GFamily{
public:
  static const int GMISS = -1;
  static const int PMISS = -1;

  int parentGeno[2]; // -1=missing here, NA in R code
  vector<int> childGeno;
  vector<int> childTrait; // int so can use perms routine, rather than bool...
  vector<double> childEnvironment;
  vector< vector<double> > childCovariate;

  int numChild() { return(childGeno.size()); }
  int numCovariate() {
    if( childCovariate.size() < 1 ) return(0);
    return( childCovariate[0].size() );
  }//GFamily::numCovariate -- SIMPLE

  vector< vector<int> > genoPerm;
  vector<double> genoPermWeight;
  vector< vector<int> > phenoPerm;
  vector<double> phenoPermWeight;

  // GFamily::clear -- SIMPLE
  void clear() {
    parentGeno[0] = parentGeno[1] = GMISS;
    childGeno.clear();
    childTrait.clear();
    childEnvironment.clear();
    childCovariate.clear();

    genoPerm.clear();
    genoPermWeight.clear();
    phenoPerm.clear();
    phenoPermWeight.clear();
  }//GFamily::clear -- SIMPLE

  // GFamily::toString
  string toString( bool extended=false ) {
    string s;
    s+="Parents g="; s+=d2s(parentGeno[0]); s+=","; s+=d2s(parentGeno[1]); s+="\n";
    for( int i=0; i<numChild(); i++ ) {
      s+=" Child ";
      s+=d2s(i); s+=" (g="; s+=d2s(childGeno[i]); s+=", trait="; s+=d2s(childTrait[i]); s+=", env="; s+=d2s(childEnvironment[i]);
      if( numCovariate() > 0 ) {
        s+=", covariate={";
        for( int c=0; c<numCovariate(); c++ )
          s+=d2s(childCovariate[i][c])+",";
        s+="}";
      }
      s+=")\n";
    }//i

    if( extended ) {
      // prints out the genoPerm (+weight) and phenoPerm (+weight)
      s+=" genoPerm = ";
      for( unsigned int gp=0; gp<genoPerm.size(); gp++ ) {
        for( unsigned int g=0; g<genoPerm[gp].size(); g++ )
          s+=d2s( genoPerm[gp][g] );
        s+="("; s+=d2s( genoPermWeight[gp] ); s+=") ";
      }//gp
      s+="\n";
      s+=" phenoPerm = ";
      for( unsigned int pp=0; pp<phenoPerm.size(); pp++ ) {
        for( unsigned int p=0; p<phenoPerm[pp].size(); p++ )
          s+=d2s( phenoPerm[pp][p] );
        s+="("; s+=d2s( phenoPermWeight[pp] ); s+=") ";
      }
      s+="\n";
    }//if(extended)

    return( s );
  }//GFamily::toString -- DEBUGGED

  // GFamily::normalizeGenoPerm
  void normalizeGenoPerm( bool size=false ) {
    int P = genoPerm.size();
    if( size ) {
      // need to create, should all be equal
      double w = 1.0 / (double)P;
      genoPermWeight.resize(P);
      for( int p=0; p<P; p++ )
        genoPermWeight[p] = w;
    }else{
      if( P != (int)genoPermWeight.size() ) {
        Rprintf("GFamily::normalizeGenoPerm error, genoPermWeight.size()=%d, but genoPerm.size()=%d.", genoPermWeight.size(), genoPerm.size());
        return;
      }

      double sumw = 0.0;
      for( int p=0; p<P; p++ )
        sumw += genoPermWeight[p];
      for( int p=0; p<P; p++ )
        genoPermWeight[p] /= sumw;
    }
  }//GFamily::normalizeGenoPerm -- DEBUGGED

  // GFamily::setGenoPerm
  void setGenoPerm() {
    // borrows from Rabinowitz and Laird
    // -- is this really what we want to do?
    //perms( childGeno, genoPerm ); // will do some reduntants...? I don't think this is enough
    genoPerm.clear();
    genoPermWeight.clear();

    // precursor -- non of the children can have missing genotypes
    for( int c=0; c<numChild(); c++ )
      if( childGeno[c] == GMISS )
        Rprintf("GFamily::setGenoPerm() cannot handle when there is missing genotype information in the offspring.\n");

    // first sort the parents
    if( parentGeno[0] > parentGeno[1] ) {
      int temp = parentGeno[0];
      parentGeno[0] = parentGeno[1];
      parentGeno[1] = temp;
    }

    if( parentGeno[0] != GMISS ) {
      // parents are both present (remember, sorted, and GMISS = -1)

      if(  ( parentGeno[0]==0 && parentGeno[1]==0 )  ||  ( parentGeno[0]==2 && parentGeno[1]==2 )  ||  ( parentGeno[0]==0 && parentGeno[1]==2 )  ) {
        // then they are just the observed
        genoPerm.push_back( childGeno );
        genoPermWeight.push_back( 1.0 );
      }else{
        vector<int> g;    // possible genotypes
        vector<double> w; // weight of the genotypes

        if( parentGeno[0]==1 && parentGeno[1]==1 ) {
          g.push_back( 0 ); w.push_back( 0.25 );
          g.push_back( 1 ); w.push_back( 0.5  );
          g.push_back( 2 ); w.push_back( 0.25 );
        }else if( parentGeno[0]==0 && parentGeno[1]==1 ) {
          g.push_back( 0 ); w.push_back( 0.5 );
          g.push_back( 1 ); w.push_back( 0.5 );
        }else if( parentGeno[0]==1 && parentGeno[1]==2 ) {
          g.push_back( 1 ); w.push_back( 0.5 );
          g.push_back( 2 ); w.push_back( 0.5 );
        }else{
          Rprintf("GFamily::setGenoPerm() family with parents fell outside of all cases, parentGeno[0]=%d, parentGeno[1]=%d\n", parentGeno[1], parentGeno[0]);
        }

        for( int c=0; c<numChild(); c++ )
          fanpermsw( g, w, genoPerm, genoPermWeight );
      }
    }else{//parentGeno[0] == GMISS
      // one or both parents are missing

      // time to whip out our Rabinowitz and Laird

      // count the number of each genotype
      int n[3] = {0,0,0};
      for( int c=0; c<numChild(); c++ )
        n[ childGeno[c] ]++;

      if(  ( n[0]==0 && n[1]==0 )  ||  ( n[0]==0 && n[2]==0 )  ||  ( n[1]==0 && n[2]==0 )  ) {
        // only one type of genotype
        // Case 1, R&L, {AA} or {AB}
        genoPerm.push_back( childGeno );
        genoPermWeight.push_back( 1.0 );
      }else if( n[0]!=0 && n[1]!=0 && n[2]==0 ) {
        // Case 2, R&L, {AA, AB}
        // random assignment of AA and AB that keeps invariant the number of both
        vector<int> curPerm; // just needs to be allocated...
        perm2categories( genoPerm, curPerm, 0, n[1], 1, 0, numChild() );
        normalizeGenoPerm( true );
      }else if( n[0]==0 && n[1]!=0 && n[2]!=0 ) {
        // Case 2, R&L, reversed {BB,AB}
        vector<int> curPerm;
        perm2categories( genoPerm, curPerm, 0, n[2], 2, 1, numChild() );
        normalizeGenoPerm( true );
      }else if( n[0]!=0 && n[2]!=0 ) {
        // Case 3, R&L, {AA,BB} or {AA,AB,BB}
        // randomly assign AA, BB, and AB with probabilities 1/2, 1/4, 1/2, independently to each sib, and discard outcomes without at least one assignment of AA and one assignment of BB

        // first randomly assign
        vector<int> g;    // possible genotypes
        vector<double> w; // weight of the genotypes

        g.push_back( 0 ); w.push_back( 0.25 );
        g.push_back( 1 ); w.push_back( 0.5  );
        g.push_back( 2 ); w.push_back( 0.25 );

        vector< vector<int> > genoPermRaw; // contains some that will not pass
        vector<double> genoPermWeightRaw; // contains some that will not pass
        for( int c=0; c<numChild(); c++ )
          fanpermsw( g, w, genoPermRaw, genoPermWeightRaw );

        // second, remove
        for( unsigned int p=0; p<genoPermRaw.size(); p++ ) {
          bool has0=false, has2=false;
          for( unsigned int g=0; g<genoPermRaw[p].size(); g++ ) {
            if( genoPermRaw[p][g]==0 ) {
              has0 = true;
            }else if( genoPermRaw[p][g]==2 ) {
              has2 = true;
            }
          }//g
          if( has0 && has2 ) {
            // then it's a valid one
            genoPerm.push_back( genoPermRaw[p] );
            genoPermWeight.push_back( genoPermWeightRaw[p] );
          }//if( has0 && has2 )
        }//p

        normalizeGenoPerm();
      }//if( n[0]==0 ... )
    }//if( parentGeno[0] != GMISS )
  }//GFamily::setGenoPerm -- DEBUGGED

  // GFamily::setGenoPermObserved
  void setGenoPermObserved() {
    genoPerm.clear();
    genoPermWeight.clear();

    genoPerm.push_back( childGeno );
    genoPermWeight.push_back( 1 );
  }//GFamily::setGenoPermObserved -- SIMPLE

  // GFamily::normalizePhenoPerm
  void normalizePhenoPerm( bool size=false ) {
    int P = phenoPerm.size();
    if( size ) {
      // need to create, should all be equal
      double w = 1.0 / (double)P;
      phenoPermWeight.resize(P);
      for( int p=0; p<P; p++ )
        phenoPermWeight[p] = w;
    }else{
      if( P != (int)phenoPermWeight.size() ) {
        Rprintf("GFamily::normalizePhenoPerm error, phenoPermWeight.size()=%d, but phenoPerm.size()=%d\n", phenoPermWeight.size(), phenoPerm.size());
        return;
      }

      // really need to normalize
      double sumw = 0.0;
      for( int p=0; p<P; p++ )
        sumw += phenoPermWeight[p];
      for( int p=0; p<P; p++ )
        phenoPermWeight[p] /= sumw;
    }
  }//GFamily::normalizePhenoPerm -- DEBUGGED

  // GFamily::setPhenoPerm
  void setPhenoPerm() {
    // Precursor -- none of the phenotypes can be missing
    for( int c=0; c<numChild(); c++ )
      if( childTrait[c] == PMISS )
        Rprintf("GFamily::setPhenoPerm() cannot handle when there is missing phenotype information in the offspring.\n");

    phenoPerm.clear();
    phenoPermWeight.clear();

    // creats set of all possible phenotype permutations
    //perms( childTrait, phenoPerm );
    //int P = phenoPerm.size();
    //phenoPermWeight.resize(P);
    //for( int p=0; p<P; p++ )
    //  phenoPermWeight[p] = 1.0; // don't need to normalize, because it's effectively in the denomerator as well...

    // can be done much more efficiently...
    // except what happens with missing phenotypes here? And what happens with missing genotypes in our other code?

    int numAffected = 0;
    for( int j=0; j<numChild(); j++ )
      numAffected += childTrait[j];
    if( numAffected == 0 ) {
      // then observed has probability 1
      phenoPerm.push_back( childTrait );
      phenoPermWeight.push_back(1.0);
      return;
    }
    vector<int> curPerm; // meant to be left empty
    perm2categories( phenoPerm, curPerm,  0, numAffected,  1, 0,  numChild() );

    // and set the weights to be all the same
    normalizePhenoPerm( true );
  }//GFamily::setPhenoPerm -- DEBUGGED

  // GFamily::setPhenoPermObserved
  void setPhenoPermObserved() {
    phenoPerm.clear();
    phenoPermWeight.clear();

    phenoPerm.push_back( childTrait );
    phenoPermWeight.push_back( 1 );
  }//GFamily::setPhenoPermObserved -- SIMPLE

  // GFamily::debugPermAddFamily
  void debugPermAddFamily( int pa, int pb, // the parents genotypes
                           int n0, int n1, int n2, // number of each offspring genotype
                           int a0, int a1, int a2, // affection of each offspring genotype
                           bool test=true ) {
    if( a0>n0 || a1>n1 || a2>n2 || a0<0 || a1<0 || a2<0 )
      Rprintf("GFamily::debugPermAddFamily invalid input.");

    clear(); // make sure everything is cleared first!

    parentGeno[0] = pa;
    parentGeno[1] = pb;

    int n[3] = {n0, n1, n2};
    int a[3] = {a0, a1, a2};

    for( int j=0; j<3; j++ ) {
      for( int i=0; i<n[j]; i++ ) {
        childGeno.push_back(j);
        if( i<a[j] ) {
          childTrait.push_back(1);
        }else{
          childTrait.push_back(0);
        }
      }
    }

    // now do the permutations, and print!
    if( test ) {
      setGenoPerm();
      setPhenoPerm();
      Rprintf("%s\n", toString(true).c_str());
    }
  }//GFamily::debugPermAddFamily -- DEBUGGED

};//GFamily


////////////////
// GPed class //
////////////////
class GPed{
public:
  vector<GFamily> families;
  enum STRATEGY { GENO, PHENO, ADAPTIVE };
  STRATEGY strategy;

  int numFamilies() { return( families.size() ); }
  void resize( unsigned int size ) { families.resize( size ); }

  // Handle communication because we need the generalized inverse of dbnuis,
  //  so store in statPrecompute, and then finish up in statCompute
  MMatrix dbge, dbnuis; // Actually averages of these derivatives
  MMatrix u; // Need to be global variables
  int M;

  // for the coding of the parameters...
  static const int X0 = 2;
  static const int X1 = 1;

  // GPed::clear
  // clear out all families
  void clear(){
    families.clear();
  }//GPed::clear -- SIMPLE

  // GPed::toString
  string toString( bool extended=false ) {
    string str;
    for( int f=0; f<numFamilies(); f++ )
      str += families[f].toString(extended); // + "\n";
    return( str );
  }//GPed::toString -- DEBUGGED

  // GPed::numCovariates
  int numCovariates() {
    for( unsigned int f=0; f<families.size(); f++ ) {
      if( families[f].numChild() > 0 ) {
        //cout <<  families[f].numCovariate() << endl;
        return( families[f].numCovariate() );
      }
    }
    return( 0 );
  }

  // GPed::fillPerms -- assumes strategy has been set...
  void fillPerms() {
    switch( strategy ) {
    case GENO:
      for( int i=0; i<numFamilies(); i++ ) {
        families[i].setGenoPerm();
        families[i].setPhenoPermObserved();
      }
      break;
    case PHENO:
      for( int i=0; i<numFamilies(); i++ ) {
        families[i].setGenoPermObserved();
        families[i].setPhenoPerm();
      }
      break;
    case ADAPTIVE:
      for( int i=0; i<numFamilies(); i++ ) {
        families[i].setGenoPerm();
        families[i].setPhenoPerm();
      }
      break;
    default:
      Rprintf("Strategy %d has not been enumerated. Likely that 'strategy' was not set before calling, or, far worse, memory is being overwritten.\n", (int)strategy);
      //exit(1);
      return;
      break;
    }

    // make sure everyone has something
    for( int i=0; i<numFamilies(); i++ ) {
      if( families[i].genoPerm.size()==0 || families[i].phenoPerm.size()==0 ) {
        Rprintf("genoPerm or phenoPerm left completely empty - should at least have the observed in it.\n");
        //exit(1);
        return;
      }
    }
  }//GPed::fillPerms -- SIMPLE

  // GPed::setStrategy
  void setStrategy( char *strat ) {
    string strategyStr = strat;

    // Set the strategy
    if( strategyStr == "geno" ) {
      strategy = GENO;
    }else if( strategyStr == "pheno" ) {
      strategy = PHENO;
    }else if( strategyStr == "adaptive" ) {
      strategy = ADAPTIVE;
    }else{
      Rprintf("GPed::setStrategy not understood, should be 'geno', 'pheno', or 'adaptive'; you supplied ' %d.\n", strat);
      //exit(1);
      return;
    }

    // AND CALL fillPerms!
    fillPerms();
  }//GPed::setStrategy -- SIMPLE

  // Unified estimating equation (works for all three approaches!)
  // GPed::estEq
  void estEq( double *beta, int betaLength,
              MMatrix &ui ) {
    // beta = bge (2), bg (2), be (1), bc (J) [[b0 is not modelled]]
    if( betaLength < 4 ) {
      // verify this precondition
      Rprintf("GPed::estEq(...) betaLength=%d, but it must be at least of length 4.\n", betaLength);
      //exit(1);
      return;
    }

    ui.resize( numFamilies(), betaLength );
    ui.fill( 0 ); // zero out ui

    // loop over families
    for( int i=0; i<numFamilies(); i++ ) {
      //double m[betaLength]; // i.e. the observed
      //double numSum[betaLength];
      vector<double> m; m.resize(betaLength);
      vector<double> numSum; numSum.resize(betaLength);
      double denSum = 0.0; // scalar
      for( int b=0; b<betaLength; b++ )
        m[b] = numSum[b] = 0.0;

      // 1) Compute the observed
      // loop over individuals
      for( int j=0; j<families[i].numChild(); j++ ) {
        if( families[i].childTrait[j] == 1 ) {
          // compute observed
          double x0 = (double)(families[i].childGeno[j]==X0);
          double x1 = (double)(families[i].childGeno[j]==X1);
          double z  = families[i].childEnvironment[j];
          m[0] += x0 * z;
          m[1] += x1 * z;
          m[2] += x0;
          m[3] += x1;
          if( betaLength > 4 )
            m[4] += z;
          if( betaLength > 5 )
            for( int c=5; c<betaLength; c++ )
              m[c] += families[i].childCovariate[j][c-5];
        }
      }

      // 2) Compute the expected
      // loop over geno perms
      for( unsigned int gp=0; gp<families[i].genoPerm.size(); gp++ ) {
        // loop over pheno perms
        for( unsigned int pp=0; pp<families[i].phenoPerm.size(); pp++ ) {
          // add the current to the expected sum
          //double mStar[betaLength];
          vector<double> mStar; mStar.resize(betaLength);
          for( int b=0; b<betaLength; b++ )
            mStar[b] = 0.0;
          // loop over individuals
          for( int j=0; j<families[i].numChild(); j++ ) {
            if( families[i].phenoPerm[pp][j] == 1 ) {
              double x0star = (double)(families[i].genoPerm[gp][j]==X0);
              double x1star = (double)(families[i].genoPerm[gp][j]==X1);
              double zstar  = families[i].childEnvironment[j]; // no perming here
              mStar[0] += x0star * zstar;
              mStar[1] += x1star * zstar;
              mStar[2] += x0star;
              mStar[3] += x1star;
              if( betaLength > 4 )
                mStar[4] += zstar;
              if( betaLength > 5 )
                for( int c=5; c<betaLength; c++ )
                  mStar[c] += families[i].childCovariate[j][c-5];
            }//fi
          }//j
          double weight = families[i].genoPermWeight[gp] * families[i].phenoPermWeight[pp];
          ////cout << "weight = " << weight << endl;
          double expStuff = 0.0;
          for( int b=0; b<betaLength; b++ )
            expStuff += beta[b] * mStar[b];
          expStuff = exp( expStuff ) * weight;

          for( int b=0; b<betaLength; b++ )
            numSum[b] += mStar[b] * expStuff;

          denSum += expStuff;
        }//pp
      }//gp

      //cout << "denSum[" << i << "]=" << denSum << endl;

      for( int b=0; b<betaLength; b++ )
        ui.m[i][b] += m[b] - numSum[b] / denSum;
    }//i

    ////cout << "(c++ ) u = ";
    ////for( int b=0; b<betaLength; b++ )
    ////  cout << u[b] << " ";
    ////cout << endl;
  }//Gped::estEqMatrix -- CRUCIAL ERROR PRONE ROUTINE, NOT DEBUGGED (theory also could have errors)
  void estEq( double *beta, int betaLength,
              double *u ) {
    MMatrix ui;
    estEq( beta, betaLength, ui );

    for( int b=0; b<betaLength; b++ ) {
      u[b] = 0.0;

      // Potentially numerically instable?
      for( int i=0; i<ui.nrows(); i++ )
        u[b] += ui.m[i][b];

      // Nope, doesn't matter
      //double temp[ui.nrows()];
      //for( int i=0; i<ui.nrows(); i++ )
      //  temp[i] = ui.m[i][b];
      //u[b] = sum( temp, ui.nrows() );
    }
  }//GPed::estEq -- overload for C++/R calling, turned out we needed the other for the derivative...

  void statPrecompute( double *beta, int betaLength,
                       double *dbnuis_ret ) { // needs to be inverted in R, and then passed back here...
    // compute the derivatives
    M = betaLength - 2;
    //MMatrix dbge, dbnuis, u; // Need to be global variables
    dbge.resize( 2, M );
    dbnuis.resize( M, M );
    dbge.fill( 0.0 );
    dbnuis.fill( 0.0 );

    // go across the families
    for( int i=0; i<numFamilies(); i++ ) {
      double Ai00 = 0.0;
      MMatrix Ai11, Ai10, Ai01t, Ai02t;
      Ai11.resize( 2, M );
      Ai10.resize( 2, 1 );
      Ai01t.resize( M, 1 ); //( 1, M );
      Ai02t.resize( M, M ); // Transpose of Ai02, a little awkward notation...
      Ai11.fill( 0.0 );
      Ai10.fill( 0.0 );
      Ai01t.fill( 0.0 );
      Ai02t.fill(0.0 );

      // 2) Compute the expected
      // loop over geno perms
      for( unsigned int gp=0; gp<families[i].genoPerm.size(); gp++ ) {
        // loop over pheno perms
        for( unsigned int pp=0; pp<families[i].phenoPerm.size(); pp++ ) {
          // add the current to the expected sum
          //double mStar[betaLength];
          vector<double> mStar; mStar.resize(betaLength);
          for( int b=0; b<betaLength; b++ )
            mStar[b] = 0.0;
          // loop over individuals
          for( int j=0; j<families[i].numChild(); j++ ) {
            if( families[i].phenoPerm[pp][j] == 1 ) {
              double x0star = (double)(families[i].genoPerm[gp][j]==X0);
              double x1star = (double)(families[i].genoPerm[gp][j]==X1);
              double zstar  = families[i].childEnvironment[j]; // no perming here
              mStar[0] += x0star * zstar;
              mStar[1] += x1star * zstar;
              mStar[2] += x0star;
              mStar[3] += x1star;
              if( betaLength > 4 )
                mStar[4] += zstar;
              if( betaLength > 5 )
                for( int c=5; c<betaLength; c++ )
                  mStar[c] += families[i].childCovariate[j][c-5];
            }//fi
          }//j

          double weight = families[i].genoPermWeight[gp] * families[i].phenoPermWeight[pp];

          double expStuff = 0.0;
          for( int b=0; b<betaLength; b++ )
            expStuff += beta[b] * mStar[b];
          expStuff = exp( expStuff ) * weight;

          MMatrix xz;
          xz.createV( (double *)mStar.data(), 2 );
          MMatrix m;
          m.createV( &mStar[2], M );
          MMatrix m_mt;
          m_mt.createVVt( &mStar[2], M );

          Ai00 += expStuff;
          Ai11.addSelfC( xz.multiplyC( m.transpose() ).multiply( expStuff ) );
          Ai10.addSelfC( xz.multiply( expStuff ) );
          Ai01t.addSelfC( m.multiply( expStuff ) );
          Ai02t.addSelfC( m_mt.multiply( expStuff ) );
        }//pp
      }//gp


      // now add these pieces on to the derivative
      dbge.addSelfC( Ai11.multiply(Ai00).subtractC( Ai10.multiplyC(Ai01t.transpose()) ).multiply( -1.0/Ai00/Ai00 ) );
      dbnuis.addSelfC( Ai02t.multiply(Ai00).subtractC( Ai01t.multiplyC( Ai01t.transpose() ) ).multiply( -1/Ai00/Ai00 ) );
    }// i

    // get the ui matrix
    estEq( beta, betaLength, u );

    // copy dbnuis into dbnuis_ret
    for( int r=0; r<dbnuis.nrows(); r++ )
      for( int c=0; c<dbnuis.ncols(); c++ )
        dbnuis_ret[ r + c*M ] = dbnuis.m[r][c]; // should be MxM

    // Now return to R to get the inverse, which will then be passed to statCompute, and we will continue
  }

  void statCompute( double *dbnuis_inv, // M x M
                    double *chisq2df, int *numInf ) { // return value follows chi-squared distribution with 2 degrees of freedom

    // Put the inverse into a MMatrix
    MMatrix dbnuisInv;
    dbnuisInv.resize( M, M );
    for( int r=0; r<M; r++ )
      for( int c=0; c<M; c++ )
        dbnuisInv.m[r][c] = dbnuis_inv[ r + c*M ];
    //cout << "statCompute_A dbnuisInv" << endl << dbnuisInv.toString() << endl;

    MMatrix uge = u.subMatrix( 0, numFamilies()-1, 0, 1 ); // rowstart, rowEnd, colStart, colEnd
    MMatrix ug  = u.subMatrix( 0, numFamilies()-1, 2, M+2-1 );
    //MMatrix w = uge.subtractC( dbge.multiplyC( dbnuisInv ).multiply( ug ) );
    MMatrix w = uge.subtractC( dbge.multiplyC( dbnuisInv ).multiplyC( ug.transpose() ).transpose() );

    MMatrix sumwi, sumwiwit;
    sumwi.resize(2,1);     sumwi.fill( 0.0 );
    sumwiwit.resize(2,2);  sumwiwit.fill( 0.0 );
    for( int i=0; i<numFamilies(); i++ ) {
      sumwi.m[0][0] += w.m[i][0];
      sumwi.m[1][0] += w.m[i][1];

      sumwiwit.m[0][0] += w.m[i][0] * w.m[i][0];
      sumwiwit.m[0][1] += w.m[i][1] * w.m[i][0];
      sumwiwit.m[1][0] += w.m[i][0] * w.m[i][1]; // horribly inefficient
      sumwiwit.m[1][1] += w.m[i][1] * w.m[i][1];
    }//i

    /*
    MMatrix sumwi, sumwiwit, newwi;
    sumwi.resize(2,1);    sumwi.fill( 0.0 );
    sumwiwit.resize(2,2); sumwiwit.fill( 0.0 );
    for( int i=0; i<numFamilies(); i++ ) {
      newwi = w.subMatrix( i, i, 0, 1 ).transpose();

      sumwi.addSelfC( newwi );
      sumwiwit.addSelfC( newwi.multiplyC( newwi.transpose() ) );
    }*/


    // Number of informative families
    *numInf = 0;
    //////*numInf += w.m[i][0]!=0 || w.m[i][1]!=0; // one of them has a contribution
    //for( int i=0; i<numFamilies(); i++ ) // could contribute to uge or to nuisance parameter
    //  *numInf += (w.m[i][0] < -ZEROTOL) || (w.m[i][0] > ZEROTOL) || (w.m[i][1] < -ZEROTOL) || (w.m[i][1] > ZEROTOL);
    for( int i=0; i<numFamilies(); i++ ) // only if contributes to uge (if only to nuisance parameter, it is not counted)
      *numInf += (uge.m[i][0] < -ZEROTOL) || (uge.m[i][0] > ZEROTOL) || (uge.m[i][1] < -ZEROTOL) || (uge.m[i][1] > ZEROTOL);

    if( sumwiwit.m[0][0]*sumwiwit.m[1][1] == sumwiwit.m[1][0] * sumwiwit.m[0][1] ) {
      *chisq2df = 0; // then set the p-value to be 1, cannot be inverted...
    }else{
      *chisq2df = sumwi.transpose().multiplyC( sumwiwit.inv2x2() ).multiply( sumwi ).m[0][0];
    }
    //cout << "CHISQ2DF " << *chisq2df << endl;
  }//GPed::statCompute


  /////////////////////////////////////////
  // The corresponding additive routines //
  /////////////////////////////////////////
  // Unified estimating equation (works for all three approaches!)
  // GPed::estEq_A
  void estEq_A( double *beta, int betaLength,
                MMatrix &ui ) {
    // ADDITIVE coding for xge
    // beta = bge (2), bg (2), be (1), bc (J) [[b0 is not modelled]]
    if( betaLength < 3 ) {
      // verify this precondition
      Rprintf("GPed::estEq_A(...) betaLength=%d, but it must be at least of length 3.\n", betaLength);
      //exit(1);
      return;
    }

    ui.resize( numFamilies(), betaLength );
    ui.fill( 0 ); // zero out ui

    // loop over families
    for( int i=0; i<numFamilies(); i++ ) {
      //double m[betaLength]; // i.e. the observed
      //double numSum[betaLength];
      vector<double> m; m.resize(betaLength);
      vector<double> numSum; numSum.resize(betaLength);
      double denSum = 0.0; // scalar
      for( int b=0; b<betaLength; b++ )
        m[b] = numSum[b] = 0.0;

      // 1) Compute the observed
      // loop over individuals
      for( int j=0; j<families[i].numChild(); j++ ) {
        if( families[i].childTrait[j] == 1 ) {
          // compute observed
          double x  = (double)(families[i].childGeno[j]); // additive coding
          double x0 = (double)(families[i].childGeno[j]==X0);  // genotype coding 0
          double x1 = (double)(families[i].childGeno[j]==X1);  // genotype coding 1
          double z  = families[i].childEnvironment[j];
          m[0] += x * z;
          m[1] += x0;
          m[2] += x1;
          if( betaLength > 3 )
            m[3] += z;
          if( betaLength > 4 )
            for( int c=4; c<betaLength; c++ )
              m[c] += families[i].childCovariate[j][c-4];
        }
      }

      // 2) Compute the expected
      // loop over geno perms
      for( unsigned int gp=0; gp<families[i].genoPerm.size(); gp++ ) {
        // loop over pheno perms
        for( unsigned int pp=0; pp<families[i].phenoPerm.size(); pp++ ) {
          // add the current to the expected sum
          //double mStar[betaLength];
          vector<double> mStar; mStar.resize(betaLength);
          for( int b=0; b<betaLength; b++ )
            mStar[b] = 0.0;
          // loop over individuals
          for( int j=0; j<families[i].numChild(); j++ ) {
            if( families[i].phenoPerm[pp][j] == 1 ) {
              double xstar  = (double)(families[i].genoPerm[gp][j]);
              double x0star = (double)(families[i].genoPerm[gp][j]==X0);
              double x1star = (double)(families[i].genoPerm[gp][j]==X1);
              double zstar  = families[i].childEnvironment[j]; // no perming here
              mStar[0] += xstar * zstar;
              mStar[1] += x0star;
              mStar[2] += x1star;
              if( betaLength > 3 )
                mStar[3] += zstar;
              if( betaLength > 4 )
                for( int c=4; c<betaLength; c++ )
                  mStar[c] += families[i].childCovariate[j][c-4];
            }//fi
          }//j
          double weight = families[i].genoPermWeight[gp] * families[i].phenoPermWeight[pp];
          double expStuff = 0.0;
          for( int b=0; b<betaLength; b++ )
            expStuff += beta[b] * mStar[b];
          expStuff = exp( expStuff ) * weight;

          for( int b=0; b<betaLength; b++ )
            numSum[b] += mStar[b] * expStuff;

          denSum += expStuff;
        }//pp
      }//gp

      for( int b=0; b<betaLength; b++ )
        ui.m[i][b] += m[b] - numSum[b] / denSum;
    }//i
  }//Gped::estEq_A - Crucial debug routine

  void statPrecompute_A( double *beta, int betaLength,
                         double *dbnuis_ret ) { // needs to be inverted in R, and then passed back here...
    // compute the derivatives
    M = betaLength - 1;
    //MMatrix dbge, dbnuis, u; // Need to be global variables
    dbge.resize( 1, M );
    dbnuis.resize( M, M );
    dbge.fill( 0.0 );
    dbnuis.fill( 0.0 );

    // go across the families
    for( int i=0; i<numFamilies(); i++ ) {
      double Ai00 = 0.0;
      MMatrix Ai11, Ai10, Ai01t, Ai02t;
      Ai11.resize( 1, M );
      Ai10.resize( 1, 1 ); // kind of wastefull...
      Ai01t.resize( M, 1 ); //( 1, M );
      Ai02t.resize( M, M ); // Transpose of Ai02, a little awkward notation...
      Ai11.fill(  0.0 );
      Ai10.fill(  0.0 );
      Ai01t.fill( 0.0 );
      Ai02t.fill( 0.0 );

      // 2) Compute the expected
      // loop over geno perms
      for( unsigned int gp=0; gp<families[i].genoPerm.size(); gp++ ) {
        // loop over pheno perms
        for( unsigned int pp=0; pp<families[i].phenoPerm.size(); pp++ ) {
          // add the current to the expected sum
          //double mStar[betaLength];
          vector<double> mStar; mStar.resize(betaLength);
          for( int b=0; b<betaLength; b++ )
            mStar[b] = 0.0;

          // loop over individuals
          for( int j=0; j<families[i].numChild(); j++ ) {
            if( families[i].phenoPerm[pp][j] == 1 ) {
              double xstar  = (double)(families[i].genoPerm[gp][j]);
              double x0star = (double)(families[i].genoPerm[gp][j]==X0);
              double x1star = (double)(families[i].genoPerm[gp][j]==X1);
              double zstar  = families[i].childEnvironment[j]; // no perming here
              mStar[0] += xstar * zstar;
              mStar[1] += x0star;
              mStar[2] += x1star;
              if( betaLength > 3 )
                mStar[3] += zstar;
              if( betaLength > 4 )
                for( int c=4; c<betaLength; c++ )
                  mStar[c] += families[i].childCovariate[j][c-4];
            }//fi
          }//j

          double weight = families[i].genoPermWeight[gp] * families[i].phenoPermWeight[pp];

          double expStuff = 0.0;
          for( int b=0; b<betaLength; b++ )
            expStuff += beta[b] * mStar[b];
          expStuff = exp( expStuff ) * weight;

          MMatrix xz;
          xz.createV( (double *)mStar.data(), 1 );
          MMatrix m;
          m.createV( &mStar[1], M );
          MMatrix m_mt;
          m_mt.createVVt( &mStar[1], M );

          Ai00 += expStuff;
          Ai11.addSelfC( xz.multiplyC( m.transpose() ).multiply( expStuff ) );
          Ai10.addSelfC( xz.multiply( expStuff ) );
          Ai01t.addSelfC( m.multiply( expStuff ) );
          Ai02t.addSelfC( m_mt.multiply( expStuff ) );
        }//pp
      }//gp


      // now add these pieces on to the derivative
      dbge.addSelfC( Ai11.multiply(Ai00).subtractC( Ai10.multiplyC(Ai01t.transpose()) ).multiply( -1.0/Ai00/Ai00 ) );
      dbnuis.addSelfC( Ai02t.multiply(Ai00).subtractC( Ai01t.multiplyC( Ai01t.transpose() ) ).multiply( -1/Ai00/Ai00 ) );
    }// i

    // get the ui matrix
    estEq_A( beta, betaLength, u );

    // copy dbnuis into dbnuis_ret
    for( int r=0; r<dbnuis.nrows(); r++ )
      for( int c=0; c<dbnuis.ncols(); c++ )
        dbnuis_ret[ r + c*M ] = dbnuis.m[r][c]; // should be MxM

    // Now return to R to get the inverse, which will then be passed to statCompute, and we will continue
  }

  void statCompute_A( double *dbnuis_inv, // M x M
                      double *chisq1df, int *numInf ) { // return value follows chi-squared distribution with 2 degrees of freedom
    // Put the inverse into a MMatrix
    MMatrix dbnuisInv;
    dbnuisInv.resize( M, M );
    for( int r=0; r<M; r++ )
      for( int c=0; c<M; c++ )
        dbnuisInv.m[r][c] = dbnuis_inv[ r + c*M ];
    //cout << "statCompute_A dbnuisInv" << endl << dbnuisInv.toString() << endl;

    MMatrix uge = u.subMatrix( 0, numFamilies()-1, 0, 0 ); // rowstart, rowEnd, colStart, colEnd
    MMatrix ug  = u.subMatrix( 0, numFamilies()-1, 1, M+1-1 );
    MMatrix w = uge.subtractC( dbge.multiplyC( dbnuisInv ).multiplyC( ug.transpose() ).transpose() );

    *numInf = 0;
    double sumwi=0.0, sumwiwit=0.0;
    for( int i=0; i<numFamilies(); i++ ) {
      double temp = w.m[i][0];
      sumwi += temp;
      sumwiwit += temp*temp;

      //////*numInf += (temp==0); // oh, you've got to be kidding me, it's backwards!!!
      //*numInf += ( (temp > ZEROTOL) | (temp < -ZEROTOL) );
      *numInf += (uge.m[i][0] < -ZEROTOL) || (uge.m[i][0] > ZEROTOL);

      // This next piece is just for debugging the dataset analysis, but nonetheless
      //cout << "Informative: " << (int)(temp==0) << endl;
      //cout << families[i].toString(true) << endl;
      // debug end
    }

    if( sumwiwit==0.0 ) {
      *chisq1df = 0.0;
    }else{
      *chisq1df = sumwi * sumwi / sumwiwit;
    }
  }//GPed::statCompute

  // Setup from an external dataset (i.e. the one we want to analyze)
  // GPed::set
  void set( int *pid, int *id, int *idfath, int *idmoth,
            int *geno, int *trait, double *env, int n ) {
    families.clear();

    int start = 0;
    vector<int> parents;
    vector<int> children;
    for( int i=0; i<n; i++ ) {
      ////cout << "set " << i << " " << idfath[i] << " " << idmoth[i] << " " << geno[i] << endl;
      if( pid[start]==pid[i] ) {
        if( idfath[i]==0 || idmoth[i]==0 ) {
          ////cout << "parental geno " << i << " " << geno[i] << endl;
          // then it's a parent, but only push on if non-missing
          if( geno[i] != GFamily::GMISS )
            parents.push_back(i);
        }else{
          children.push_back(i);
        }//fi
      }//fi

      if( pid[start]!=pid[i] || i==n-1 ) {
        // don't push on cases with no offspring, or cases with 1 offspring but no parents
        if( children.size() > 0 && !( children.size()==1 && parents.size()<2 ) ) {
          // push into gped structure
          GFamily newFam;

          newFam.parentGeno[0] = newFam.parentGeno[1] = GFamily::GMISS;
          for( unsigned int p=0; p<parents.size(); p++ )
            newFam.parentGeno[p] = geno[ parents[p] ];

          for( unsigned int c=0; c<children.size(); c++ ) {
            newFam.childGeno.push_back( geno[ children[c] ] );
            newFam.childTrait.push_back( trait[ children[c] ] );
            newFam.childEnvironment.push_back( env[ children[c] ] );
          }//c

          families.push_back( newFam );
        }//fi

        // and clear out the old, regardless of whether pushed on or not
        start = i;
        parents.clear();
        children.clear();

        // We didn't push on the individual who changed the pids, so we need to set our iterator _back_ one!
        if( i!=n-1 )
          i--;
      }//fi
    }//i
  }//GPed::set
  // Setup from an external dataset (i.e. the one we want to analyze)
  // GPed::set_C
  void set_C( int *pid, int *id, int *idfath, int *idmoth,
              int *geno, int *trait, double *env, double *cov, int nCov, int n ) {
    families.clear();

    int start = 0;
    vector<int> parents;
    vector<int> children;
    for( int i=0; i<n; i++ ) {
      if( pid[start]==pid[i] ) {
        if( idfath[i]==0 || idmoth[i]==0 ) {
          if( geno[i] != GFamily::GMISS )
            parents.push_back(i);
        }else{
          children.push_back(i);
        }//fi
      }//fi

      if( pid[start]!=pid[i] || i==n-1 ) {
        // don't push on cases with no offspring, or cases with 1 offspring but no parents
        if( children.size() > 0 && !( children.size()==1 && parents.size()<2 ) ) {
          // push into gped structure
          GFamily newFam;

          newFam.parentGeno[0] = newFam.parentGeno[1] = GFamily::GMISS;
          for( unsigned int p=0; p<parents.size(); p++ )
            newFam.parentGeno[p] = geno[ parents[p] ];

          for( unsigned int c=0; c<children.size(); c++ ) {
            newFam.childGeno.push_back( geno[ children[c] ] );
            newFam.childTrait.push_back( trait[ children[c] ] );
            newFam.childEnvironment.push_back( env[ children[c] ] );

            if( nCov > 0 ) {
              vector<double> cCov;
              for( int v=0; v<nCov; v++ )
                cCov.push_back( cov[ children[c] + n*v ] ); // index is row first, then mult by cols
              newFam.childCovariate.push_back( cCov );
            }
          }//c

          families.push_back( newFam );
        }//fi

        // and clear out the old, regardless of whether pushed on or not
        start = i;
        parents.clear();
        children.clear();

        // We didn't push on the individual who changed the pids, so we need to set our iterator _back_ one!
        if( i!=n-1 )
          i--;
      }//fi
    }//i
  }//GPed::set
};//GPed

// Is there really any need to have the reference, like we did last time?
// -- I don't think so, so eliminate it!
GPed gped; // GPed global object used to communicate with R

extern "C" {
  void cpp_gped_clear() { gped.clear(); }
  void cpp_gped_print( int *extended ) {
    Rprintf("%s\n", gped.toString( *extended != 0 ).c_str());
  }
  void cpp_gped_numCovariates( int *ret ) { *ret = gped.numCovariates(); }
  void cpp_gped_setStrategy( char **strat ) { gped.setStrategy( *strat ); }
  void cpp_gped_estEq( double *beta, int *betaLength, double *u ) {
    gped.estEq( beta, *betaLength, u );
  }
  void cpp_gped_statPrecompute( double *beta, int *betaLength,
                                double *dbnuis_ret ){
    gped.statPrecompute( beta, *betaLength, dbnuis_ret );
  }
  void cpp_gped_statCompute( double *dbnuis_inv, double *chisq2df, int *numInf ) {
    gped.statCompute( dbnuis_inv, chisq2df, numInf );
  }

  void cpp_gped_statPrecompute_A( double *beta, int *betaLength,
                                  double *dbnuis_ret ){
    gped.statPrecompute_A( beta, *betaLength, dbnuis_ret );
  }
  void cpp_gped_statCompute_A( double *dbnuis_inv, double *chisq1df, int *numInf ) {
    gped.statCompute_A( dbnuis_inv, chisq1df, numInf );
  }

  void cpp_gped_set( int *pid, int *id, int *idfath, int*idmoth,
                     int *geno, int *trait, double *env, int *n ) {
    gped.set( pid, id, idfath, idmoth,
              geno, trait, env, *n );
  }
  void cpp_gped_set_C( int *pid, int *id, int *idfath, int*idmoth,
                       int *geno, int *trait, double *env, double *cov, int *nCov, int *n ) {
    gped.set_C( pid, id, idfath, idmoth,
                geno, trait, env, cov, *nCov, *n );
  }
}

////////////////////
// GESimSub class //
////////////////////
class GESimSub {
public:
  enum LINK { LOG, LOGISTIC }; // see character strings in toString() that correspond
  enum ENV { DICHOTOMOUS, NORMAL }; // ""
  enum GENETIC { ADDITIVE, DOMINANT, RECESSIVE }; // ""

  int numParents, numOffspring, numFam;
  int minAffected, maxAffected;
  double afreq; GENETIC geneticModel;
  LINK link;
  vector<double> beta; // b0, bge, bg, be
  ENV env; double envCutoff; Random r;

  vector< vector<int> > perm;

  double IMPORTANCE_MAX;

  // for misspecified covariate
  double betaCov;
  ENV distCov;

  // for arbitrary phenotypic correlation
  double phenoCor;
  double phenoCutoff;
  Random r2;

  // for a DSL marker in LD with our marker
  double markerCor;
  double markerAfreq;

  // for arbitrary phenotypic correlation - better
  double phenoOR;

  //GESimSub::toString
  string toString() {
    const char* LINK_STR[] = {"log","logistic"};
    const char* ENV_STR[] = {"dichotomous","normal"};
    const char* GENETIC_STR[] = {"additive","dominant","recessive"};

    // creates a string of this object (mostly just for debug purposes)
    string str;
    str = "numParents="+d2s(numParents) + " numOffspring="+d2s(numOffspring)+" numFam="+d2s(numFam) + "\n";
    str += "minAffected="+d2s(minAffected)+ " maxAffected="+d2s(maxAffected) + "\n";
    str += "afreq="+d2s(afreq) + " genetic_model=" + d2s( (int)geneticModel ) + "(" + GENETIC_STR[geneticModel] + ")" + "\n";
    str += "link="+d2s(link) + "(" + LINK_STR[link] + ")" + "\n";
    str += "beta="; for( unsigned int b=0; b<beta.size(); b++ ) str+=d2s(beta[b])+",";
    str += "\n";
    str += "env=" + d2s(env) + "(" + ENV_STR[env] + ")" + " envCutoff=" + d2s(envCutoff) + "\n";

    str += "betaCov=" + d2s(betaCov) + " distCov=" + ENV_STR[distCov] + "\n";
    str += "phenoCor=" + d2s(phenoCor) + " phenoCutoff=" + d2s(phenoCutoff) + "\n";
    str += "markerCor=" + d2s(markerCor) + "\n";

    if( perm.size() > 0 ) {
      str += "perm = ";
      for( unsigned int p=0; p<perm.size(); p++ ) {
        for( unsigned int j=0; j<perm[p].size(); j++ )
          str += d2s(perm[p][j]);
        str += " ";
      }
      str += "\n";
    }else{
      str += "perm EMPTY\n";
    }

    str += "IMPORTANCE_MAX=" + d2s(IMPORTANCE_MAX) + "\n";

    return( str );
  }//GESimSub::toString -- SIMPLE

  // GESimSub::setPossiblePerms
  void setPossiblePerms() {
    // based on ascertainment strategy

    //perm.clear();
    //vector<int> curPerm;
    //for( int naff=minAffected; naff<=maxAffected; naff++ )
    //  perm2categories( perm, curPerm, 0, naff, 1, 0, numOffspring );

    // proband _must_ be affected, so we need to alter the above just a little...?
    perm.clear();

    // setup the curPerm so that the first offspring is affected
    vector<int> curPerm;
    curPerm.resize(numOffspring);
    curPerm[0] = 1; // proband must be affected
    for( unsigned int j=1; j<curPerm.size(); j++ )
      curPerm[j] = 0;

    for( int naff=minAffected-1; naff<=maxAffected-1; naff++ ) {
      if( naff==0 ) { // will only happen potentially on the first one, so ok
        // then need to just push curPerm on there
        perm.push_back( curPerm );
      }else{
        perm2categories( perm, curPerm, 1, naff, 1, 0, 0 );
      }
    }//naff
  }//GESimSub::setPossiblePerms -- not debugged?

  // GESimSub::xcode
  int xcode( int g ) {
    if( geneticModel==ADDITIVE ) {
      return( g ); // pretty much already set then! that's the new idea
    }else if( geneticModel==DOMINANT ) {
      //return( g==1 || g==2 );
      return( g != 0 );
    }else if( geneticModel==RECESSIVE ) {
      return( g == 2 );
    }
    Rprintf("GESimSub::xcode not ADDITIVE, DOMINANT, or RECESSIVE.\n");
    //exit(1);
    return( -999 );
  }//GESimSub::xcode -- SIMPLE

  // GESimSub::pd
  double pd( int g, double z ) {
    double x = (double)xcode(g);
    double mu = beta[0] + beta[1]*x*z + beta[2]*x + beta[3]*z;
    if( link==LOG ) {
      return( exp( mu ) );
    }else if( link==LOGISTIC ) {
      return( exp( mu ) / ( 1 + exp( mu ) ) );
    }
    Rprintf("GESimSub::pd link function incorrect(%d).\n", link);
    //exit( 1 );
    return( 0 );
  }//GESimSub::pd -- SIMPLE
  double pd_cov( int g, double z, double cov ) {
    double x = (double)xcode(g);
    double mu = beta[0] + beta[1]*x*z + beta[2]*x + beta[3]*z + betaCov*cov;
    if( link==LOG ) {
      return( exp( mu ) );
    }else if( link==LOGISTIC ) {
      return( exp( mu ) / ( 1 + exp( mu ) ) );
    }
    Rprintf("GESimSub::pd link function incorrect(%d).\n", link);
    //exit( 1 );
    return( 0 );
  }//GESimSub::pd -- SIMPLE

  //GESimSub::setImportanceSampling
  void setImportanceSampling() {
    IMPORTANCE_MAX = 1.0;
    //return; // DEBUG ONLY DEBUG ONLY DEBUG ONLY

    double pr[4]; pr[0] = 0; pr[1] = 0; pr[2] = 0; pr[3] = 0; // to values to avoid pedantic warnings, even though it would have to in the code below
    if( env==DICHOTOMOUS ) {
      // all possible combinations for the max/min (depends on sign of beta, which is why this is as complicated as it is)
      pr[0] = pd( 2, 1 );
      pr[1] = pd( 0, 1 );
      pr[2] = pd( 2, 0 );
      pr[3] = pd( 0, 0 );
    }else if( env==NORMAL ) {
      pr[0] = pd( 2, envCutoff );
      pr[1] = pd( 0, envCutoff );
      pr[2] = pd( 2, -envCutoff );
      pr[3] = pd( 0, -envCutoff );
    }else{
      Rprintf("GESimSub::setImportanceSampling, env type does not exist.\n");
    }

    // compute the max (for affected) / min (for unaffected) of an offspring
    double mmax = pr[0];
    double mmin = pr[0];

    for( int i=1; i<4; i++ ) {
      if( pr[i] < mmin )
        mmin = pr[i];

      if( pr[i] > mmax )
        mmax = pr[i];
    }

    // now compute the maximum possible under the ascertainment scheme
    double mmax_pd = 0;
    for( unsigned int p=0; p<perm.size(); p++ ) {
      double mmax_pdi = 1.0;
      for( unsigned int j=0; j<perm[p].size(); j++ ) {
        if( perm[p][j] == 1 ) {
          // affected
          mmax_pdi *= mmax;
        }else if( perm[p][j] == 0 ) {
          mmax_pdi *= (1-mmin);
        }
      }//j

      //if( mmax_pdi > mmax_pd )
      //  mmax_pd = mmax_pdi;
      mmax_pd += mmax_pdi;
    }//p

    IMPORTANCE_MAX = mmax_pd;
  }//GESimSub::setImportanceSampling --

  //GESimSub::set
  void set( int numParents, int numOffspring, int numFam,
            int minAffected, int maxAffected,
            double afreq, char* geneticModel,
            char* link,
            double* beta, int betaLength, // b0, bge, bg, be
            char* env, double envCutoff, double *chol, int p, // envCutoff is either the qtl number to dichotomize, or the maximum a normal can be (before it will be truncated)
            double betaCov, char* distCov, // NEW
            double *phenoCor, double phenoCutoff,
            double markerCor, double markerAfreq,
            double phenoOR ) {
    // save all of the arguments
    this->numParents = numParents;
    this->numOffspring = numOffspring;
    this->numFam = numFam;

    this->minAffected = minAffected;
    this->maxAffected = maxAffected;

    this->afreq = afreq;
    string geneticModelString = geneticModel;
    if( geneticModelString == "additive" ) {
      this->geneticModel = ADDITIVE;
    }else if( geneticModelString == "dominant" ) {
      this->geneticModel = DOMINANT;
    }else if( geneticModelString == "recessive" ) {
      this->geneticModel = RECESSIVE;
    }else{
      Rprintf("GESimSub genetic model must be 'additive', 'dominant', or 'recessive'. You supplied '%d'.\n", geneticModel);
      //exit(1);
      return;
    }

    string linkString = link;
    if( linkString == "log" ) {
      this->link = LOG;
    }else if( linkString=="logit" || linkString == "logistic" ) {
      this->link = LOGISTIC;
    }else{
      Rprintf("GESimSub::set link function must be 'log' or 'logit' ('logistic' also accepted); you supplied '%d'.\n", link);
      //exit(0);
      return;
    }

    this->beta.resize( betaLength );
    for( int b=0; b<betaLength; b++ )
      this->beta[b] = beta[b];

    string envString = env;
    if( envString == "dichotomous" ) {
      this->env = DICHOTOMOUS;
    }else if( envString == "normal" ) {
      this->env = NORMAL;
    }
    this->envCutoff = envCutoff;
    if( this->env==DICHOTOMOUS || this->env==NORMAL ) {
      this->r.setNormalSigma( chol, p );
    }

    // new things to set for model misspecification
    this->betaCov = betaCov;
    string distCovString = distCov;
    if( distCovString == "dichotomous" ) {
      this->distCov = DICHOTOMOUS;
    }else if( distCovString == "normal" ) {
      this->distCov = NORMAL;
    }else{
      Rprintf("GESimSub::set distCov must be either 'normal' or 'dichotomous', not '%s'.\n", distCovString.c_str());
      //exit(1);
      return;
    }

    this->phenoCutoff = phenoCutoff;
    if( phenoCutoff != 0 )
      this->r2.setNormalSigma( phenoCor, p );

    this->markerCor = markerCor;
    this->markerAfreq = markerAfreq;

    this->phenoOR = phenoOR;

    // create all possible phenotype according to the ascertainment
    setPossiblePerms();

    // setup the maximum for importance sampling!
    setImportanceSampling();
  }//GESimSub::set -- SIMPLE

  // GESimSub::draw
  void draw( GFamily &fam ) {
    int MAX_TRY = 10000;
    for( int t=0; t<MAX_TRY; t++ ) {
      fam.clear();

      // draw up the parents from allele frequencies
      fam.parentGeno[0] = (int)(r.runif()<afreq) + (int)(r.runif()<afreq);
      fam.parentGeno[1] = (int)(r.runif()<afreq) + (int)(r.runif()<afreq);

      // draw up the children using Mendel's laws
      for( int j=0; j<numOffspring; j++ ) {
        int childGeno = 0;
        for( int parent=0; parent<2; parent++ ) {
          if( fam.parentGeno[parent]==1 ) { // 50/50 passes it on
            childGeno += (int)(r.runif()<0.5);
          }else if( fam.parentGeno[parent]==2 ) { // has to pass on the allele
            childGeno++;
          } // otherwise parent doesn't have allele, can't pass it on
        }
        fam.childGeno.push_back( childGeno );
      }

      // draw up the correlated? environment (childEnvironment)
      if( env==DICHOTOMOUS ) {
        r.mvrnorm( fam.childEnvironment );
        for( int j=0; j<numOffspring; j++ )
          fam.childEnvironment[j] = (double)( fam.childEnvironment[j] < envCutoff );
      }else if( env==NORMAL ) {
        r.mvrnormTrunc( fam.childEnvironment, envCutoff );
        //r.mvrnorm( fam.childEnvironment );
        //for( int j=0; j<numOffspring; j++ )
        //  if( fam.childEnvironment[j] > envCutoff )
        //    fam.childEnvironment[j] = envCutoff;
      }

      // compute the probability of diseased of each offspring
      //double py[ numOffspring ];
      vector<double> py; py.resize(numOffspring);
      for( int j=0; j<numOffspring; j++ )
        py[j] = pd( fam.childGeno[j], fam.childEnvironment[j] );

      // for each possible phenotype combination, compute P( Y | ... )
      // sum these, the result is P( ascertainment criterion | ... )
      double sumpy = 0.0;
      //double prd[ perm.size() ];
      vector<double> prd; prd.resize(perm.size());
      for( unsigned int p=0; p<perm.size(); p++ ) {
        prd[p] = 1.0;
        for( unsigned int j=0; j<perm[p].size(); j++ ) {
          if( perm[p][j] == 1 ) {
            prd[p] *= py[j];
          }else{
            prd[p] *= (1-py[j]);
          }
        }
        sumpy += prd[p];
      }

      // divide by importance sampling max, and keep the result if a random uniform is less than this result, otherwise repeat
      if( sumpy/IMPORTANCE_MAX > 1.00001 ) { // allow for a little numerical error to perhaps propagate through here, I guess, although I don't really expect it to?
        Rprintf("Importance sampling mistake -- one of the probabilities goes above 1 (%f)\n", sumpy/IMPORTANCE_MAX);
        //exit(1);
        return;
      }
      if( r.runif() < sumpy/IMPORTANCE_MAX ) {
        // if found a good one then randomly choose the phenotype

        // normalize prd
        for( unsigned int p=0; p<perm.size(); p++ )
          prd[p] /= sumpy;

        // now draw one of them
        double curp = 0.0;
        double rnd = r.runif();
        for( unsigned int p=0; p<perm.size(); p++ ) {
          curp += prd[p];
          if( rnd <= curp )
            fam.childTrait = perm[p];
        }
        if( fam.childTrait.size() == 0 ) // wasn't set in the above, really shouldn't happen (measure zero)
          fam.childTrait = perm[perm.size()-1]; // shouldn't hit this often idealy...

        // lastly, do we need to erase the parents?
        if( numParents != 2 )
          fam.parentGeno[0] = fam.parentGeno[1] = GFamily::GMISS;
        return;
      }

      //cout << "Need to do another try!" << endl; // DEBUG only, to make sure this is really working, and not just accepting everyone
    }// t

    // past the number of tries
    Rprintf("GESimSub::draw() exceeded maximum possible tries to get a good family, terminating.\n");
    //exit(1);
  }//GESimSub::draw -- SEMI-DEBUGGED

  void inefficientDraw( GFamily &fam ) {
    int MAX_TRY = 100000;
    for( int t=0; t<MAX_TRY; t++ ) {
      fam.clear();

      // draw up the parents from allele frequencies
      fam.parentGeno[0] = (int)(r.runif()<afreq) + (int)(r.runif()<afreq);
      fam.parentGeno[1] = (int)(r.runif()<afreq) + (int)(r.runif()<afreq);

      // draw up the children using Mendel's laws
      for( int j=0; j<numOffspring; j++ ) {
        int childGeno = 0;
        for( int parent=0; parent<2; parent++ ) {
          if( fam.parentGeno[parent]==1 ) { // 50/50 passes it on
            childGeno += (int)(r.runif()<0.5);
          }else if( fam.parentGeno[parent]==2 ) { // has to pass on the allele
            childGeno++;
          } // otherwise parent doesn't have allele, can't pass it on
        }
        fam.childGeno.push_back( childGeno );
      }

      // draw up the correlated? environment (childEnvironment)
      if( env==DICHOTOMOUS ) {
        r.mvrnorm( fam.childEnvironment );
        for( int j=0; j<numOffspring; j++ )
          fam.childEnvironment[j] = (double)( fam.childEnvironment[j] < envCutoff );
      }else if( env==NORMAL ) {
        r.mvrnormTrunc( fam.childEnvironment, envCutoff );
        //r.mvrnorm( fam.childEnvironment );
        //for( int j=0; j<numOffspring; j++ )
        //  if( fam.childEnvironment[j] > envCutoff )
        //    fam.childEnvironment[j] = envCutoff;
      }

      // draw up Y
      fam.childTrait.resize(numOffspring);
      for( int j=0; j<numOffspring; j++ )
        fam.childTrait[j] = (int)( r.runif() < pd( fam.childGeno[j], fam.childEnvironment[j] ) );

      // does Y satisfy the ascertainment criterion?
      if( fam.childTrait[0] != 1 )
        continue; // first child must be affected
      int nAffected = 0;
      for( int j=0; j<numOffspring; j++ )
        nAffected += fam.childTrait[j];
      if( nAffected>=minAffected && nAffected<=maxAffected ) {
        // lastly, do we need to erase the parents?
        if( numParents != 2 )
          fam.parentGeno[0] = fam.parentGeno[1] = GFamily::GMISS;
        return; // SUCCESS!
      }
    }

    // hit the maximum number of tries (i.e. it's just going to be way to inefficient)
    Rprintf("GeSimSub::inefficientDraw() hit maximum number of tries.\n");
    //exit(0);
  }

  void inefficientDraw_missedCovariate( GFamily &fam ) {
    //cout << "MISSED COVARIATE " << "betaCov=" << betaCov << " distCov=" << (int)distCov << endl; exit(1);
    // mostly like the previous routine...

    int MAX_TRY = 100000;
    for( int t=0; t<MAX_TRY; t++ ) {
      fam.clear();

      // draw up the parents from allele frequencies
      fam.parentGeno[0] = (int)(r.runif()<afreq) + (int)(r.runif()<afreq);
      fam.parentGeno[1] = (int)(r.runif()<afreq) + (int)(r.runif()<afreq);

      // draw up the children using Mendel's laws
      for( int j=0; j<numOffspring; j++ ) {
        int childGeno = 0;
        for( int parent=0; parent<2; parent++ ) {
          if( fam.parentGeno[parent]==1 ) { // 50/50 passes it on
            childGeno += (int)(r.runif()<0.5);
          }else if( fam.parentGeno[parent]==2 ) { // has to pass on the allele
            childGeno++;
          } // otherwise parent doesn't have allele, can't pass it on
        }
        fam.childGeno.push_back( childGeno );
      }

      // draw up the correlated? environment (childEnvironment)
      if( env==DICHOTOMOUS ) {
        r.mvrnorm( fam.childEnvironment );
        for( int j=0; j<numOffspring; j++ )
          fam.childEnvironment[j] = (double)( fam.childEnvironment[j] < envCutoff );
      }else if( env==NORMAL ) {
        r.mvrnormTrunc( fam.childEnvironment, envCutoff );
      }

      // draw up the correlated? missing covariate [[NEW]]
      //double missingCovariate[numOffspring];
      vector<double> missingCovariate; missingCovariate.resize(numOffspring);
      if( distCov==DICHOTOMOUS ) {
        Rprintf("missing dichotomous covariate not yet supported\n");
        //exit(1);
        return;
      }else if( distCov==NORMAL ) {
        for( int j=0; j<numOffspring; j++ )
          missingCovariate[j] = r.rnormTrunc(1.64);
      }

      // draw up Y
      fam.childTrait.resize(numOffspring);
      for( int j=0; j<numOffspring; j++ )
        fam.childTrait[j] = (int)( r.runif() < pd_cov( fam.childGeno[j], fam.childEnvironment[j], missingCovariate[j] ) );

      // does Y satisfy the ascertainment criterion?
      if( fam.childTrait[0] != 1 )
        continue; // first child must be affected
      int nAffected = 0;
      for( int j=0; j<numOffspring; j++ )
        nAffected += fam.childTrait[j];
      if( nAffected>=minAffected && nAffected<=maxAffected ) {
        // lastly, do we need to erase the parents?
        if( numParents != 2 )
          fam.parentGeno[0] = fam.parentGeno[1] = GFamily::GMISS;
        return; // SUCCESS!
      }
    }

    // hit the maximum number of tries (i.e. it's just going to be way to inefficient)
    Rprintf("GeSimSub::inefficientDraw() hit maximum number of tries.\n");
    //exit(0);
  }
  void inefficientDraw_phenoCor( GFamily &fam ) {
    //cout << "PHENOTYPIC CORRELATION " << "phenoCor=" << phenoCor << "phenoCutoff=" << phenoCutoff << endl; exit(1);

    int MAX_TRY = 100000;
    for( int t=0; t<MAX_TRY; t++ ) {
      fam.clear();

      // draw up the parents from allele frequencies
      fam.parentGeno[0] = (int)(r.runif()<afreq) + (int)(r.runif()<afreq);
      fam.parentGeno[1] = (int)(r.runif()<afreq) + (int)(r.runif()<afreq);

      // draw up the children using Mendel's laws
      for( int j=0; j<numOffspring; j++ ) {
        int childGeno = 0;
        for( int parent=0; parent<2; parent++ ) {
          if( fam.parentGeno[parent]==1 ) { // 50/50 passes it on
            childGeno += (int)(r.runif()<0.5);
          }else if( fam.parentGeno[parent]==2 ) { // has to pass on the allele
            childGeno++;
          } // otherwise parent doesn't have allele, can't pass it on
        }
        fam.childGeno.push_back( childGeno );
      }

      // draw up the correlated? environment (childEnvironment)
      if( env==DICHOTOMOUS ) {
        r.mvrnorm( fam.childEnvironment );
        for( int j=0; j<numOffspring; j++ )
          fam.childEnvironment[j] = (double)( fam.childEnvironment[j] < envCutoff );
      }else if( env==NORMAL ) {
        r.mvrnormTrunc( fam.childEnvironment, envCutoff );
      }

      // draw up Y
      fam.childTrait.resize(numOffspring);
      vector<double> traitTemp;
      r2.mvrnorm( traitTemp );
      for( int j=0; j<numOffspring; j++ ) {
        double x = fam.childGeno[j];
        double z = fam.childEnvironment[j];
        //fam.childTrait[j] = traitTemp[j] + beta[1]*x*z + beta[2]*x + beta[3]*z;
        //fam.childTrait[j] = (double)(fam.childTrait[j] >= phenoCutoff);
        traitTemp[j] += beta[1]*x*z + beta[2]*x + beta[3]*z;
        fam.childTrait[j] = (int)(traitTemp[j] >= phenoCutoff);
      }

      // does Y satisfy the ascertainment criterion?
      if( fam.childTrait[0] != 1 )
        continue; // first child must be affected
      int nAffected = 0;
      for( int j=0; j<numOffspring; j++ )
        nAffected += fam.childTrait[j];
      if( nAffected>=minAffected && nAffected<=maxAffected ) {
        // lastly, do we need to erase the parents?
        if( numParents != 2 )
          fam.parentGeno[0] = fam.parentGeno[1] = GFamily::GMISS;
        return; // SUCCESS!
      }
    }

    // hit the maximum number of tries (i.e. it's just going to be way to inefficient)
    Rprintf("GeSimSub::inefficientDraw_phenoCor() hit maximum number of tries.\n");
    //exit(0);

  }
  void inefficientDraw_markerCor( GFamily &fam ) {
    //cout << "MARKER CORRELATION markerCor=" << markerCor << " markerAfreq=" << markerAfreq << endl; exit(1);

    int MAX_TRY = 100000;

    double afreqMarker = markerAfreq; // ah, well, then...

    // preliminary calculations
    double pmc[2][2];
    double d = sqrt( markerCor * afreq * (1-afreq) * afreqMarker * (1-afreqMarker) );
    pmc[0][0] = afreq * afreqMarker + d;
    pmc[0][1] = afreq * ( 1 - afreqMarker ) - d;
    pmc[1][0] = ( 1 - afreq ) * afreqMarker - d;
    pmc[1][1] = ( 1 - afreq ) * (1 - afreqMarker ) + d;

    for( int a=0; a<2; a++ ) {
      for( int b=0; b<2; b++ ) {
        if( pmc[a][b] < 0 || pmc[a][b] > 1 ) {
          Rprintf("Marker correlation too high, cannot be a distribution. Reduce markerCorrelation.\n");
          //exit(1);
          return;
        }
      }
    }

    for( int t=0; t<MAX_TRY; t++ ) {
      // draw up the parents
      int parentsM[4], parentsC[4];
      for( int i=0; i<4; i++ ) {
        double u = r.runif();
        if( u < pmc[0][0] ) {
          parentsM[i] = 0;
          parentsC[i] = 0;
        }else if( u < pmc[0][0] + pmc[0][1] ) {
          parentsM[i] = 0;
          parentsC[i] = 1;
        }else if( u < pmc[0][0] + pmc[0][1] + pmc[1][0] ) {
          parentsM[i] = 1;
          parentsC[i] = 0;
        }else{
          parentsM[i] = 1;
          parentsC[i] = 1;
        }

        fam.parentGeno[0] = parentsM[0] + parentsM[1];
        fam.parentGeno[1] = parentsM[2] + parentsM[3];
      }

      // draw up the children using Mendel's laws
      //int childM1[numOffspring], childC1[numOffspring];
      //int childM2[numOffspring], childC2[numOffspring];
      vector<double> childM1, childC1, childM2, childC2;
      childM1.resize(numOffspring);
      childC1.resize(numOffspring);
      childM2.resize(numOffspring);
      childC2.resize(numOffspring);
      fam.childGeno.resize(numOffspring);
      for( int j=0; j<numOffspring; j++ ) {
        if( r.runif() < 0.5 ) {
          childM1[j] = parentsM[0];
          childC1[j] = parentsC[0];
        }else{
          childM1[j] = parentsM[1];
          childC1[j] = parentsC[1];
        }

        if( r.runif() < 0.5 ) {
          childM2[j] = parentsM[2];
          childC2[j] = parentsC[2];
        }else{
          childM2[j] = parentsM[3];
          childC2[j] = parentsC[3];
        }

        fam.childGeno[j] = childM1[j] + childM2[j];
      }

      // draw up the correlated? environment
      if( env==DICHOTOMOUS ) {
        r.mvrnorm( fam.childEnvironment );
        for( int j=0; j<numOffspring; j++ )
          fam.childEnvironment[j] = (double)( fam.childEnvironment[j] < envCutoff );
      }else if( env==NORMAL ) {
        r.mvrnormTrunc( fam.childEnvironment, envCutoff );
      }

      // draw up Y
      fam.childTrait.resize(numOffspring);
      for( int j=0; j<numOffspring; j++ )
        fam.childTrait[j] = (int)( r.runif() < pd( childC1[j]+childC2[j], fam.childEnvironment[j] ) );

      // does Y satisfy the ascertainment criterion?
      if( fam.childTrait[0] != 1 )
        continue; // first child must be affected
      int nAffected = 0;
      for( int j=0; j<numOffspring; j++ )
        nAffected += fam.childTrait[j];
      if( nAffected>=minAffected && nAffected<=maxAffected ) {
        // IT'S A GOOD ONE!
        // lastly, do we need to erase the parents?
        if( numParents != 2 )
          fam.parentGeno[0] = fam.parentGeno[1] = GFamily::GMISS;
        return; // SUCCESS!
      }
    }

    // hit the maximum number of tries (i.e. it's just going to be way to inefficient)
    Rprintf("GeSimSub::inefficientDraw() hit maximum number of tries.\n");
    //exit(0);
  }//

  double solveORc( double w, double y, double p ) {
    // a b | y
    // c d | z
    // -------
    // w x |
    //
    // OR = p = ad/bc

    double x = 1-w;
    double z = 1-y;

    double A = p - 1;
    double B = z - w - x*p - z*p;
    double C = p*x*z;

    double d = (-B - sqrt(B*B-4*A*C)) / 2.0 / A;
    //double b = x - d;
    double c = z - d;
    //double a = 1 - b - c - d;

    return( c );
  }//originally debugged in R Code

  void inefficientDraw_phenoOR( GFamily &fam ) {
    //cout << "phenoOR = " << phenoOR << endl;
    //exit(1);

    int MAX_TRY = 100000;
    if( numOffspring!=2 || minAffected!=1 || maxAffected!=1 ) {
      Rprintf("Currently phenoOR only works for numOffspring==2, minAffected==1, maxAffected==1.\n");
      //exit(1);
      return;
    }

    for( int t=0; t<MAX_TRY; t++ ) {
      fam.clear();

      // draw up the parents from allele frequencies
      fam.parentGeno[0] = (int)(r.runif()<afreq) + (int)(r.runif()<afreq);
      fam.parentGeno[1] = (int)(r.runif()<afreq) + (int)(r.runif()<afreq);

      // draw up the children using Mendel's laws
      for( int j=0; j<numOffspring; j++ ) {
        int childGeno = 0;
        for( int parent=0; parent<2; parent++ ) {
          if( fam.parentGeno[parent]==1 ) { // 50/50 passes it on
            childGeno += (int)(r.runif()<0.5);
          }else if( fam.parentGeno[parent]==2 ) { // has to pass on the allele
            childGeno++;
          } // otherwise parent doesn't have allele, can't pass it on
        }
        fam.childGeno.push_back( childGeno );
      }

      // draw up the correlated? environment (childEnvironment)
      if( env==DICHOTOMOUS ) {
        r.mvrnorm( fam.childEnvironment );
        for( int j=0; j<numOffspring; j++ )
          fam.childEnvironment[j] = (double)( fam.childEnvironment[j] < envCutoff );
      }else if( env==NORMAL ) {
        r.mvrnormTrunc( fam.childEnvironment, envCutoff );
      }

      double py0 = pd( fam.childGeno[0], fam.childEnvironment[0] );
      double py1 = pd( fam.childGeno[1], fam.childEnvironment[1] );
      double p = solveORc( py0, py1, phenoOR );
      if( r.runif() < p ) {
        fam.childTrait.resize(2);
        fam.childTrait[0] = 1;
        fam.childTrait[1] = 0;

        // lastly, do we need to erase the parents?
        if( numParents != 2 )
          fam.parentGeno[0] = fam.parentGeno[1] = GFamily::GMISS;
        return; // SUCCESS!
      }
    }

    // hit the maximum number of tries (i.e. it's just going to be way to inefficient)
    Rprintf("GeSimSub::inefficientDraw() hit maximum number of tries.\n");
    //exit(0);
  }


  // GeSimSub::draw
  void draw() {
    // essentially loops the above draw...
    int start = gped.numFamilies();
    int end = start + numFam-1;
    gped.resize(end+1);
    //for( int f=start; f<=end; f++ ) // '<=' intended here
    //  inefficientDraw( gped.families[f] );
    //  //draw( gped.families[f] );

    if( betaCov != 0 ) {
      for( int f=start; f<=end; f++ ) // '<=' intended here
        inefficientDraw_missedCovariate( gped.families[f] );
    }else if( phenoCutoff != 0 ) {
      for( int f=start; f<=end; f++ )
        inefficientDraw_phenoCor( gped.families[f] );
    }else if( markerCor != 0 ) {
      for( int f=start; f<=end; f++ )
        inefficientDraw_markerCor( gped.families[f] );
    }else if( phenoOR != 1 ) {
      for( int f=start; f<=end; f++ )
        inefficientDraw_phenoOR( gped.families[f] );
    }else{
      for( int f=start; f<=end; f++ )
        inefficientDraw( gped.families[f] );
    }
  }//GeSimSub::draw -- SIMPLE
};

/////////////////
// GESim class //
/////////////////
class GESim {
public:
  vector<GESimSub> simSub;

  // Completely clears it out
  void clear() {
    simSub.clear();
  }//GEsim::clear -- SIMPLE

  string toString() {
    string str;
    for( unsigned int s=0; s<simSub.size(); s++ )
      str += simSub[s].toString() + "\n";
    return( str );
  }//GESim::toString -- SIMPLE

  //
  // First set, before things are simulated
  // copied from GESimSub
  void set( int numParents, int numOffspring, int numFam,
            int minAffected, int maxAffected,
            double afreq, char* geneticModel,
            char* link,
            double* beta, int betaLength, // b0, bge, bg, be
            char* env, double envCutoff, double *chol, int p,
            double betaCov, char* distCov, // NEW
            double *phenoCor, double phenoCutoff,
            double markerCor, double markerAfreq,
            double phenoOR ) {
    GESimSub simSubNew;
    simSubNew.set( numParents, numOffspring, numFam,
                   minAffected, maxAffected,
                   afreq, geneticModel,
                   link,
                   beta, betaLength,
                   env, envCutoff, chol, p,
                   betaCov, distCov,
                   phenoCor, phenoCutoff,
                   markerCor, markerAfreq,
                   phenoOR );
    simSub.push_back( simSubNew );
  }//GESim::set -- SIMPLE

  // for simulations
  void draw() {
    // empty out the gped
    gped.clear();

    // push on each type of simulation (i.e. could be a mixture of 50/50 trios and sibpairs...)
    for( unsigned int s=0; s<simSub.size(); s++ )
      simSub[s].draw();
  }//GESim::draw -- SIMPLE

  // GeSim::setOnlyFirstAffected
  void setOnlyFirstAffected() {
    for( int i=0; i<gped.numFamilies(); i++ )
      for( int j=1; j<gped.families[i].numChild(); j++ )
        gped.families[i].childTrait[j] = 0; // unaffected
  }
};

GESim gesim; // GESim global object - used to communicate with R
extern "C" {
  void cpp_gesim_print() {
    Rprintf("%s\n", gesim.toString().c_str());
  }
  // set up for data generation
  // multiple calls allow for multiple different family designs
  //  in each simulation
  void cpp_gesim_set( int *numParents, int *numOffspring, int *numFam,
                      int *minAffected, int *maxAffected,
                      double *afreq, char** geneticModel,
                      char** link,
                      double* beta, int *betaLength, // b0, bge, bg, be
                      char** env, double *envCutoff, double *chol, int *p,
                      double* betaCov, char** distCov, // NEW
                      double* phenoCor, double* phenoCutoff,
                      double* markerCor, double *markerAfreq,
                      double* phenoOR ) {
    gesim.set( *numParents, *numOffspring, *numFam,
               *minAffected, *maxAffected,
               *afreq, *geneticModel,
               *link,
               beta, *betaLength,
               *env, *envCutoff, chol, *p,
               *betaCov, *distCov,
               phenoCor, *phenoCutoff,
               *markerCor, *markerAfreq,
               *phenoOR );

  }

  // draw up a new dataset (automatically clears old dataset)
  void cpp_gesim_draw() {
    rn.attach(); // do we want this here?
    gesim.draw();
    rn.detach();
  }

  void cpp_gesim_setOnlyFirstAffected() {
    gesim.setOnlyFirstAffected();
  }

  void cpp_gesim_clear() {
    gesim.clear();
  }
}



// Old debugging routines
/*
int main() {
  vector< vector<int> > perm;
  vector<int> permAddi;
  permAddi.push_back( 1 );
  permAddi.push_back( 2 );

  fanperms( permAddi, perm );
  permAddi.push_back( 3 );
  fanperms( permAddi, perm );
  permAddi.push_back( 4 );
  fanperms( permAddi, perm );
  printperms( perm );

  cout << "Second perm test." << endl;
  perm.clear();
  perms( permAddi, perm );
  printperms( perm );

  cout << "Third perm test." << endl;
  perm.clear();
  permAddi.clear();
  vector<double> weight, weightAddi;
  // resulting weight should be the columns of the permutation multiplied together
  // -- hack, permAddi must be an integer
  permAddi.push_back( 10 ); weightAddi.push_back( 0.1 );
  permAddi.push_back( 40 ); weightAddi.push_back( 0.4 );
  permAddi.push_back( 50 ); weightAddi.push_back( 0.5 ); // sum to 1
  fanpermsw( permAddi, weightAddi, perm, weight );
  fanpermsw( permAddi, weightAddi, perm, weight );
  printpermsw( perm, weight );


  cout << "Fourth permutation test (a)." << endl;
  perm.clear();
  permAddi.clear();
  perm2categories( perm, permAddi, 0, 2, 9, 0, 5 );
  printperms( perm );
  cout << "Fourth permutation test (b)." << endl;
  perm.clear();
  permAddi.clear();
  perm2categories( perm, permAddi, 0, 3, 9, 0, 7 );
  printperms( perm );


  cout << "GFamily tests, SS and affection perms..." << endl;
  GFamily f;
  f.debugPermAddFamily( 0,2, 2,3,4, 1,1,2, false );
  cout << f.toString() << endl;
  f.debugPermAddFamily( 0,0, 1,0,0, 1,0,0 ); // trio

  cout << "##############################" << endl;
  cout << "GFamily tests, PARENTS PRESENT" << endl;
  f.debugPermAddFamily( 0,0, 1,0,0, 1,0,0 );
  f.debugPermAddFamily( 0,1, 1,1,0, 1,0,0 );
  f.debugPermAddFamily( 0,2, 0,1,0, 1,0,0 );
  f.debugPermAddFamily( 1,1, 1,1,1, 1,0,0 );
  f.debugPermAddFamily( 1,2, 0,1,1, 0,0,1 );
  f.debugPermAddFamily( 2,2, 0,0,1, 0,0,1 );

  int GMISS = GFamily::GMISS;

  cout << "#######################" << endl;
  cout << "GFamily tests, sibpairs" << endl;
  f.debugPermAddFamily( GMISS,GMISS, 2,0,0, 1,0,0 );
  f.debugPermAddFamily( GMISS,GMISS, 1,1,0, 1,1,0 );
  f.debugPermAddFamily( GMISS,GMISS, 1,0,1, 1,0,0 );
  f.debugPermAddFamily( GMISS,GMISS, 0,2,0, 0,1,0 );
  f.debugPermAddFamily( GMISS,GMISS, 0,1,1, 0,0,1 );
  f.debugPermAddFamily( GMISS,GMISS, 0,0,2, 0,0,2 );

  cout << "#################" << endl;
  cout << "GFamily, sibtrios" << endl;
  f.debugPermAddFamily( GMISS,GMISS, 3,0,0, 2,0,0 );
  f.debugPermAddFamily( GMISS,GMISS, 2,1,0, 1,1,0 );
  f.debugPermAddFamily( GMISS,GMISS, 2,0,1, 1,0,0 );
  f.debugPermAddFamily( GMISS,GMISS, 1,2,0, 1,0,0 );
  f.debugPermAddFamily( GMISS,GMISS, 1,1,1, 1,0,0 );
  f.debugPermAddFamily( GMISS,GMISS, 1,0,2, 0,0,2 );
  f.debugPermAddFamily( GMISS,GMISS, 0,3,0, 0,2,0 );
  f.debugPermAddFamily( GMISS,GMISS, 0,2,1, 0,1,0 );
  f.debugPermAddFamily( GMISS,GMISS, 0,1,2, 0,0,0 );
  f.debugPermAddFamily( GMISS,GMISS, 0,0,3, 0,0,3 );

  cout << "#################" << endl;
  cout << "GFamily, sibquads" << endl;
  f.debugPermAddFamily( GMISS,GMISS, 4,0,0, 4,0,0 );
  f.debugPermAddFamily( GMISS,GMISS, 3,1,0, 3,0,0 );
  f.debugPermAddFamily( GMISS,GMISS, 3,0,1, 2,0,0 );
  f.debugPermAddFamily( GMISS,GMISS, 2,2,0, 1,0,0 );
  f.debugPermAddFamily( GMISS,GMISS, 2,1,1, 1,1,0 );
  f.debugPermAddFamily( GMISS,GMISS, 2,0,2, 1,0,2 );
  f.debugPermAddFamily( GMISS,GMISS, 1,3,0, 1,0,0 );
  f.debugPermAddFamily( GMISS,GMISS, 1,2,1, 1,0,0 );
  f.debugPermAddFamily( GMISS,GMISS, 1,1,2, 1,0,0 );
  f.debugPermAddFamily( GMISS,GMISS, 1,0,3, 1,0,0 );
  f.debugPermAddFamily( GMISS,GMISS, 0,4,0, 0,3,0 );
  f.debugPermAddFamily( GMISS,GMISS, 0,3,1, 0,2,1 );
  f.debugPermAddFamily( GMISS,GMISS, 0,2,2, 0,0,0 );
  f.debugPermAddFamily( GMISS,GMISS, 0,1,3, 0,0,2 );
  f.debugPermAddFamily( GMISS,GMISS, 0,0,4, 0,0,2 );

  // debugged up to here -- works flawlessly
}
*/
