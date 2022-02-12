/* Thomas Hoffmann
 * fbatDist.cpp
 * Created: 5/19/2007
 * Last modified: 2/11/2008
 */

/* This has been _completely_ verified via routines in 'powerSim.cpp'
 *  that compared this directly with the FBAT routines,
 *  at least for bi-allelic markers, which is the current target!
 * Hence this option is no longer necessary, and is removed,
 *  but could be put in if you suspect any abberant behaviour.
 */
//#define DEBUG_FBATDIST

#include <R.h>

#include <vector>
using namespace std;

#include "fbatdist.h"

bool first( int a, int b, int c ) {
  return( a>0 && b==0 && c==0 );
}
bool second( int b, int a, int c ) {
  return( a>0 && b==0 && c==0 );
}
bool third( int c, int b, int a ) {
  return( a>0 && b==0 && c==0 );
}

bool first_second( int a, int b, int c ) {
  return( a>0 && b>0 && c==0 );
}
bool first_third( int a, int c, int b ) {
  return( a>0 && b>0 && c==0 );
}
bool second_third( int c, int b, int a ) {
  return( a>0 && b>0 && c==0 );
}
bool all( int a, int b, int c ) {
  return( a>0 && b>0 && c>0 );
}

void printFamily( int *p1, int *p2,
                  int *ca, int *cb,
                  int numSibs ) {
  //cout << "P: " << p1[0] << " "  << p1[1] << ", " << p2[0] << " " << p2[1] << endl << "C: ";
  //for( int i=0; i<numSibs; i++ )
  //  cout << ca[i] << " " << cb[i] << ", ";
  //cout << endl;
  Rprintf("P: %d %d, %d %d\nC: ", p1[0], p1[1], p2[0], p2[1]);
  for(int i=0; i<numSibs; i++)
    Rprintf("%d %d, ", ca[i], cb[i]);
  Rprintf("\n");
}

bool pG( int gP1, int gP2, // parental mating type
         int nAA, int nAB, int nBB, // number of offspring with said genotype
         double pg[3] ) // output parameter, prob of each genotype
{
  // swap gP1 and gP2 if gP1 is missing
  if( gP1==gMiss ) {
    gP1 = gP2;
    gP2 = gMiss;
  }

  double nSum = nAA + nAB + nBB;

  // neither parent is missing
  if( gP1!=gMiss ) {
    switch( gP1 ){
    case gAA:
      switch( gP2 ) {
      case gAA:
	pg[gAA] = 1.00;
	pg[gAB] = 0.00;
	pg[gBB] = 0.00;
        return(true);
      case gAB:
	pg[gAA] = 0.50;
	pg[gAB] = 0.50;
	pg[gBB] = 0.00;
        return(true);
      case gBB:
	pg[gAA] = 0.00;
	pg[gAB] = 1.00;
	pg[gBB] = 0.00;
        return(true);
      }

    case gAB:
      switch( gP2 ) {
      case gAA:
	pg[gAA] = 0.50;
	pg[gAB] = 0.50;
	pg[gBB] = 0.00;
        return(true);
      case gAB:
	pg[gAA] = 0.25;
	pg[gAB] = 0.50;
	pg[gBB] = 0.25;
        return(true);
      case gBB:
	pg[gAA] = 0.00;
	pg[gAB] = 0.50;
	pg[gBB] = 0.50;
        return(true);
      }

    case gBB:
      switch( gP2 ) {
      case gAA:
	pg[gAA] = 0.00;
	pg[gAB] = 1.00;
	pg[gBB] = 0.00;
        return(true);
      case gAB:
	pg[gAA] = 0.00;
	pg[gAB] = 0.50;
	pg[gBB] = 0.50;
        return(true);
      case gBB:
	pg[gAA] = 0.00;
	pg[gAB] = 0.00;
	pg[gBB] = 1.00;
        return(true);
      }
    }
  }

  // one parent is missing

  switch( gP1 ) {
  case gAA:
  case gBB:
    // one homozygous parent
    if( first(nAA,nAB,nBB) ) {
      pg[gAA] = 1.00;
      pg[gAB] = 0.00;
      pg[gBB] = 0.00;
      return(true);
    }
    if( second(nAA,nAB,nBB) ) {
      pg[gAA] = 0.00;
      pg[gAB] = 1.00;
      pg[gBB] = 0.00;
      return(true);
    }
    if( third(nAA,nAB,nBB) ) {
      pg[gAA] = 0.00;
      pg[gAB] = 0.00;
      pg[gBB] = 1.00;
      return(true);
    }

    if( first_second(nAA,nAB,nBB) ) {
      pg[gAA] = 0.50;
      pg[gAB] = 0.50;
      pg[gBB] = 0.00;
      return(true);
    }
    if( first_third(nAA,nAB,nBB) ) {
      pg[gAA] = 0.50;
      pg[gAB] = 0.00;
      pg[gBB] = 0.50;
      return(true);
    }
    if( second_third(nAA,nAB,nBB) ) {
      pg[gAA] = 0.00;
      pg[gAB] = 0.50;
      pg[gBB] = 0.50;
      return(true);
    }
  }

  // otherwise both parents are missing,
  //  or the parent is heterozygous
  // either are the same
  if( gP1==gAB || gP1==gMiss ) {
    // one heterozygous parent
    if( first(nAA,nAB,nBB) ) {
      pg[gAA] = 1.00;
      pg[gAB] = 0.00;
      pg[gBB] = 0.00;
      return(true);
    }
    if( second(nAA,nAB,nBB) ) {
      pg[gAA] = 0.00;
      pg[gAB] = 1.00;
      pg[gBB] = 0.00;
      return(true);
    }
    if( third(nAA,nAB,nBB) ) {
      pg[gAA] = 0.00;
      pg[gAB] = 0.00;
      pg[gBB] = 1.00;
      return(true);
    }

    if( !first_third(nAA,nAB,nBB) && !all(nAA,nAB,nBB) ){
      pg[gAA] = (double)nAA / nSum;
      pg[gAB] = (double)nAB / nSum;
      pg[gBB] = (double)nBB / nSum;
      return(true);
    }

    pg[gAA] = ( pow(4,nSum-1) - pow(3,nSum-1) ) / ( pow(4,nSum) - 2*pow(3,nSum) + pow(2,nSum) );
    pg[gBB] = pg[gAA];
    pg[gAB] = 1.0 - pg[gAA] - pg[gBB];
    return(true);
  }

  if( (gP1==gAA||gP1==gBB) && gP2==gMiss
      && all(nAA,nAB,nBB) ) {
#ifndef LOOKUP_COMPARE
    Rprintf("WARNING: impossible genotype in file.\n");
    //printFamily();
#endif
    return(true);
  }

  // This really only happens when genotypes of all the children are missing
  //  in which case they will be tossed anyway,
  //  thus there is no reason to warn anymore.
  //cout << "failed all cases!" << endl;
  //cout << gP1 << " " << gP2 << ";" << nAA << " " << nAB << " " << nBB << endl;
  return(false);
}

// rewriting code
/*
bool pGG( int gP1, int gP2,
          int nAA, int nAB, int nBB,
          double pg[3], double pgg[3][3] ) {
  // swap gP1 and gP2 if gP1 is missing
  if( gP1==gMiss ) {
    gP1 = gP2;
    gP2 = gMiss;
  }

  // zero everything
  pg[0] = pg[1] = pg[2] = 0.0;
  pgg[0][0] = pgg[0][1] = pgg[0][2] = pgg[1][0] = pgg[1][1] = pgg[1][2] = pgg[2][0] = pgg[2][1] = pgg[2][2] = 0.0;

  // fill in the univariate
  if( !pG( gP1, gP2,  nAA, nAB, nBB,  pg ) )
    return( false );

  int n = nAA + nAB + nBB;

  if( gP1!=gMiss ) {
    // neither parent is missing
    for( int g1=gAA; g1<=gBB; g1++ )
      for( int g2=gAA; g2<=gBB; g2++ )
        pgg[g1][g2] = pg[g1] * pg[g2];
    return( true );
  }else{
    // one or both parents are missing (doesn't matter)
    if( first(nAA,nAB,nBB) ) {
      pgg[gAA][gAA] = 1.0;
      return( true );
    }
    if( second(nAA,nAB,nBB) ) {
      pgg[gAB][gAB] = 1.0;
      return( true );
    }
    if( third(nAA,nAB,nBB) ) {
      pgg[gBB][gBB] = 1.0;
      return( true );
    }
    if( first_second(nAA,nAB,nBB) ) {
      pgg[gAA][gAA] = nAA*(nAA-1) / (n*(n-1));
      pgg[gAB][gAB] = nAB*(nAB-1) / (n*(n-1));
      pgg[gAA][gAB] = pgg[gAB][gAA] = nAA*nAB/(n*(n-1)); // ? div by 2 ??
      return( true );
    }
    if( second_third(nAA,nAB,nBB) ) {
      pgg[gBB][gBB] = nBB*(nBB-1) / (n*(n-1));
      pgg[gAB][gAB] = nAB*(nAB-1) / (n*(n-1));
      pgg[gBB][gAB] = pgg[gAB][gBB] = nBB*nAB/(n*(n-1)); // ? div by 2 ??
      return( true );
    }
    if( first_third(nAA,nAB,nBB) || all(nAA,nAB,nBB) ) {
      pgg[gAA][gAA] = pgg[gBB][gBB] = ( pow(4,n-1) - pow(3,n-1) ) / ( pow(4,n)-2*pow(3,n)+pow(2,n) );
      pgg[gAA][gAB] = pgg[gAB][gAA] = pgg[gBB][gAB] = pgg[gAB][gBB] = 2 * pgg[gAA][gAA];
      pgg[gAA][gBB] = pgg[gBB][gAA] = pow(4,n-2) / ( pow(4,n) - 2*pow(3,n) + pow(2,n) );
      pgg[gAB][gAB] = ( pow(4,n-1) - 8*pow(3,n-2) + pow(2,n) ) / ( pow(4,n) - 2*pow(3,n) + pow(2,n) );
    }
  }

  return( true );
}
*/


// allele code
int gCode( int a, int b )
{
  if( a==0 || b==0 ) return( gMiss );
  if( a==ALLELE_A && b==ALLELE_A ) return( gAA ); // alteration used to be 1 and two
  if( a==ALLELE_B && b==ALLELE_B ) return( gBB );
  return( gAB );
}

// NOTE: assumes ALLELE_A is the disease allele
int xCode( int a, int b, int MODEL ) {
  switch( MODEL ){
  case( MODEL_ADDITIVE ):
    return( (int)(a==ALLELE_A) + (int)(b==ALLELE_A) );
  case( MODEL_DOMINANT ):
    return( (int)(a==ALLELE_A || b==ALLELE_A) );
  case( MODEL_RECESSIVE ):
    return( (int)(a==ALLELE_A && b==ALLELE_A) );
  }
  Rprintf("xCode (1) out of bounds! %d %d\n", a, b);
  return( -1 ); // should never get here
}
int xCode( int g, int MODEL ) {
  switch( g ) {
  case gAA:
    return( xCode( ALLELE_A, ALLELE_A, MODEL ) );
  case gAB:
    return( xCode( ALLELE_A, ALLELE_B, MODEL ) );
  case gBB:
    return( xCode( ALLELE_B, ALLELE_B, MODEL ) );
  }
  Rprintf("xCode (2) out of bounds! %d\n", g);
  return( -1 ); // should never get here
}

int ggConvert( int g1, int g2 ) {
  return( g1*3 + g2 );
}

bool pG( int n,
         int *p1, int *p2, // parental alleles
         int *ca, int *cb, // childrens alleles
         double pg[3] ) // output parameter, prob of each genotype
{
  int nG[3] = {0,0,0};
  for( int i=0; i<n; i++ ) {
    int gcodefme = gCode( ca[i], cb[i] );
    if( gcodefme!=-1 )
      nG[ gcodefme ]++;
  }

  return( pG( gCode(p1[0],p1[1]), gCode(p2[0],p2[1] ),
          nG[0], nG[1], nG[2],
          pg ) );
}

// addition for fbati extension (gXe work)
// returns the hash to the group which is given by
//  <xcode(p1)+1><xcode(p2)+1><nAA><nAB><nBB>, where p1<=p2 for missing parents,
//    since the number of each genotype affects the strata of sufficient statistics
// or when there is parents, then it should just be
//  <xcode(p1)+1><xcode(p2)+1>000000
int pG_group( int n,
              int *p1, int *p2, // parental alleles
              int *ca, int *cb, // childrens alleles
              double pg[3] ) // output parameter, prob of each genotype
{
  // copying from above
  int nG[3] = {0,0,0};
  for( int i=0; i<n; i++ )
    nG[ gCode( ca[i], cb[i] ) ]++;

  int gP1 = gCode(p1[0],p1[1]);
  int gP2 = gCode(p2[0],p2[1]);

  // fill in the genotypes
  pG( gP1, gP2,  nG[0], nG[1], nG[2],  pg );

  /*
  cout << " nAA=" << nG[gAA]
       << " nAB=" << nG[gAB]
       << " nBB=" << nG[gBB]
       << " gP1=" << gP1
       << " gP2=" << gP2
       << "  group=";
  */

  // gMiss=-1, so need to add 1 to parents code to hash correctly
  int pCode = (int)( (gP2+1)*1e6 + (gP1+1)*1e7 );
  if( gP1>gP2 )
    pCode = (int)( (gP1+1)*1e6 + (gP2+1)*1e7 );
  if( gP2 != gMiss  &&  gP1 != gMiss )  // update 08/06/07 for one missing parent...
    return( pCode ); // both parents are available!
  return( (int)(nG[0]*1e0 + nG[1]*1e2 + nG[2]*1e4) + pCode );
}

int extractDigitRHS( int number, int digit ) {
  for( int i=0; i<digit; i++ )
    number /= 10;
  int temp = number / 10;
  number -= temp*10;

  return( number );
}
void pG_group_dehash_gstr( int g, char *str ) {
  if( g==(gAA+1) )
    sprintf( str, "AA" );
  else if( g==(gAB+1) )
    sprintf( str, "AB" );
  else if( g==(gBB+1) )
    sprintf( str, "BB" );
  else
    sprintf( str, "?" );
}
extern "C" {
  void pG_group_dehash( int* num, char **str ) {
    int number = *num;

    int p1 = extractDigitRHS( number, 7 );
    int p2 = extractDigitRHS( number, 6 );

    //cout << "p1, p2 " << p1 << " " << p2 << endl;

    int n[3];
    for( int j=0; j<3; j++ )
      n[j] = extractDigitRHS(number,j*2) + 10*extractDigitRHS(number,j*2+1);

    //cout << "n" << n[0] << " " << n[1] << " " << n[2] << endl;
    //cout << number << endl << endl;

    char p1str[3], p2str[3];
    pG_group_dehash_gstr( p1, p1str );
    pG_group_dehash_gstr( p2, p2str );

    if( p1!=0 && p2!=0 )
      sprintf( *str, "%s,%s", p1str,p2str );
    else
      sprintf( *str, "%s,%s - AA%i AB%i BB%i", p1str,p2str, n[gAA], n[gAB], n[gBB] );
  }
}


bool pGG( int gP1, int gP2, // parental mating type
          int nAA, int nAB, int nBB, // number of offspring with said genotype
          double pgg[9] ) // output parameter, prob of each genotype
{
  // swap gP1 and gP2 if gP1 is missing
  if( gP1==gMiss ) {
    gP1 = gP2;
    gP2 = gMiss;
  }

  double n = nAA + nAB + nBB;
  //cout << "pGG n=" << n << " nAA=" << nAA << " nAB=" << nAB << " nBB=" << nBB << endl;

  // NOTE: unspecified genotypes are 0!
  for( int i=0; i<9; i++ )
    pgg[i] = 0.0;

  if( gP2!=gMiss ) {
    // neither parent is missing
    double pg[3];
    pG( gP1, gP2,  nAA, nAB, nBB,  pg );
    pgg[gAAgAA] = pg[gAA]*pg[gAA];
    pgg[gAAgAB] = pg[gAA]*pg[gAB];
    pgg[gAAgBB] = pg[gAA]*pg[gBB];
    pgg[gABgAA] = pg[gAB]*pg[gAA];
    pgg[gABgAB] = pg[gAB]*pg[gAB];
    pgg[gABgBB] = pg[gAB]*pg[gBB];
    pgg[gBBgAA] = pg[gBB]*pg[gAA];
    pgg[gBBgAB] = pg[gBB]*pg[gAB];
    pgg[gBBgBB] = pg[gBB]*pg[gBB];
    return(true);
  }

  if( gP1!=gMiss ) {
    // just one parent is missing
    switch( gP1 ) {
    case gBB: // they are the same...
    case gAA:
      // one homozygous parent - all good under BB
      if( first(nAA,nAB,nBB) ) {
        // only happens under BB
        pgg[gAAgAA] = 1.00;
        return(true);
      }
      if( second(nAA,nAB,nBB) ) {
        pgg[gABgAB] = 1.00;
        return(true);
      }
      if( third(nAA,nAB,nBB) ) {
        // only happens under AA
        pgg[gBBgBB] = 1.00;
        return(true);
      }

      if( first_second(nAA,nAB,nBB) ) {
        // only happens under AA
        pgg[gAAgAB] = pow(2,n-2)/(pow(2,n)-2);
        pgg[gAAgAA] = (pow(2,n-2)-1)/(pow(2,n)-2);
        pgg[gABgAB] = pgg[gAAgAA];

        pgg[gABgAA] = pgg[gAAgAB];
        return(true);
      }
      if( first_third(nAA,nAB,nBB) ) {
        Rprintf("Impossible genotypes, 1 missing parent.\n");
        return(false);
      }
      if( second_third(nAA,nAB,nBB) ) {
        // only happens under BB
        pgg[gBBgAB] = pow(2,n-2)/(pow(2,n)-2);
        pgg[gBBgBB] = (pow(2,n-2)-1)/(pow(2,n)-2);
        pgg[gABgAB] = pgg[gBBgBB];

        pgg[gABgBB] = pgg[gBBgAB];
        return(true);
      }
    }
  }

  // otherwise both parents are missing,
  //  or the parent is heterozygous
  // either are the same
  if( gP1==gAB || gP1==gMiss ) {
    // this is pretty much a copy of the lower part of the above...
    if( first(nAA,nAB,nBB) ) {
      pgg[gAAgAA] = 1.00;
      return(true);
    }
    if( second(nAA,nAB,nBB) ) {
      pgg[gABgAB] = 1.00;
      return(true);
    }
    if( third(nAA,nAB,nBB) ) {
      pgg[gBBgBB] = 1.00;
      return(true);
    }

    if( first_second(nAA,nAB,nBB) ) {
      pgg[gAAgAA] = ( nAA * (nAA-1) ) / (n * (n-1));
      pgg[gABgAB] = ( nAB * (nAB-1) ) / (n * (n-1));
      pgg[gAAgAB] = ( nAA * nAB ) / (n * (n-1));

      pgg[gABgAA] = pgg[gAAgAB];
      return(true);
    }
    if( first_third(nAA,nAB,nBB) || all(nAA,nAB,nBB) ) {
      pgg[gAAgAA] = (pow(4,n-2)-pow(3,n-2)) / (pow(4,n)-2*pow(3,n)+pow(2,n));
      pgg[gBBgBB] = pgg[gAAgAA];
      pgg[gAAgAB] = 2.0 * pgg[gAAgAA];
      pgg[gBBgAB] = 2.0 * pgg[gAAgAA];

      pgg[gAAgBB] = pow(4,n-2) / (pow(4,n)-2*pow(3,n)+pow(2,n));
      pgg[gABgAB] = ( pow(4,n-1) - 8*pow(3,n-2) + pow(2,n) ) / (pow(4,n)-2*pow(3,n)+pow(2,n));

      pgg[gABgAA] = pgg[gAAgAB];
      pgg[gABgBB] = pgg[gBBgAB];
      pgg[gBBgAA] = pgg[gAAgBB];
      return(true);
    }
    if( second_third(nAA,nAB,nBB) ) {
      pgg[gBBgBB] = ( nBB * (nBB-1) ) / (n * (n-1));
      pgg[gABgAB] = ( nAB * (nAB-1) ) / (n * (n-1));
      pgg[gBBgAB] = ( nBB * nAB ) / (n * (n-1));;

      pgg[gABgBB] = pgg[gBBgAB];
      return(true);
    }
  }

  return(false);
}

bool pGG( int n,
          int *p1, int *p2, // parental alleles
          int *ca, int *cb, // childrens alleles
          double pgg[9] ) { // output parameter, prob of each genotype

  //cout << "pGG n=" << n << " p1=" << p1[0] << p1[1] << " p2=" << p2[0] << p2[1];
  //for( int i=0; i<n; i++ ) cout << "c" << i << "=" << ca[i] << cb[i];
  //cout << endl;

  int nG[3] = {0,0,0};
  for( int i=0; i<n; i++ )
    if( ca[i]!=0 && cb[i]!=0 ) // NEW 06.24.2009
      nG[ gCode( ca[i], cb[i] ) ]++;

  return( pGG( gCode(p1[0],p1[1]), gCode(p2[0],p2[1]),
          nG[0], nG[1], nG[2],
          pgg ) );
}

// Addition 02/11/2008
//  computes E(X|S) [[which is constant for a family]]
double fbat_EXS( int n,
                 int *p1, int *p2, // parental alleles
                 int *ca, int *cb, // childrens alleles
                 double *y,           // childrens trait
                 int model        // genetic model (a/d/r)
               ) {
  /*
  cout << "fbat_EXS n=" << n << " model=" << model << endl;
  cout << " p1[0]=" << p1[0] << " p1[1]=" << p1[1] << endl;
  cout << " p2[0]=" << p2[0] << " p2[1]=" << p2[1] << endl;
  int i;
  for( i=0; i<n; i++ ) cout << " ca[" << i << "]=" << ca[i];
  cout << endl;
  for( i=0; i<n; i++ ) cout << " cb[" << i << "]=" << cb[i];
  cout << endl;
  for( i=0; i<n; i++ ) cout << " y[" << i << "]=" << y[i];
  cout << endl;
  */

  double pg[3];
  if( !pG( n,  p1, p2,  ca, cb,  pg ) ) {
    //cout << "really did fail..." << endl;
    //printFamily( p1, p2, ca, cb, n );
    // Happens when the kids are all missing
    return( 0.0 );
  }

  double exj = 0;
  int g;
  for( g=0; g<3; g++ )
    exj += pg[g] * xCode( g, model );

  return( exj );
}

// computes sum_j X-E(X|S)
// - i.e. family-wise contribution to fbat statistic
double fbat_Si( int n,
                int *p1, int *p2, // parental alleles
                int *ca, int *cb, // childrens alleles
                double *y,           // childrens trait
                int model,        // genetic model (a/d/r)
                double &fbat_Vi, // variance of what calculating
                double offset,
                int nPhenotyped ) { // hack for power
  // get rid of individuals who are missing a trait (y) or have bad genotypes; this is the easy way to do this!:
  int ngenopheno = 0;
  //double y_genopheno[n]; int ca_genopheno[n], cb_genopheno[n];
  vector<double> y_genopheno; y_genopheno.resize(n);
  vector<int> ca_genopheno; ca_genopheno.resize(n);
  vector<int> cb_genopheno; cb_genopheno.resize(n);
  int ngeno = 0;
  //int ca_geno[n], cb_geno[n];
  vector<int> ca_geno; ca_geno.resize(n);
  vector<int> cb_geno; cb_geno.resize(n);
  for( int j=0; j<n; j++ ) {
    if( ca[j]!=0 && cb[j]!=0 ) {
      if( !ISNAN(y[j]) ) {
        y_genopheno[ngenopheno] = y[j];
        ca_genopheno[ngenopheno] = ca[j];
        cb_genopheno[ngenopheno] = cb[j];
        ngenopheno++;
      }else{
        ca_geno[ngeno] = ca[j];
        cb_geno[ngeno] = cb[j];
        ngeno++;
      }
    }
  }//j
  if( ngenopheno==0 ) {
    fbat_Vi = 0.0;
    return(0.0);
  }
  for( int j=0; j<ngenopheno; j++ ) {
    y[j] = y_genopheno[j];
    ca[j] = ca_genopheno[j];
    cb[j] = cb_genopheno[j];
  }
  for( int j=0; j<ngeno; j++ ) {
    y[j+ngenopheno] = -99999; //offset; // debatable.... -- no, really screw things up if used
    ca[j+ngenopheno] = ca_geno[j];
    cb[j+ngenopheno] = cb_geno[j];
  }
  n = ngenopheno + ngeno;
// end of get rid of bad individuals

  double pg[3];
  if( !pG( n,  p1, p2,  ca, cb,  pg ) ) {
    Rprintf("really did fail...\n");
    fbat_Vi = 0.0;
    return( 0.0 );
  }

  // otherwise we've got multiple offspring
  // -- ugh!!! This needs to be before we then chop off the pheno... ahhhHHHH!!!
  double pgg[9];
  if( n>1 && nPhenotyped>1 )
    pGG( n,  p1, p2,  ca, cb,  pgg );

  n = ngenopheno; // NEW NEW OY VAY!!!!

  double exj = 0;
  int g;
  for( g=0; g<3; g++ )
    exj += pg[g] * xCode( g, model );

  double Si = 0.0;
  int j;
  for( j=0; j<n && j<nPhenotyped; j++ ) {
    // (X-E(X|S))*Y
    Si += ( xCode( ca[j], cb[j], model ) - exj ) * (y[j]-offset);
    //cout << "Si=" << Si << endl;
  }

  // now the variance
  // we've got a special case if there is only one child... (power hack)
  if( n==1 || nPhenotyped==1 ) {
    fbat_Vi = 0.0;
    // E(X^2)
    for( int g=0; g<3; g++ ) {
      double x = xCode(g,model);
      fbat_Vi += x*x*pg[g];
    }
    // - E(X)^2
    fbat_Vi -= exj*exj;
    fbat_Vi *= (y[0]-offset)*(y[0]-offset); // and don't loose the trait!

    if( n==1 || nPhenotyped==1 ) {
      return( Si );
    }
  }

  // used to calculate Pgg here, but that needed to be moved to the top!!!
  //cout << "Pgg" << endl
  //  << pgg[0] << " " << pgg[1] << " " << pgg[2] << endl
  //  << pgg[3] << " " << pgg[4] << " " << pgg[5] << endl
  //  << pgg[6] << " " << pgg[7] << " " << pgg[8] << endl;

#ifdef DEBUG_FBATDIST
  double pgg_sum=0.0;
  int gg;
  for( gg=0; gg<9; gg++ ) pgg_sum += pgg[gg];
  if( pgg_sum<0.99 || pgg_sum>1.01 ) {
    printFamily( p1, p2,  ca, cb,  n );
    Rprintf("pgg_sum = %f\n", pgg_sum);
    for( gg=0; gg<9; gg++ )
      Rprintf(" P[%g]=%g\n", gg, pgg[gg]);
    exit(1);
  }
#endif

  // The coding from the paper works fine for the variance
  //  so long as there is _more_ than one offspring!
  double sumTj=0;
  for( j=0; j<n && j<nPhenotyped; j++ )
    sumTj += (y[j]-offset);

  // the first term
  fbat_Vi = 0.0;
  int g1, g2;
  for( g1=0; g1<3; g1++ )
    for( g2=0; g2<3; g2++ )
      fbat_Vi += xCode(g1,model)*xCode(g2,model)
        *( pgg[ggConvert(g1,g2)] - pg[g1]*pg[g2] );
  fbat_Vi *= sumTj*sumTj;


  // the second term
  for( j=0; j<n && j<nPhenotyped; j++ ) {
    double sum = 0.0; // jth term sum over g, g~
    for( g1=0; g1<3; g1++ ) {
      sum += pow((double)xCode(g1,model),2)*pg[g1];
      for( g2=0; g2<3; g2++ ) {
        sum -= xCode(g1,model)*xCode(g2,model) * pgg[ggConvert(g1,g2)];
      }
    }
    sum *= (y[j]-offset) * (y[j]-offset);

    fbat_Vi += sum; // ah... there was the rub, this coding works just fine...
  }

  // and that computes the variance!
  //fbat_Vi = Si * Si; // just for a test...

  //cout << "fbatdist.cpp Si=" << Si << " Vi=" << fbat_Vi << endl;

  return( Si );
}


// A rewrite specially formatted to calculate the joint G, GxE case
void fbat_Si_joint_G_GE( int n,
                         int *p1, int *p2, // parental alleles
                         int *ca, int *cb, // childrens alleles
                         double *y, double *z,  // childrens trait, environmental information
                         int model,        // genetic model (a/d/r)
                         double &Si0, double &Si1, // sum Tj(Xj - E[Xj|s])
                         double &Vi00, double &Vi01, double &Vi10, double &Vi11, // variance of what calculating
                         double offset,
                         int nPhenotyped ) { // hack for power
  Si0 = Si1 = Vi00 = Vi01 = Vi10 = Vi11 = 0.0;

  // get rid of individuals who are missing a trait (y) or have bad genotypes; this is the easy way to do this!:
  int ngenopheno = 0;
  //double y_genopheno[n], z_genopheno[n]; int ca_genopheno[n], cb_genopheno[n];
  vector<double> y_genopheno; y_genopheno.resize(n);
  vector<double> z_genopheno; z_genopheno.resize(n);
  vector<int> ca_genopheno; ca_genopheno.resize(n);
  vector<int> cb_genopheno; cb_genopheno.resize(n);
  int ngeno = 0;
  //int ca_geno[n], cb_geno[n];
  vector<int> ca_geno; ca_geno.resize(n);
  vector<int> cb_geno; cb_geno.resize(n);
  for( int j=0; j<n; j++ ) {
    if( ca[j]!=0 && cb[j]!=0 ) {
      if( !ISNAN(y[j]) && !ISNAN(z[j]) ) {
        y_genopheno[ngenopheno] = y[j];
        z_genopheno[ngenopheno] = z[j];
        ca_genopheno[ngenopheno] = ca[j];
        cb_genopheno[ngenopheno] = cb[j];
        ngenopheno++;
      }else{
        ca_geno[ngeno] = ca[j];
        cb_geno[ngeno] = cb[j];
        ngeno++;
      }
    }
  }//j
  if( ngenopheno==0 ) {
    // everything was already zeroed out on entry
    return;
  }
  for( int j=0; j<ngenopheno; j++ ) {
    y[j] = y_genopheno[j];
    z[j] = z_genopheno[j];
    ca[j] = ca_genopheno[j];
    cb[j] = cb_genopheno[j];
  }
  for( int j=0; j<ngeno; j++ ) {
    y[j+ngenopheno] = -99999; //offset; // debatable.... -- no, really screw things up if used
    ca[j+ngenopheno] = ca_geno[j];
    cb[j+ngenopheno] = cb_geno[j];
  }
  n = ngenopheno + ngeno;
  // end of get rid of bad individuals

  double pg[3];
  if( !pG( n,  p1, p2,  ca, cb,  pg ) ) {
    Rprintf("really did fail...\n");
    // everything was already zeroed out
    return;
  }

  // otherwise we've got multiple offspring
  // -- ugh!!! This needs to be before we then chop off the pheno... ahhhHHHH!!!
  double pgg[9];
  if( n>1 && nPhenotyped>1 )
    pGG( n,  p1, p2,  ca, cb,  pgg );

  n = ngenopheno; // NEW NEW OY VAY!!!!

  double exj = 0;
  int g;
  for( g=0; g<3; g++ )
    exj += pg[g] * xCode( g, model );

  int j;
  for( j=0; j<n && j<nPhenotyped; j++ ) {
    // (X-E(X|S))*Y
    double temp = ( xCode( ca[j], cb[j], model ) - exj );
    Si0 += (y[j]-offset) * temp;
    Si1 += (y[j]-offset) * z[j] * temp;
  }

  // now the variance
  // we've got a special case if there is only one child... (power hack)
  if( n==1 || nPhenotyped==1 ) {
    double Vi = 0.0;
    // E(X^2)
    for( int g=0; g<3; g++ ) {
      double x = xCode(g,model);
      Vi += x*x*pg[g];
    }
    // - E(X)^2
    Vi -= exj*exj;
    Vi00 = Vi * (y[0]-offset)*(y[0]-offset);
    Vi01 = Vi * (y[0]-offset)*(y[0]-offset) * z[0];
    Vi10 = Vi01;
    Vi11 = Vi * (y[0]-offset)*(y[0]-offset) * z[0] * z[0];

    if( n==1 || nPhenotyped==1 ) {
      return;
    }
  }

  // Removed DEBUG_FBAT piece (it's debugged in the main effects test...)
/*
  // OK, completely rewrite this piece -- ignore the formula in the paper,
  //  and just do it the brute force way...
  double VXi = 0.0;
  double covXi = 0.0;
  double E_Xi = 0.0;
  double E_Xisq = 0.0;

  for( int g=0; g<3; g++ ) {
    double x = xCode(g,model);
    E_Xi += x * pg[g];
    E_Xisq += x * x * pg[g];
  }
  VXi = E_Xisq - E_Xi*E_Xi;

  double E_XiXj = 0.0;
  for( int g1=0; g1<3; g1++ )
    for( int g2=0; g2<3; g2++ )
      E_XiXj += xCode(g1,model) * xCode(g2,model) * pgg[ggConvert(g1,g2)];
  covXi = E_XiXj - E_Xi*E_Xi;

  for( int j0=0; j0<n && j0<nPhenotyped; j0++ ) {
    for( int j1=0; j1<n && j1<nPhenotyped; j1++ ) {
      double temp = 0.0;
      if( j0==j1 ) {
        temp = VXi;
      }else{
        temp = covXi;
      }

      double Tij0 = (y[j0] - offset) * (y[j1] - offset);
      double Tij1 = Tij0 * z[j0] * z[j1];

      Vi00 += Tij0 * Tij0 * temp;
      Vi01 += Tij0 * Tij1 * temp;
      Vi10 += Tij1 * Tij0 * temp;
      Vi11 += Tij1 * Tij1 * temp;
    }
  }

  // and that should do it???
*/

  // The coding from the paper works fine for the variance
  //  so long as there is _more_ than one offspring!
  /*
  double sumTj00=0;
  double sumTj01=0;
  double sumTj10=0;
  double sumTj11=0;
  for( j=0; j<n && j<nPhenotyped; j++ ) {
    double temp = (y[j]-offset) * (y[j]-offset);
    sumTj00 += temp;
    sumTj01 += temp * z[j];
    sumTj10 += temp * z[j];
    sumTj11 += temp * z[j] * z[j];
  }

  // the first term
  double Vi = 0.0;
  int g1, g2;
  for( g1=0; g1<3; g1++ )
    for( g2=0; g2<3; g2++ )
      Vi += xCode(g1,model)*xCode(g2,model)
        *( pgg[ggConvert(g1,g2)] - pg[g1]*pg[g2] );
  //fbat_Vi *= sumTj*sumTj;
  Vi00 = sumTj00 * Vi;
  Vi01 = sumTj01 * Vi;
  Vi10 = sumTj10 * Vi;
  Vi11 = sumTj11 * Vi;
  */



  // The coding ... redo!
  double sumTj0 = 0.0;
  double sumTj1 = 0.0;
  for( j=0; j<n && j<nPhenotyped; j++ ) {
    sumTj0 += (y[j]-offset);
    sumTj1 += (y[j]-offset) * z[j];
  }
  double Vi = 0.0;
  int g1, g2;
  for( g1=0; g1<3; g1++ )
    for( g2=0; g2<3; g2++ )
      Vi += xCode(g1,model)*xCode(g2,model)
        *( pgg[ggConvert(g1,g2)] - pg[g1]*pg[g2] );
  Vi00 = Vi * sumTj0 * sumTj0;
  Vi01 = Vi * sumTj0 * sumTj1;
  Vi10 = Vi * sumTj1 * sumTj0;
  Vi11 = Vi * sumTj1 * sumTj1;

  // the second term
  for( j=0; j<n && j<nPhenotyped; j++ ) {
    double sum = 0.0; // jth term sum over g, g~ -- really x-exs
    for( g1=0; g1<3; g1++ ) {
      sum += pow((double)xCode(g1,model),2)*pg[g1];
      for( g2=0; g2<3; g2++ ) {
        sum -= xCode(g1,model)*xCode(g2,model) * pgg[ggConvert(g1,g2)];
      }
    }

    double temp = (y[j]-offset) * (y[j]-offset);
    double sum00 = sum * temp;
    double sum01 = sum * temp * z[j];
    double sum10 = sum * temp * z[j];
    double sum11 = sum * temp * z[j] * z[j];

    Vi00 += sum00;
    Vi01 += sum01;
    Vi10 += sum10;
    Vi11 += sum11;
  }


  // and that computes the variance!
}




#ifdef LOOKUP_COMPARE
bool fuzzyEqual( double a, double b ) {
  return( fabs(a-b) <= 0.0001 );
}

// Make sure this works with our lookup table created from fbat
void recursiveFillLookupCompare( int curSib, int numSibs,
                          int *p1, int *p2,
                          int *ca, int *cb )
{
  if( curSib >=1 ){
    // need to do a recursive call

    for( int cai=1; cai<3; cai++ ){
      for( int cbi=1; cbi<3; cbi++ ){
        ca[curSib-1] = cai;
        cb[curSib-1] = cbi;

        recursiveFillLookupCompare( curSib-1, numSibs, p1, p2, ca, cb );
      }
    }
  }else{
    // needs to be computed

    // first by the new code
    double pg[3];
    if( !pG( numSibs,  p1, p2,  ca, cb,  pg ) )
      printFamily( p1, p2, ca, cb, numSibs );

    // now by the lookup
    int index = indexLookup( p1, p2,  ca, cb,  numSibs );
    if( chart_g1(index) != -1 ) {
      // then it's an informative family
      if( !fuzzyEqual( chart_g1(index), pg[0] ) ||
          !fuzzyEqual( chart_g2(index), pg[1] ) ||
          !fuzzyEqual( chart_g3(index), pg[2] ) ) {
        //cout << "Lookup failure! You: " << pg[0] << " " << pg[1] << " " << pg[2] << "Fbat: " << chart_g1(index) << " " << chart_g2(index) << " " << chart_g3(index) << endl;
        Rprintf("Lookup failure! You: %g %g %g, FBAT: %g %g %g\n", pg[0], pg[1], pg[2], chart_g1(index), chart_g2(index), chart_g3(index));
        printFamily( p1, p2,  ca, cb,  numSibs );
      }

    }else{
      // need to think about this, get above first...
    }
  }

  // fell out of all of the normal cases...
  if( (gP1==gAA||gP1==gBB) && gP2==gMiss
      && all(nAA,nAB,nBB) ) {
#ifndef LOOKUP_COMPARE
    Rprintf("WARNING: impossible genotype in file.\n");
    printFamily();
#endif
    return(true);
  }

  Rprintf("failed all cases!\n");
  Rprintf(" %g %g\n", gP1, gP2);
  return(false);
}

void fillLookupCompare()
{
  int i;

  int N=4; // temporary...
  Rprintf("Table size = %d\n", N);

  // now begin the recursive fill
  // - setup
  int p1[2], p2[2], ca[N], cb[N];
  for( int numSibs=1; numSibs<=N; numSibs++ ){
    for( int p1ai=0; p1ai<3; p1ai++ ){
      for( int p1bi=0; p1bi<3; p1bi++ ){
        if( p1ai==0 && p1bi!=0 ) continue;
        if( p1ai!=0 && p1bi==0 ) continue;

        p1[0] = p1ai;
        p1[1] = p1bi;
        for( int p2ai=0; p2ai<3; p2ai++ ){
          for( int p2bi=0; p2bi<3; p2bi++ ){
            if( p2ai==0 && p2bi!=0 ) continue;
            if( p2ai!=0 && p2bi==0 ) continue;

            p2[0] = p2ai;
            p2[1] = p2bi;

            recursiveFillLookupCompare( numSibs, numSibs, p1, p2, ca, cb );
          }
        }
      }
    }
  }
}


int main()
{
  MODEL model = ADDITIVE;
  setupLookup( model );
  fillLookupCompare();
}
#endif
