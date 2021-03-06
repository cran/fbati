/* Thomas Hoffmann
 * fbatDist.cpp
 * Created: 5/28/2007
 * Modified: 2/11/2008
 */

/** A is the disease allele **/

// biallelic markers only

#ifndef _fbatDist_h_
#define _fbatDist_h_

// compile the lookup compare with g++ *.cpp
//#define LOOKUP_COMPARE

#ifdef LOOKUP_COMPARE
#include "tableLookup.h"
#endif

#include "random.h"

#include <math.h>
#include <iostream>
#include <cstdio>
using namespace std;

// constants
//const int gA = 0;
//const int gB = 1;

const int gMiss = -1;
const int gAA = 0;
const int gAB = 1;
const int gBB = 2;

const int gAAgAA = 0;
const int gAAgAB = 1;
const int gAAgBB = 2;
const int gABgAA = 3;
const int gABgAB = 4;
const int gABgBB = 5;
const int gBBgAA = 6;
const int gBBgAB = 7;
const int gBBgBB = 8;

const int MODEL_ADDITIVE=0;
const int MODEL_DOMINANT=1;
const int MODEL_RECESSIVE=2;
const char MODEL_FBAT_CHARS[] = {'a', 'd', 'r'};

const int DA_MA = 0;
const int DA_MB = 1;
const int DB_MA = 2;
const int DB_MB = 3;

const int ALLELE_A=2; // coding the disease allele
const int ALLELE_B=1;



// functions

int gCode( int a, int b );
int xCode( int a, int b, int MODEL );

bool pG( int gP1, int gP2, // parental mating type
         int nAA, int nAB, int nBB, // number of offspring with said genotype
         double pg[3] ); // output parameter, prob of each genotype
bool pG( int n,
         int *p1, int *p2, // parental alleles
         int *ca, int *cb, // childrens alleles
         double pg[3] ); // output parameter, prob of each genotype

bool pGG( int gP1, int gP2,
          int nAA, int nAB, int nBB,
          double pg[3], double pgg[3][3] ); // ouput parameters, prob of each genotype, and joint prob of the genotypes

// addition for fbati extension (gXe work)
// returns the hash to the group which is given by
//  <xcode(p1)+1><xcode(p2)+1><nAA><nAB><nBB>, where p1<=p2
int pG_group( int n,
              int *p1, int *p2, // parental alleles
              int *ca, int *cb, // childrens alleles
              double pg[3] ); // output parameter, prob of each genotype
extern "C" {
  void pG_group_dehash( int* num, char** str );
}

bool pGG( int gP1, int gP2, // parental mating type
          int nAA, int nAB, int nBB, // number of offspring with said genotype
          double pgg[9] ); // output parameter, prob of each genotype
bool pGG( int n,
          int *p1, int *p2, // parental alleles
          int *ca, int *cb, // childrens alleles
          double pgg[9] ); // output parameter, prob of each genotype

double fbat_EXS( int n,
                 int *p1, int *p2, // parental alleles
                 int *ca, int *cb, // childrens alleles
                 double *y,           // childrens trait
                 int model        // genetic model (a/d/r)
               );

double fbat_Si( int n,
                int *p1, int *p2, // parental alleles
                int *ca, int *cb, // childrens alleles
                double *y,           // childrens trait
                int model,        // genetic model (a/d/r)
                double &fbat_Vi,   // variance of what calculating
                double offset,
                int nPhenotyped  ); // hack for pbatR power sims

void fbat_Si_joint_G_GE( int n,
                         int *p1, int *p2, // parental alleles
                         int *ca, int *cb, // childrens alleles
                         double *y, double *z,  // childrens trait, environmental information
                         int model,        // genetic model (a/d/r)
                         double &Si0, double &Si1, // sum Tj(Xj - E[Xj|s])
                         double &Vi00, double &Vi01, double &Vi10, double &Vi11, // variance of what calculating
                         double offset,
                         int nPhenotyped ); // hack for power

#endif
