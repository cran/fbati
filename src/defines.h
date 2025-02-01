/* This file contains all the major options that can be set
 *  so that we maintain a central repository of them.
 * defines.cpp auto-generates a compatible R-file based on this.
 */

#include <R.h>

#ifndef _ge_defines_h_
#define _ge_defines_h_

// New addition here because sunCC is a bitch
#ifndef NAN
#define NAN numeric_limits<float>::quiet_NaN()
#endif

/** For debugging **/
#define DEBUG_RMATRIX
#define DEBUG_RVECTOR


/** setup for globals **/
/*
#ifdef _genall_cpp_
#define EXTERN
#define DATA(x) = x
#else
#define EXTERN extern
#define DATA(x)
#endif
*/


// Using GSL typically allows us to accelerate some
//  functionality (root finding especially),
//  but we can't assume an R user has it
// Woohoo! We figured out how to throw this in the makefile!
//#define USE_GSL


// an example: const double FOO = 13;

/** rmatrix.h - for debug only (but fantastic for that)! **/

//#define DEBUG_RMATRIX

#define GEN_TOL 1e7


// fbatdist.h

/** gen*.h **/
#define GEN_NORM_TRUNC 3.0

// under Witte's logit model, E(Z_1) \ne E(Z_2) ?
// - this doesn't really make sense,
//    and could be causing some oddities with multiple affected?
// - but if you want to see what it does instead...
//#define COVARIATE_MODEL_WITTE
// Bahadur's model isn't really a 'fix' -- doesn't work for highly corr vars...
//#define COVARIATE_MODEL_FIX_USE_BAHADUR
// However, going back to using the odds ratio
//  really should be a true fix to a lot of our problems..
// SUGGESTED
#define COVARIATE_MODEL_FIX_USE_OR

// similarly, use it for the ddata generation
// (as opposed to using the bahadur fix...)
#define FIX_USE_OR

// defined to use the logistic probability of disease rather
//  than log (as per Witte); requires GSL for integration (although
//  I've written this in R before, what's the point of
//  not just using GSL?)
//#define WITTE_LOGIT

// defined to alter the cutpoint for BINARYNORMAL
// -- half affetcted, half unaffected
#define BINARYNORMAL_CUT 0
#define BINARYNORMAL_P 0.5
// -- 90% unaffected, 10% affected
//> options( digits=22 )
//> qnorm( .9 )
//[1] 1.281551565544601
//#define BINARYNORMAL_CUT 1.281551565544601
//#define BINARYNORMAL_P   0.1


enum COVGENMODEL { BINARY, NORMAL, EXPNORMAL, UNIF, CHI, BINARYNORMAL };
enum CODEC { CODE_C, CODE_MEDIANC, CODE_EXPC };
enum MODEL { ADDITIVE, DOMINANT, RECESSIVE };


#endif
