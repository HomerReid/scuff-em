/*
 * scuff-caspol.h    -- header file for scuff-caspol
 *
 * homer reid        -- 10/2006 -- 2/2012
 *
 */
#include <stdio.h>
#include <math.h>
#include <stdarg.h>

#include <libhrutil.h>
#include <libIncField.h>
#include <libhmat.h>

#include "libscuff.h"

// default tolerances for numerical summation and integration
#define ABSTOL 1.0e-8
#define RELTOL 1.0e-3

// the minimum frequency at which we do calculations. 
#define XIMIN 1.0e-3

/***************************************************************/
/* prototype for user-supplied routine to compute              */
/* polarizability.                                             */
/*                                                             */
/* on entry, Xi is the imaginary (angular) frequency in units  */
/* of 3e14 radians/second.                                     */
/*                                                             */
/* the routine fills in the Alpha array with the components of */
/* the polarizability tensor, as follows:                      */
/*                                                             */
/*  Alpha[ i + 3*j ] = Alpha_{ij}                              */
/*                                                             */
/* for i,j = 0,1,2.                                            */
/*                                                             */
/* thus, on return,                                            */
/*  Alpha[0] is the xx polarizability                          */
/*  Alpha[1] is the xy polarizability                          */
/*  Alpha[2] is the xz polarizability                          */
/*  Alpha[3] is the yz polarizability                          */
/* etc.                                                        */
/* (only nonzero entries need be filled in.)                   */
/***************************************************************/
typedef void (*AlphaFuncType)(double Xi, double *Alpha);

/***************************************************************/
/* SCPData ('scuff-caspol-data') is the basic structure passed */
/* around among the various routines; it contains everything   */
/* needed to compute the casimir-polder potential              */
/***************************************************************/
typedef struct SCPData
 {
   RWGGeometry *G;
   HMatrix *M;
   HVector *KN;

   MatProp *AlphaMP[9];
   AlphaFuncType AlphaFunc;

   HMatrix *EPList;

   int nThread;
   double AbsTol, RelTol;
   FILE *ByXiFile;

 } SCPData; 
