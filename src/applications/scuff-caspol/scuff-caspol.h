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
   int nThread;

 } SCPData; 
