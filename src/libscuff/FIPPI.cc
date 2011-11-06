/*
 * FIPPI.cc     -- libscuff routines for working with frequency-independent
 *                 panel-panel integrals
 * 
 * homer reid   -- 11/2005 -- 11/2011
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <libhrutil.h>
#include <libTriInt.h>

#include "libscuff.h"
#include "libscuffInternals.h"

/***************************************************************/
/* implementation of the FIPPIDT class *************************/
/***************************************************************/

/*--------------------------------------------------------------*/
/*- class constructor ------------------------------------------*/
/*--------------------------------------------------------------*/
FIPPIDataTable()
{
}

/*--------------------------------------------------------------*/
/*- class destructor  ------------------------------------------*/
/*--------------------------------------------------------------*/
~FIPPIDataTable()
{
} 

/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
void ComputeFIPPIDataRecord_Cubature(double **VA, double *QA,
                                     double **VB, double *QB,
                                     int NeedDerivatives,
                                     FIPPIDataRecord *FDR)
{
  FDR->HaveDerivatives=NeedDerivatives;
}

/*--------------------------------------------------------------*/
/*- routine for computing frequency-independent panel-panel     */
/*- integrals                                                   */
/*-                                                             */
/*- inputs:                                                     */
/*-                                                             */
/*-  Va[i][j] = jth cartesian coord of ith vertex of panel A    */
/*-  Qa[i]    = ith cartesian coordinate of current source/sink */
/*-             vertex of panel A                               */
/*-  Vb, Qb   = similarly for panel B                           */
/*-  NeedDerivatives = set to 1 or 0                            */
/*-                                                             */
/*--------------------------------------------------------------*/
void ComputeFIPPIDataRecord(double **VA, double *QA, 
                            double **VB, double *QB,
                            int NeedDerivatives
                            FIPPIDataRecord *FDR)
{ 
  ncv=AssessPanelPair(Va, Vb);

  /*--------------------------------------------------------------*/
  /*- if there are no common vertices, then use 4-dimensional    -*/
  /*- adaptive cubature over both triangles to compute the FIPPIs-*/
  /*--------------------------------------------------------------*/
  if (ncv==0)
   { ComputeFIPPIDataRecord_Cubature(VA, QA, VB, QB,
                                     NeedDerivatives,
                                     FDR);
     return;
   };

  FDR->HaveDerivatives=0; // no derivatives for common-vertex cases 
  /*--------------------------------------------------------------*/
  /*- otherwise (there are common vertices) compute the FIPPIs   -*/
  /*- using the taylor-duffy method                              -*/
  /*--------------------------------------------------------------*/
  int WhichCase;
  switch (ncv)
   { case 1: WhichCase=TM_COMMONVERTEX; break;
     case 2: WhichCase=TM_COMMONEDGE;   break;
     case 3: WhichCase=TM_COMMONPANEL;  break; // should never happen
   };
  
  /*--------------------------------------------------------------*/
  /*- hDot integrals ---------------------------------------------*/
  /*--------------------------------------------------------------*/
  FDR->hDotRM1=real( TaylorMaster(WhichCase, TM_RP, TM_DOT, -1.0, 
                                  VA[0], VA[1], VA[2], 
                                  VB[1], VB[2], QA, QB);
                   );

  /* FIXME this one can be computed analytically */ 
  FDR->hDotR0 =real( TaylorMaster(WhichCase, TM_RP, TM_DOT, 0.0, 
                                  VA[0], VA[1], VA[2], 
                                  VB[1], VB[2], QA, QB);
                   );

  FDR->hDotR1 =real( TaylorMaster(WhichCase, TM_RP, TM_DOT, 1.0, 
                                  VA[0], VA[1], VA[2], 
                                  VB[1], VB[2], QA, QB);
                   );

  /* FIXME this one can be computed analytically */ 
  FDR->hDotR2 =real( TaylorMaster(WhichCase, TM_RP, TM_DOT, 2.0, 
                                  VA[0], VA[1], VA[2], 
                                  VB[1], VB[2], QA, QB);
                   );
  
}
