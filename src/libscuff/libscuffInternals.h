/*
 * libscuffInternals.h -- some definitions and prototypes that are used
 *                     -- internally within libscuff code files, but which
 *                     -- probably don't need to be exported as part of the API
 *
 * homer reid          -- 10/2005 -- 10/2011
 */

#ifndef LIBSCUFFINTERNALS_H
#define LIBSCUFFINTERNALS_H

#include "libscuff.h"

/***************************************************************/
/* 1. argument structures for routines whose input/output      */
/*    interface is so complicated that an ordinary C++         */
/*    function prototype would be too unwieldy.                */
/*                                                             */
/*    note: in each case, we provide an Init() function that   */
/*    fills in default values for the lesser-used fields       */
/*    in the argument structure; however, this leaves many     */
/*    fields uninitialized, and it is up to the caller to      */
/*    fill those in.                                           */
/***************************************************************/

/*--------------------------------------------------------------*/
/*- GetPanelPanelInteractions() --------------------------------*/
/*--------------------------------------------------------------*/
typedef struct GPPIArgStruct
 { 
   // input fields to be filled in by caller
   RWGObject *Oa, *Ob;
   int npa, npb;
   int iQa, iQb;
   cdouble k;

   int NumGradientComponents;
   int NumTorqueAxes; 
   double *GammaMatrix;

   // output fields filled in by routine
   // note: GC[0,1] = G, C
   // note: GradGC[0,1,2,...] = dG/dx, dC/dx, dG/dy, ... 
   // note: dGCdT[0,1,2,...] = dG/dTheta_1, dC/dTheta_1, dG/dTheta_2, ... 
   cdouble GC[2];
   cdouble GradGC[6];
   cdouble dGCdT[6];

 } GPPIArgStruct;

void InitGPPIArgs(GPPIArgStruct *Args);
void GetPanelPanelInteractions(GEEIArgStruct *Args);
void GetPanelPanelInteractions(GEEIArgStruct *Args,
                               cdouble *GC, 
                               cdouble *GradGC, 
                               cdouble *dGCdT);

/*--------------------------------------------------------------*/
/*- GetEdgeEdgeInteractions() ----------------------------------*/
/*--------------------------------------------------------------*/
typedef struct GEEIArgStruct
 { 
   // input fields to be filled in by caller
   RWGObject *Oa, *Ob;
   int nea, neb;
   cdouble k;

   int NumGradientComponents;
   int NumTorqueAxes; 
   double *GammaMatrix;

   // output fields filled in by routine
   // note: GC[0,1] = G, C
   // note: GradGC[0,1,2,...] = dG/dx, dC/dx, dG/dy, ... 
   // note: dGCdT[0,1,2,...] = dG/dTheta_1, dC/dTheta_1, dG/dTheta_2, ... 
   cdouble GC[2];
   cdouble GradGC[6];
   cdouble dGCdT[6];

 } GEEIArgStruct;

void InitGEEIArgs(GEEIArgStruct *Args);
void GetEdgeEdgeInteractions(GEEIArgStruct *Args);

/*--------------------------------------------------------------*/
/*- AssembleBEMMatrixBlock() -----------------------------------*/
/*--------------------------------------------------------------*/
typedef struct ABMBArgStruct
 {
   // input fields to be filled in by caller
   RWGGeometry *G;
   RWGObject *Oa, *Ob;
   cdouble Frequency;
   int nThread;

   int NumTorqueAxes;
   double *GammaMatrix;
  
   int RowOffset, ColOffset;

   int Symmetric;

   // output fields filled in by routine
   HMatrix *B;
   HMatrix **GradB;
   HMatrix **dBdTheta;

   // additional fields used internally that may be ignored by 
   // the called both before and after the call
   double Sign;
   cdouble EpsA, EpsB; 
   double MuA, MuA;
   int OaIsPEC, ObIsPEC;

 } ABMBArgStruct;

void InitABMBArgs(ABMBArgStruct *Args);
void AssembleBEMMatrixBlock(ABMBArgStruct *Args);

#endif //LIBSCUFFINTERNALS_H
