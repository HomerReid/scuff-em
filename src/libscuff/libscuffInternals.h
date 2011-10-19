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
   RWGObject *O1, *O2;
   int ne1, ne2;
   cdouble Frequency;

   int NumGradientComponents;
   int NumTorqueAxes; 
   double *GammaMatrix;

   // output fields filled in by routine
   cdouble L[3]; 
   cdouble GradL[3][3];
   cdouble dLdT[3][3];

 } GPPIArgStruct;

void InitGPPIArgs(GPPIArgStruct *Args);
void GetPanelPanelInteractions(GEEIArgStruct *Args);

/*--------------------------------------------------------------*/
/*- GetEdgeEdgeInteractions() ----------------------------------*/
/*--------------------------------------------------------------*/
typedef struct GEEIArgStruct
 { 
   // input fields to be filled in by caller
   RWGObject *Oa, *Ob;
   int nea, neb;
   cdouble K;

   int NumGradientComponents;
   int NumTorqueAxes; 
   double *GammaMatrix;

   // output fields filled in by routine
   cdouble GInt, CInt;
   cdouble GradGInt[3], GradCInt[3];
   cdouble dGIntdTheta[3], dCIntdTheta[3];

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
