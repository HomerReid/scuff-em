/* Copyright (C) 2005-2011 M. T. Homer Reid
 *
 * This file is part of SCUFF-EM.
 *
 * SCUFF-EM is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * SCUFF-EM is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

/*
 * TaylorDuffy.cc:  implementation of the Taylor-Duffy method for 
 *                      evaluating panel-panel integrals over pairs
 *                      of panels with common vertices
 *
 * homer reid           10/2012
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <complex>

#include <libhrutil.h>
#include <libSGJC.h>

#include "libscuff.h"
#include "libscuffInternals.h"
#include "TaylorDuffy.h"

#ifdef USE_PTHREAD
#  include <pthread.h>
#endif

namespace scuff {

#define DEFABSTOL 0.0
#define DEFRELTOL 1.0e-10
#define MAXFEVALS 10000
#define INTERVALS (MAXFEVALS/15)

#define II cdouble(0,1)

// ordering of monomials 
#define Y1TERM        1
#define Y12TERM       2
#define Y2TERM        3
#define Y1Y2TERM      4
#define Y22TERM       5
#define Y13TERM       6
#define Y12Y2TERM     7
#define Y1Y22TERM     8
#define Y23TERM       9
#define Y12Y22TERM    10 
#define Y13Y2TERM     11 
#define Y13Y22TERM    12 

#define Y3TERM        13
#define Y1Y3TERM      14
#define Y2Y3TERM      15
#define Y12Y3TERM     16
#define Y1Y2Y3TERM    17
#define Y22Y3TERM     18
#define Y12Y2Y3TERM   19

#define Y1Y22Y32TERM  20
#define Y22Y32TERM    21 

#define NUMMONOMIALS  22

// the maximum number of subregions is 6
#define NUMREGIONS 6

// the maximum power of w is actually p=4, so including the 
// p=0 term that's 5 possible powers
#define NUMWPOWERS  5

//
//
#define NUMYPOWERS  3

/***************************************************************/
/* Data structure containing various data passed back and      */
/* forth among taylor-duffy routines.                          */
/***************************************************************/
typedef struct TDWorkspace
 {
   /* */
   int WhichCase;
   int NumPKs;
   int *PIndex;
   int *KIndex;
   cdouble *KParam;
   int TwiceIntegrable;

   /* geometric data on triangles needed to compute X functions */
   double A2, B2, AP2, BP2, L2;
   double AdB, AdAP, AdBP, AdL;
   double BdAP, BdBP, APdBP, BPdL;

   /* coefficients of powers of y in script P functions */
   double Upsilon[NUMPS][NUMREGIONS][NUMWPOWERS][NUMMONOMIALS];
   int nMin[NUMPS], nMax[NUMPS], MaxMonomial[NUMPS];

   /* */
   bool NeedP[NUMPS];
   bool NeedK[NUMKS];

   int nCalls;

 } TDWorkspace;

/********************************************************************/
/* nitty-gritty subroutines, implemented at the bottom of this file */
/********************************************************************/
void GetAlphaBetaGamma2(TDWorkspace *TDW, const double *y, 
                        double *A, double *B, double *G2);

void GetX(TDWorkspace *TDW, const double *yVector, double *X);

void GetScriptP(TDWorkspace *TDW, int WhichP, const double *yVector, 
                double P[NUMREGIONS][NUMWPOWERS][NUMYPOWERS]);

void GtScriptJL(int WhichK, cdouble KParam,
                 double Alpha, double Beta, double Gamma,
                 int nMin, int nMax, 
                 cdouble JVector[NUMREGIONS], cdouble LVector[NUMREGIONS]);

void GetScriptK(int WhichK, cdouble KParam, double X,
                int nMin, int nMax, cdouble KVector[NUMREGIONS]);

void CMVStoUpsilon(int WhichCase, 
                   double *C, double *M, double *V, double S, 
                   double Upsilon[NUMREGIONS][NUMWPOWERS][NUMMONOMIALS]);

void ComputeGeometricParameters(TaylorDuffyArgStruct *Args, 
                                TDWorkspace *TDW);

/***************************************************************/
/***************************************************************/
/***************************************************************/
int TaylorDuffySum(unsigned ndim, const double *yVector, void *parms, 
                   unsigned nfun, double *f)
{
  (void) nfun; // unused
  (void) ndim; // unused

  /*--------------------------------------------------------------*/
  /*- extract parameters from data structure ---------------------*/
  /*--------------------------------------------------------------*/
  TDWorkspace *TDW = (TDWorkspace *)parms;
  int WhichCase       = TDW->WhichCase;
  int TwiceIntegrable = TDW->TwiceIntegrable;
  int NumPKs          = TDW->NumPKs;
  int *PIndex         = TDW->PIndex;
  int *KIndex         = TDW->KIndex;
  cdouble *KParam     = TDW->KParam;
  TDW->nCalls++;

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  int NumRegions, nOffset;
  double y, Jacobian;
  if (WhichCase==TD_COMMONTRIANGLE)
   { NumRegions=3;
     nOffset=1;
     y = (TwiceIntegrable ? 0.0 : yVector[0]);
     Jacobian=1.0;
   }
  else if (WhichCase==TD_COMMONEDGE)
   { NumRegions=6;
     nOffset=2;
     y = yVector[1];
     Jacobian=yVector[0];
   }
  else // (WhichCase==TD_COMMONVERTEX)
   { NumRegions=2;
     nOffset=3;
     y = yVector[2];
     Jacobian=yVector[1];
   };

  if (Jacobian==0.0)
   { memset(f,0,nfun*sizeof(double));
     return 0;
   };

  /*--------------------------------------------------------------*/
  /*- prefetch values of the X function (once-integrable case) or */
  /*- the Alpha, Beta, Gamma coefficients (twice-integrable case) */
  /*- for all subregions.                                         */
  /*--------------------------------------------------------------*/
  double X[NUMREGIONS], A[NUMREGIONS], B[NUMREGIONS], G2[NUMREGIONS];
  if (TwiceIntegrable)
   GetAlphaBetaGamma2(TDW, yVector, A, B, G2);
  else
   GetX(TDW, yVector, X);

  /*--------------------------------------------------------------*/
  /*- prefetch values of the scriptP vector for all P functions  -*/
  /*- we will need                                               -*/
  /*--------------------------------------------------------------*/
  double P[NUMPS][NUMREGIONS][NUMWPOWERS][NUMYPOWERS];
  for(int np=0; np<NUMPS; np++)
   if (TDW->NeedP[np])
    GetScriptP(TDW, np, yVector, P[np]);

  /*--------------------------------------------------------------*/
  /*- assemble the integrand vector by adding all subregions and  */
  /*- all n-values                                                */
  /*--------------------------------------------------------------*/
  cdouble J[NUMREGIONS][7], L[NUMREGIONS][7], K[NUMREGIONS][7];
  cdouble *Sum=(cdouble *)f;
  for(int npk=0; npk<NumPKs; npk++)
   { 
     int np = PIndex[npk];
     int nMin = TDW->nMin[ np ];
     int nMax = TDW->nMax[ np ];

     if (TwiceIntegrable)
      for(int d=0; d<NumRegions; d++)
       GetScriptJL( KIndex[npk], KParam[npk], A[d], B[d], G2[d], 
                    nMin+nOffset, nMax+nOffset, J[d], L[d]);
     else // once integrable
      for(int d=0; d<NumRegions; d++)
       GetScriptK( KIndex[npk], KParam[npk], X[d], 
                   nMin+nOffset, nMax+nOffset, K[d]);

     Sum[npk]=0.0;
     if (TwiceIntegrable)
      for(int n=nMin; n<=nMax; n++)
       for(int d=0; d<NumRegions; d++)
        Sum[npk] += P[np][d][n][0]*J[d][n+nOffset] + P[np][d][n][1]*L[d][n+nOffset];
     else // once integrable
      for(int n=nMin; n<=nMax; n++)
       for(int d=0; d<NumRegions; d++)
        Sum[npk] += (P[np][d][n][0] + y*P[np][d][n][1]) * K[d][n+nOffset];

     Sum[npk] *= Jacobian/(4.0*M_PI);

   };

  return 0;

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void InitTaylorDuffyArgs(TaylorDuffyArgStruct *Args)
{
  Args->Q=0;
  Args->QP=0;  
  Args->nHat=0;
  Args->AbsTol=0.0;     // DEFABSTOL;
  Args->RelTol=1.0e-10; // DEFRELTOL;
  Args->MaxEval=1000;  // DEFRELTOL;
  Args->ForceOnceIntegrable=0;
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void TaylorDuffy(TaylorDuffyArgStruct *Args)
{
  /***************************************************************/
  /* unpack fields from argument structure ***********************/
  /***************************************************************/
  int WhichCase    = Args->WhichCase;

  int NumPKs       = Args->NumPKs;
  int *PIndex      = Args->PIndex;
  int *KIndex      = Args->KIndex;
  cdouble *KParam  = Args->KParam;

  double AbsTol    = Args->AbsTol;
  double RelTol    = Args->RelTol;
  double MaxEval   = Args->MaxEval;

  /***************************************************************/
  /* initialize TDW structure to pass data to integrand routines */
  /***************************************************************/
  TDWorkspace MyTDW, *TDW=&MyTDW;
  TDW->WhichCase = WhichCase;
  TDW->NumPKs    = NumPKs;
  TDW->PIndex    = PIndex;
  TDW->KIndex    = KIndex;
  TDW->KParam    = KParam;

  memset(TDW->NeedP, 0, NUMPS*sizeof(bool));
  memset(TDW->NeedK, 0, NUMKS*sizeof(bool));
  for(int npk=0; npk<NumPKs; npk++)
   { if ( PIndex[npk]<0 || PIndex[npk]>=NUMPS )
      ErrExit("invalid PIndex (%i) in TaylorDuffy",PIndex[npk]);
     if ( KIndex[npk]<0 || KIndex[npk]>=NUMKS )
      ErrExit("invalid KIndex (%i) in TaylorDuffy",KIndex[npk]);
     TDW->NeedP[PIndex[npk]] = true;
     TDW->NeedK[KIndex[npk]] = true;
   };

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  bool NeednHat =    TDW->NeedP[TD_RNORMAL] 
                  || TDW->NeedP[TD_NMULLERG1]
                  || TDW->NeedP[TD_NMULLERG2]
                  || TDW->NeedP[TD_NMULLERC];
  if ( NeednHat && (Args->nHat)==0 )
   ErrExit("TaylorDuffy() called with nHat unspecified"); 

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  ComputeGeometricParameters(Args,TDW);
  
  /***************************************************************/
  /* assume we are twice integrable and check for otherwise      */
  /***************************************************************/
  int TwiceIntegrable=1;
  if (Args->ForceOnceIntegrable) TwiceIntegrable=0;
  for(int npk=0; TwiceIntegrable==1 && npk<NumPKs; npk++)
   if (KIndex[npk]==TD_HELMHOLTZ || KIndex[npk]==TD_GRADHELMHOLTZ)
    TwiceIntegrable=0;
  TDW->TwiceIntegrable=TwiceIntegrable;

  /***************************************************************/
  /* evaluate the 1-, 2-, or 3- dimensional cubature (for once-  */
  /* integrable kernels) or the 0-, 1-, or 2- dimensional        */
  /* cubature (for twice-integrable kernels).                    */
  /***************************************************************/
  static double Lower[3]={0.0, 0.0, 0.0};
  static double Upper[3]={1.0, 1.0, 1.0};
  int fDim=2*NumPKs;
  double *dResult=(double *)(Args->Result);
  double *dError=(double *)(Args->Error);
  TDW->nCalls=0;
  int IntegralDimension = 4 - WhichCase - TwiceIntegrable;

  if (IntegralDimension==0)
   TaylorDuffySum(0, 0, (void *)TDW, 0, dResult);
  else
   pcubature(fDim, TaylorDuffySum, (void *)TDW, IntegralDimension, 
             Lower, Upper, MaxEval, AbsTol, RelTol, 
             ERROR_INDIVIDUAL, dResult, dError);

  Args->nCalls = TDW->nCalls;

}

/***************************************************************/
/* convert a one-variable quadratic expression into a new form:*/
/*  Px^2 + 2Qx + R -> A^2 [ (x+B)^2 + G^2 ]                    */
/*                                                             */
/* note: the quantity returned as G2 is the quantity G^2.      */
/***************************************************************/
static void PQRtoABG2(double P, double Q, double R, double *A, double *B, double *G2)
{ 
  *A=sqrt(P);
  *B=Q/P;
  //*G2 = R/P - (*B)*(*B);

  double ROverP = R/P;
  *G2 = ROverP - (*B)*(*B);

  if ( fabs(*G2) < 1.0e-6*ROverP )
   { *B=0.0;
     *G2=ROverP;
   };
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GetAlphaBetaGamma2(TDWorkspace *TDW, const double *yVector,
                               double *AVector, double *BVector, double *G2Vector)
{
  if (TDW->WhichCase==TD_COMMONTRIANGLE) 
   { 
     double A2  = TDW->A2; 
     double AdB = TDW->AdB;
     double B2  = TDW->B2; 

     double C2 = A2 + B2 + 2.0*AdB;
     double CdB = AdB + B2;

     AVector[0] = sqrt(B2);
     AVector[1] = sqrt(C2);
     AVector[2] = sqrt(A2);

     BVector[0] = AdB/B2;
     BVector[1] = -CdB/C2;
     BVector[2] = AdB/A2;

     G2Vector[0] = A2/B2 - (AdB*AdB)/(B2*B2);
     G2Vector[1] = B2/C2 - (CdB*CdB)/(C2*C2);
     G2Vector[2] = B2/A2 - (AdB*AdB)/(A2*A2);
   }
  else if (TDW->WhichCase==TD_COMMONEDGE) 
   { 
     double y1 = yVector[0], y12=y1*y1;

     double A2   = TDW->A2;
     double BP2  = TDW->BP2;
     double L2   = TDW->L2;
     double AdBP = TDW->AdBP;
     double AdL  = TDW->AdL;
     double BPdL = TDW->BPdL;

     PQRtoABG2( (BP2 - 2.0*BPdL + L2)*y12, 
                y1*(L2 + BPdL*(y1-1.0) - (AdL + L2-AdBP)*y1),
                L2 - 2*AdL*y1 - 2*L2*y1 + A2*y12 + 2*AdL*y12 + L2*y12,
                AVector+0, BVector+0, G2Vector+0);

     PQRtoABG2( BP2*y12,
                (BPdL + AdBP*y1 - BPdL*y1)*y1,
                L2 + 2*AdL*y1 - 2*L2*y1 + (A2 - 2*AdL + L2)*y12,
                AVector+1, BVector+1, G2Vector+1);

     PQRtoABG2( (A2 + 2*AdBP + BP2)*y12, 
                y1*(AdL*(-1 + y1) + BPdL*(-1 + y1) - (AdBP + BP2)*y1),
                L2 + 2*BPdL*y1 - 2*L2*y1 + BP2*y12 - 2*BPdL*y12 + L2*y12,
                AVector+2, BVector+2, G2Vector+2);

     PQRtoABG2( (A2 + 2*AdBP - 2*AdL + BP2 - 2*BPdL + L2)*y12,
                y1*(AdL - L2 - (AdBP + BP2)*y1 + BPdL*(1 + y1)),
                L2 - 2*BPdL*y1 + BP2*y12,
                AVector+3, BVector+3, G2Vector+3);

     PQRtoABG2( A2*y12, 
                (AdBP*y1-AdL)*y1, 
                L2 - 2*BPdL*y1 + BP2*y12,
                AVector+4, BVector+4, G2Vector+4);

     PQRtoABG2( A2*y12,
                y1*(AdL + AdBP*y1 - AdL*y1),
                L2 + 2*BPdL*y1 - 2*L2*y1 + BP2*y12 - 2*BPdL*y12 + L2*y12,
                AVector+5, BVector+5, G2Vector+5);
   }
  else // (WhichCase==TD_COMMONVERTEX) 
   {
     double y1 = yVector[0];
     double y2 = yVector[1], y22=y2*y2;

     double A2    = TDW->A2;
     double B2    = TDW->B2;
     double AP2   = TDW->AP2;
     double BP2   = TDW->BP2;
     double AdB   = TDW->AdB;
     double AdAP  = TDW->AdAP;
     double AdBP  = TDW->AdBP;
     double BdAP  = TDW->BdAP;
     double BdBP  = TDW->BdBP;
     double APdBP = TDW->APdBP;

     PQRtoABG2( BP2*y22, 
                -AdBP*y2 -BdBP*y1*y2 +APdBP*y22,
                A2 + B2*y1*y1 + AP2*y2*y2 + 2.0*AdB*y1 - 2.0*AdAP*y2 -2.0*BdAP*y1*y2,
                AVector+0, BVector+0, G2Vector+0);

     PQRtoABG2( B2*y22,
                AdB*y22 - BdAP*y2 - BdBP*y1*y2,
                A2*y2*y2 + AP2 + BP2*y1*y1 - 2*AdAP*y2 - 2*AdBP*y2*y1 + 2*APdBP*y1,
                AVector+1, BVector+1, G2Vector+1);

   };

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GetX(TDWorkspace *TDW, const double *yVector, double *X)
{
  if (TDW->WhichCase==TD_COMMONTRIANGLE)
   { double y=yVector[0], u1, u2;
     u1=1.0; u2=y;     X[0]=sqrt( TDW->A2*u1*u1 + 2.0*TDW->AdB*u1*u2 + TDW->B2*u2*u2);
     u1=y;   u2=(y-1); X[1]=sqrt( TDW->A2*u1*u1 + 2.0*TDW->AdB*u1*u2 + TDW->B2*u2*u2);
     u1=y;   u2=1.0;   X[2]=sqrt( TDW->A2*u1*u1 + 2.0*TDW->AdB*u1*u2 + TDW->B2*u2*u2);
   }
  else if (TDW->WhichCase==TD_COMMONEDGE)
   { 
     double y1=yVector[0], y2=yVector[1], u1, u2, xi2;
     u1=-y1;    u2=-y1*y2;       xi2=(1.0-y1+y1*y2);
     X[0]=sqrt(       u1*u1*TDW->A2 + u2*u2*TDW->BP2 + xi2*xi2*TDW->L2 
                + 2.0*(u1*(u2*TDW->AdBP + xi2*TDW->AdL) + u2*xi2*TDW->BPdL) );

     u1=y1;     u2=y1*y2;        xi2=(1.0-y1);
     X[1]=sqrt(       u1*u1*TDW->A2 + u2*u2*TDW->BP2 + xi2*xi2*TDW->L2 
                + 2.0*(u1*(u2*TDW->AdBP + xi2*TDW->AdL) + u2*xi2*TDW->BPdL) );

     u1=-y1*y2; u2=y1*(1.0-y2);  xi2=(1.0-y1);
     X[2]=sqrt(       u1*u1*TDW->A2 + u2*u2*TDW->BP2 + xi2*xi2*TDW->L2 
                + 2.0*(u1*(u2*TDW->AdBP + xi2*TDW->AdL) + u2*xi2*TDW->BPdL) );

     u1=y1*y2;  u2=-y1*(1.0-y2); xi2=(1.0-y1*y2);
     X[3]=sqrt(       u1*u1*TDW->A2 + u2*u2*TDW->BP2 + xi2*xi2*TDW->L2 
                + 2.0*(u1*(u2*TDW->AdBP + xi2*TDW->AdL) + u2*xi2*TDW->BPdL) );

     u1=-y1*y2; u2=-y1;          xi2=1.0;
     X[4]=sqrt(       u1*u1*TDW->A2 + u2*u2*TDW->BP2 + xi2*xi2*TDW->L2 
                + 2.0*(u1*(u2*TDW->AdBP + xi2*TDW->AdL) + u2*xi2*TDW->BPdL) );

     u1=y1*y2;  u2=y1;           xi2=(1.0-y1);
     X[5]=sqrt(       u1*u1*TDW->A2 + u2*u2*TDW->BP2 + xi2*xi2*TDW->L2 
                + 2.0*(u1*(u2*TDW->AdBP + xi2*TDW->AdL) + u2*xi2*TDW->BPdL) );

   }
  else // (TDW->WhichCase==TD_COMMONVERTEX)
   { double y1=yVector[0], y2=yVector[1], y3=yVector[2], xi1, xi2, eta1, eta2;

     xi1=1.0; xi2=y1; eta1=y2; eta2=y2*y3;
     X[0] = sqrt(     xi1*xi1*TDW->A2    +   xi2*xi2*TDW->B2 
                    + eta1*eta1*TDW->AP2   + eta2*eta2*TDW->BP2
                    + 2.0*xi1*(xi2*TDW->AdB - eta1*TDW->AdAP - eta2*TDW->AdBP)
                    - 2.0*xi2*(eta1*TDW->BdAP + eta2*TDW->BdBP) 
                    + 2.0*eta1*eta2*TDW->APdBP
                 );

     xi1=y2; xi2=y2*y3; eta1=1.0; eta2=y1;
     X[1] = sqrt(     xi1*xi1*TDW->A2    +   xi2*xi2*TDW->B2 
                    + eta1*eta1*TDW->AP2   + eta2*eta2*TDW->BP2
                    + 2.0*xi1*(xi2*TDW->AdB - eta1*TDW->AdAP - eta2*TDW->AdBP)
                    - 2.0*xi2*(eta1*TDW->BdAP + eta2*TDW->BdBP) 
                    + 2.0*eta1*eta2*TDW->APdBP
                 );

   };
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GetScriptP(TDWorkspace *TDW, int WhichP, const double *yVector, 
                double P[NUMREGIONS][NUMWPOWERS][NUMYPOWERS])
{
  int nMin = TDW->nMin[WhichP]; 
  int nMax = TDW->nMax[WhichP]; 

  if ( TDW->WhichCase == TD_COMMONTRIANGLE)
   { 
     for(int d=0; d<3; d++)
      for(int n=nMin; n<=nMax; n++)
       { P[d][n][0] = TDW->Upsilon[WhichP][d][n][0];
         P[d][n][1] = TDW->Upsilon[WhichP][d][n][Y1TERM];
         P[d][n][2] = TDW->Upsilon[WhichP][d][n][Y12TERM];
       };
   }
  else if ( TDW->WhichCase == TD_COMMONEDGE )
   { 
     double y1=yVector[0], y12=y1*y1;
     for(int d=0; d<6; d++)
      for(int n=nMin; n<=nMax; n++)
       { 
         P[d][n][0]=  TDW->Upsilon[WhichP][d][n][0]
                    + TDW->Upsilon[WhichP][d][n][Y1TERM]     * y1
                    + TDW->Upsilon[WhichP][d][n][Y12TERM]    * y12;

         P[d][n][1]=  TDW->Upsilon[WhichP][d][n][Y2TERM]     
                    + TDW->Upsilon[WhichP][d][n][Y1Y2TERM]   * y1
                    + TDW->Upsilon[WhichP][d][n][Y12Y2TERM]  * y12;

         P[d][n][2]=  TDW->Upsilon[WhichP][d][n][Y22TERM]
                    + TDW->Upsilon[WhichP][d][n][Y1Y22TERM]  * y1
                    + TDW->Upsilon[WhichP][d][n][Y12Y22TERM] * y12;
       };
   }
  else // ( TDW->WhichCase == TD_COMMONVERTEX )
   { 
     double y1=yVector[0], y12=y1*y1;
     double y2=yVector[1], y22=y2*y2;

     for(int d=0; d<2; d++)
      for(int n=nMin; n<=nMax; n++)
       { 
         P[d][n][0] = TDW->Upsilon[WhichP][d][n][0]
                     +TDW->Upsilon[WhichP][d][n][Y1TERM]    * y1
                     +TDW->Upsilon[WhichP][d][n][Y2TERM]    * y2
                     +TDW->Upsilon[WhichP][d][n][Y12TERM]   * y12
                     +TDW->Upsilon[WhichP][d][n][Y1Y2TERM]  * y1*y2
                     +TDW->Upsilon[WhichP][d][n][Y22TERM]   * y22
                     +TDW->Upsilon[WhichP][d][n][Y12Y2TERM] * y12*y2
                     +TDW->Upsilon[WhichP][d][n][Y1Y22TERM] * y1*y22;

         P[d][n][1] = TDW->Upsilon[WhichP][d][n][Y2Y3TERM]    * y2
                     +TDW->Upsilon[WhichP][d][n][Y1Y2Y3TERM]  * y1*y2
                     +TDW->Upsilon[WhichP][d][n][Y22Y3TERM]   * y22
                     +TDW->Upsilon[WhichP][d][n][Y12Y2Y3TERM] * y12*y2;

         P[d][n][2] = TDW->Upsilon[WhichP][d][n][Y22Y32TERM]   * y22
                     +TDW->Upsilon[WhichP][d][n][Y1Y22Y32TERM] * y1*y22;
       };
   }
  
}

/***************************************************************/
/* This routine computes the two quantities                    */
/*                                                             */
/* IntQFP(P,Q2)  = \int_0^1 [(y+P)^2 + Q2]^{p/2} dy            */
/* IntyQFP(P,Q2) = \int_0^1 y[(y+P)^2 + Q2]^{p/2} dy           */
/*                                                             */
/* Note: QFP stands for 'quadratic form to the power P.'       */
/*                                                             */
/* (slightly confusing: the exponent is p/2, not p. i should   */
/*  call the parameter 'p' something like 'TwoP' but that's    */
/*  unwieldy.)                                                 */
/***************************************************************/
void GetQFPIntegral(double P, double Q2, int p, double *IntQFP, double *IntyQFP)
{ 
  double P2     = P*P; 
  double PP1    = P+1.0, PP12=PP1*PP1;
  double Q4     = Q2*Q2;
  double S2     = P*P + Q2, S=sqrt(S2); 
  double T2     = PP1*PP1 + Q2, T=sqrt(T2); 

  double Q, Q6, P4, PP14, T3, T4, T5, T6, T7, T8, S3, S4, S5, S6, S7, S8;
  double LogFac, atanFac;
  double SignP = (P>0.0) ? 1.0 : -1.0;

  switch(p)
   { 
      case  0: IntQFP[0] = 1.0;
               IntyQFP[0] = 0.5;
               return;

      case  1: S3=S2*S; T3=T2*T; 
               if ( Q2 < 1.0e-8 )
                IntQFP[0] = SignP*(T2-S2)/2.0;
               else 
                { LogFac = log( (PP1+T) / (P+S) );
                  IntQFP[0]  = (T+P*(T-S) + Q2*LogFac) / 2.0;
                };
               IntyQFP[0] = -P*IntQFP[0] + (T3-S3)/3.0;
               return;

      case -1: if ( Q2<1.0e-8 )
                IntQFP[0] = SignP*log( 1.0 + 1.0/P );
               else
                IntQFP[0] = log( (PP1+T) / (P+S) );
               IntyQFP[0] = -P*IntQFP[0] + (T-S);
               return;

      case  2: S4=S2*S2; T4=T2*T2; 
               IntQFP[0]  = (T2+S2)/2.0 - 1.0/6.0;
               IntyQFP[0] = -P*IntQFP[0] + (T4-S4)/4.0;
               return;

      case -2: if ( Q2 < 1.0e-8 ) 
                IntQFP[0] = 1.0/(P*PP1);
               else 
                { Q=sqrt(Q2);
                  IntQFP[0]  = ( atan(PP1/Q) - atan(P/Q) ) / Q;
                };
               IntyQFP[0] = -P*IntQFP[0] + log(T/S);
               return;

      case  3: S4=S2*S2; S5=S4*S; 
               T4=T2*T2; T5=T4*T; 
               if ( Q2 < 1.0e-8 )
                IntQFP[0] = SignP*(T4-S4)/4.0;
               else
                { LogFac = log( (PP1+T) / (P+S) );
                  IntQFP[0]  = ( PP1*(2.0*PP12+5.0*Q2)*T 
                                 - P*(2.0*P2  +5.0*Q2)*S 
                                 + 3.0*Q4*LogFac)/8.0;
                };
               IntyQFP[0] = -P*IntQFP[0]  + (T5-S5) / 5.0;
               return;

      case -3: if ( Q2<1.0e-8) 
                IntQFP[0]  = SignP*( 1.0/(2.0*S2) - 1.0/(2.0*T2));
               else
                IntQFP[0]  = ( PP1/T - P/S ) / Q2;
               IntyQFP[0] = -P*IntQFP[0] + 1.0/S - 1.0/T;
               return;

      case  4: T4=T2*T2; T6=T4*T2;
               S4=S2*S2; S6=S4*S2;
               IntQFP[0]  = (T4+S4+T2*S2)/3.0 - (T2+S2)/6.0 + 1.0/30.0;
               IntyQFP[0]  = -P*IntQFP[0]  + (T6-S6) / 6.0;
               return;

      case -4: if ( Q2<1.0e-6 )
                { S3=S*S2; T3=T*T2;
                  IntQFP[0] = SignP*(1.0/(3.0*S3) - 1.0/(3.0*T3));
                }
               else
                { Q=sqrt(Q2);
                  atanFac=( atan(PP1/Q) - atan(P/Q) ) / Q;
                  IntQFP[0]  = ((Q2-P*PP1)/(S2*T2) + atanFac) / (2.0*Q2);
                }
               IntyQFP[0] = -P*IntQFP[0] + 1.0/(2.0*S2) - 1.0/(2.0*T2);
               return;

      case  5: 
               P4=P2*P2; PP14=PP12*PP12; Q6=Q4*Q2;
               S6=S2*S2*S2; S7=S6*S; 
               T6=T2*T2*T2; T7=T6*T; 
               if ( Q2<1.0e-8 )
                IntQFP[0] = SignP*(T6-S6)/6.0;
               else
                { 
                  LogFac = log( (PP1+T) / (P+S) );
                  IntQFP[0]  = (PP1*(8.0*PP14+26.0*PP12*Q2+33.0*Q4)*T 
                                 -P*(8.0*P4  +26.0*P2*Q2  +33.0*Q4)*S 
                                 + 15.0*Q6*LogFac) / 48.0;
                };
               IntyQFP[0] = -P*IntQFP[0]  + (T7-S7) / 7.0;
               return;

      case -5: if ( Q2<1.0e-4 ) // note: lower thresholds for Q2 for p==-5, -6
                { S3=S2*S;  S4=S3*S;
                  T3=T2*T;  T4=T3*T;
                  IntQFP[0] = SignP*(1.0/(4.0*S4) - 1.0/(4.0*T4));
                }
               else
                { S3 = S*S2; T3 = T*T2;
                  IntQFP[0] = (PP1*(2.0*PP12 + 3.0*Q2)*S3 - 
                                 P*(2.0*P2   + 3.0*Q2)*T3)/(3.0*Q4*S3*T3);
                };
               IntyQFP[0] = -P*IntQFP[0] + 1.0/(3.0*S3) - 1.0/(3.0*T3);
               return;

      case 6:  T4 = T2*T2; T6 = T4*T2; T8=T6*T2; 
               S4 = S2*S2; S6 = S4*S2; S8=S6*S2;
               IntQFP[0]  =  (T6 + S6 + S2*T4 + T2*S4)/4.0 
                           - (3.0*T4+3.0*S4+4.0*S2*T2)/20.0
                           + (T2+S2)/20.0 - 1.0/140.0;
               IntyQFP[0]  = -P*IntQFP[0]  + (T8-S8) / 8.0;
               return;

      case -6: if ( Q2<1.0e-4 )
                { S4=S2*S2;  S5=S4*S;
                  T4=T2*T2;  T5=T4*T;
                  IntQFP[0] = SignP*(1.0/(5.0*S5) - 1.0/(5.0*T5));
                }
               else
                { Q=sqrt(Q2);
                  Q4=Q2*Q2; 
                  T4 = T2*T2; T6 = T4*T2; T8=T6*T2; 
                  S4 = S2*S2; S6 = S4*S2; S8=S6*S2;
                  atanFac=( atan(PP1/Q) - atan(P/Q) ) / Q;
                  IntQFP[0] = (  3.0*PP1/T2 + 2.0*Q2*PP1/T4  
                                -3.0*P/S2   - 2.0*Q2*P/S4 + 3*atanFac)/(8.0*Q4);
                };
               IntyQFP[0] = -P*IntQFP[0] + 1.0/(4.0*S4) - 1.0/(4.0*T4);
               return;

     default:  ErrExit("%s:%i: internal error (%i)",__FILE__,__LINE__,p);

   }; // switch(p)

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
static double FactorialTable[7]={1.0, 1.0, 2.0, 6.0, 24.0, 120.0, 620.0};
void GetScriptJL(int WhichK, cdouble KParam,
                 double Alpha, double Beta, double Gamma2,
                 int nMin, int nMax,
                 cdouble JVector[7], cdouble LVector[7])
{ 
  double IntQFP, IntyQFP;

  if (WhichK==TD_RP)
   { 
     double p=real(KParam);
     double AlphaP=pow(Alpha, p);
     if (p==0) 
      { IntQFP=1.0; IntyQFP=0.5; }
     else
      GetQFPIntegral(Beta, Gamma2, (int) p, &IntQFP, &IntyQFP);
     for(int n=nMin; n<=nMax; n++)
      { JVector[n] = AlphaP*IntQFP  / (1.0 + n + p);
        LVector[n] = AlphaP*IntyQFP / (1.0 + n + p);
      };
   }
  else if (WhichK==TD_HIGHK_HELMHOLTZ)
   { 
     for(int n=nMin; n<=nMax; n++)
      { cdouble ikAlphaMN=pow(-II*KParam*Alpha, -(double)n);
        GetQFPIntegral(Beta, Gamma2, -(n+1), &IntQFP, &IntyQFP);
        JVector[n] = FactorialTable[n-1]*ikAlphaMN*IntQFP/Alpha;
        LVector[n] = FactorialTable[n-1]*ikAlphaMN*IntyQFP/Alpha;
      };
   }
  else if (WhichK==TD_HIGHK_GRADHELMHOLTZ)
   { double Alpha3=Alpha*Alpha*Alpha;
     for(int n=nMin; n<=nMax; n++)
      { cdouble ikAlphaMNM2=pow(-II*KParam*Alpha, -(double)(n-2));
        GetQFPIntegral(Beta, Gamma2, -(n+1), &IntQFP, &IntyQFP);
        JVector[n] = -(n-1)*FactorialTable[n-3]*ikAlphaMNM2*IntQFP/Alpha3;
        LVector[n] = -(n-1)*FactorialTable[n-3]*ikAlphaMNM2*IntyQFP/Alpha3;
      };
   };

}


/***************************************************************/
/* this is what the GSL calls the 'relative exponential'       */
/* function. it is equal to exp(Z), minus the first n terms in */
/* the taylor expansion for exp(Z), divided by the (n+1)th     */
/* term in the taylor expansion for exp(Z). (this latter       */
/* normalization ensures that ExpRel(n,0) = 1 for all n.)      */
/***************************************************************/
#define EXPRELTOL  1.0e-8
#define EXPRELTOL2 EXPRELTOL*EXPRELTOL
cdouble ExpRelV3P0(int n, cdouble Z)
{
  int m;

  /*--------------------------------------------------------------*/
  /*- purely real case                                           -*/
  /*--------------------------------------------------------------*/
  if ( imag(Z)==0.0 )
   { 
     double Term, Sum;
     double X=real(Z);

     //////////////////////////////////////////////////
     // small-Z expansion
     //////////////////////////////////////////////////
     if ( fabs(X) < 0.1 )
      { Sum=1.0;
        for(Term=1.0, m=1; m<100; m++)
         { Term*=X/((double)(m+n));
           Sum+=Term;
           if ( fabs(Term) < EXPRELTOL*fabs(Sum) )
            break;
         };
        return Sum;
      }
     else
      { Sum=exp(X);
        for(Term=1.0, m=0; m<n; m++)
         { 
          Sum-=Term;
           Term*=X/((double)(m+1));
         };
        return Sum / Term;
      };
   }
  else
   { 
     cdouble Term, Sum;
     /*--------------------------------------------------------------*/
     /*- small-Z expansion                                          -*/
     /*--------------------------------------------------------------*/
     if ( abs(Z) < 0.1 )
      { for(Sum=Term=1.0, m=1; m<100; m++)
         { Term*=Z/((double)(m+n));
           Sum+=Term;
           if ( norm(Term) < EXPRELTOL2*norm(Sum) )
            break;
         };
        return Sum;
      }
     else
      {
        Sum=exp(Z);
        for(Term=1.0, m=0; m<n; m++)
         { 
           Sum-=Term; 
           Term*=Z/((double)(m+1));
         };
        return Sum / Term;
      };
   };

} 

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GetScriptK(int WhichK, cdouble KParam, double X,
                       int nMin, int nMax, cdouble KVector[7])
{
  if (WhichK==TD_RP)
   { 
     double P=real(KParam);
     double XP=pow(X,P);
     for(int n=nMin; n<=nMax; n++)
      KVector[n] = XP / ( 1.0 + P + (double)n );
   }
  else if (WhichK==TD_HELMHOLTZ)
   { cdouble IK = II*KParam, IKX = IK*X, eIKX = exp(IKX);
     for(int n=nMin; n<=nMax; n++)
      KVector[n] = eIKX * ExpRelV3P0(n,-IKX) / (n*X);
   }
  else if (WhichK==TD_GRADHELMHOLTZ)
   { cdouble IK = II*KParam, IKX = IK*X, eIKX = exp(IKX);
     cdouble ExpRelTable[10]; 

     if ( (nMin-2) < 0 || (nMax-1) >=10)
      { Warn("%s:%i: internal inconsistency (%i,%i)",__FILE__,__LINE__,nMin,nMax);
        for(int n=0; n<nMin; n++)
         KVector[n]=0.0;
        nMin=2;
      }; 

      //ErrExit("%s:%i: internal error",__FILE__,__LINE__);
 
     // precompute ExpRel for all values we need; for better efficiency,
     // maybe implement a routine that computes these values all at once?
     for(int n=nMin-2; n<=nMax-1; n++)
      ExpRelTable[n] = ExpRelV3P0(n,-IKX);

     for(int n=nMin; n<=nMax; n++)
      KVector[n] = eIKX * ( IK*ExpRelTable[n-1] / ((double)n-1.0)
                              -ExpRelTable[n-2] / (((double)n-2.0)*X)
                          ) / (X*X);
   };
  
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void CMVStoUpsilon(int WhichCase,
                   double *C, double *M, double *V, double S,
                   double Upsilon[NUMREGIONS][NUMWPOWERS][NUMMONOMIALS])
{
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  for(int nr=0; nr<NUMREGIONS; nr++)
   for(int nwp=0; nwp<NUMWPOWERS; nwp++)
    memset(Upsilon[nr][nwp],0,NUMMONOMIALS*sizeof(double));

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  double C2010 = C[0], C2001 = C[1], C0210 = C[2], C0201 = C[3];
  double M11 = M[0 + 4*0], M12 = M[0 + 4*1], M13 = M[0 + 4*2], M14 = M[0 + 4*3];
  double                   M22 = M[1 + 4*1], M23 = M[1 + 4*2], M24 = M[1 + 4*3];
  double                                     M33 = M[2 + 4*2], M34 = M[2 + 4*3];
  double                                                       M44 = M[3 + 4*3];
  double V1=V[0], V2=V[1], V3=V[2], V4=V[3];

  if (WhichCase == TD_COMMONTRIANGLE)
   { 
      Upsilon[0][0][0]       =  M11/4 + M12/4 + M13/2 + M14/4 + M22/12 + M23/4 + M24/6 + M33/4
                               +M34/4 + M44/12 +  (2*V1)/3 + V2/3 + (2*V3)/3 + V4/3 + S;

      Upsilon[0][1][Y1TERM]  =  M12/3 + M14/3 + M22/6 + M23/3 + M24/3 + M34/3 + M44/6 + V2/2 + V4/2;
      Upsilon[0][1][0]       = -(2*M11)/3 - (5*M12)/6 - (4*M13)/3 - (5*M14)/6 - M22/3 - (5*M23)/6
                               -(2*M24)/3 - (2*M33)/3 - (5*M34)/6 - M44/3 - (3*V1)/2 - V2 -(3*V3)/2 
                               - V4 - 2*S;

      Upsilon[0][2][Y12TERM] =  M22/4 + M44/4;
      Upsilon[0][2][Y1TERM]  = -M12/2 - M14 - M22/2 - M23 - M24 - M34/2 - M44/2 - V2 - V4;
      Upsilon[0][2][0]       =  (3*M11)/4 + M12 + M13 + M14 + M22/2 + M23 + M24 
                               +(3*M33)/4 + M34 + M44/2 + V1 + V2 + V3 + V4 + S;

      Upsilon[0][3][Y12TERM] = -M22/2 - M44/2; 
      Upsilon[0][3][Y1TERM]  =  M14 + M22/2 + M23 + M24 + M44/2 + V2/2 + V4/2;
      Upsilon[0][3][0]       = -M11/2 - M12/2 - M14/2 - M22/3 - M23/2 - (2*M24)/3 
                               -M33/2 - M34/2 - M44/3 - V1/6 - V2/3 - V3/6 - V4/3;

      Upsilon[0][4][Y12TERM] =  M22/4 + M44/4;
      Upsilon[0][4][Y1TERM]  =  M12/6 - M14/3 - M22/6 - M23/3 - M24/3 + M34/6 - M44/6;
      Upsilon[0][4][0]       =  M11/6 + M12/12 - M13/6 + M14/12 + M22/12 + M23/12 + M24/6 
                               +M33/6 + M34/12 + M44/12;

      Upsilon[1][0][0]       =  M11/4 + M12/4 + M13/2 + M14/4 + M22/12 + M23/4 + M24/6 + M33/4
                               +M34/4 + M44/12 +  (2*V1)/3 + V2/3 + (2*V3)/3 + V4/3 + S;

      Upsilon[1][1][Y1TERM]  = -M11/3 - M12/2 - (2*M13)/3 - M14/2 - M22/6 - M23/2 - M24/3 - M33/3
                               -M34/2 - M44/6 - V1/2 - V2/2 - V3/2 - V4/2;
      Upsilon[1][1][0]       = -M11/3 - M12/3 - (2*M13)/3 - M14/3 - M22/6 - M23/3 - M24/3 - M33/3 
                               -M34/3 - M44/6 - 2*S - V1 - V2/2 - V3 - V4/2;

      Upsilon[1][2][Y12TERM] =  M11/4 + M12/2 + M22/4 + M33/4 + M34/2 + M44/4;
      Upsilon[1][2][Y1TERM]  =  M11/2 + M12/2 + M13 + M14 + M23 + M24 + M33/2 + M34/2 
                               +V1 + V2 + V3 + V4;
      Upsilon[1][2][0]       =  M22/4 + M44/4 + S;

      Upsilon[1][3][Y12TERM] = -M11/2 - M12 - M22/2 - M33/2 - M34 - M44/2;
      Upsilon[1][3][Y1TERM]  =  M12/2 - M14/2 + M22/2 - M23/2 - M24 + M34/2 + M44/2
                               -V1/2 - V2/2 - V3/2 - V4/2;
      Upsilon[1][3][0]       = -M22/3 + M24/3 - M44/3 + V1/3 + V2/6 + V3/3 + V4/6;

      Upsilon[1][4][Y12TERM] =  M11/4 + M12/2 + M22/4 + M33/4 + M34/2 + M44/4;
      Upsilon[1][4][Y1TERM]  = -M11/6 - M12/2 - M13/3 - M22/3 + M24/3 - M33/6 - M34/2 - M44/3;
      Upsilon[1][4][0]       =  M11/12 + M12/12 + M13/6 + M14/12 + M22/6 + M23/12 - M24/6 
                               +M33/12 + M34/12 + M44/6;

      Upsilon[2][0][0]       =  M11/4 + M12/4 + M13/2 + M14/4 + M22/12 + M23/4 + M24/6 + M33/4
                               +M34/4 + M44/12 + S + (2*V1)/3 + V2/3 + (2*V3)/3 + V4/3;

      Upsilon[2][1][Y1TERM]  = -M11/3 - M12/6 - (2*M13)/3 - M14/6 - M23/6 - M33/3 - M34/6 - V1/2 - V3/2;
      Upsilon[2][1][0]       = -M11/3 - M12/3 - (2*M13)/3 - M14/3 - M22/6 - M23/3 - M24/3 
                               -M33/3 - M34/3 - M44/6 - 2*S - V1 - V2/2 - V3 - V4/2;

      Upsilon[2][2][Y12TERM] =  M11/4 + M33/4;
      Upsilon[2][2][Y1TERM]  =  M11/2 + M12/2 + M13 + M33/2 + M34/2 + V1 + V3;
      Upsilon[2][2][0]       =  M22/4 + M44/4 + S;

      Upsilon[2][3][Y12TERM] = -M11/2 - M33/2; 
      Upsilon[2][3][Y1TERM]  = -M12/2 + M14/2 + M23/2 - M34/2 - V1/2 - V3/2;
      Upsilon[2][3][0]       = -M22/3 + M24/3 - M44/3 + V1/3 + V2/6 + V3/3 + V4/6;

      Upsilon[2][4][Y12TERM] =  M11/4 + M33/4;
      Upsilon[2][4][Y1TERM]  = -M11/6 + M12/6 - M13/3 - M14/3 - M23/3 - M33/6 + M34/6;
      Upsilon[2][4][0]       =  M11/12 + M12/12 + M13/6 + M14/12 + M22/6 + M23/12 - M24/6 
                               +M33/12 + M34/12 + M44/6;

   }
  else if (WhichCase == TD_COMMONEDGE)
   {
      Upsilon[0][0][0]          = C2010/4 + M11/6 + M13/3 + M33/6 + S + V1/2 + V3/2; 
     
      Upsilon[0][1][Y1Y2TERM]   = M12/2 + M23/2 + V2;
      Upsilon[0][1][Y1TERM]     = -C2001/3 - C2010/3 - M12/2 - M13/2 - M14/2 - M23/2 - M33/2 - M34/2 - V2 - V3 - V4;
      Upsilon[0][1][0]          = C2001/3 + M12/2 + M14/2 + M23/2 + M34/2 - S + V2 + V4;
    
      Upsilon[0][2][Y12Y22TERM] = C0210/2 + M22/2;
      Upsilon[0][2][Y12Y2TERM]  = -C0210 - M22 - M23 - M24;
      Upsilon[0][2][Y12TERM]    = C0210/2 + M22/2 + M23 + M24 + M33/2 + M34 + M44/2;
      Upsilon[0][2][Y1Y2TERM]   = C0210 + M22 + M24 - V2;
      Upsilon[0][2][Y1TERM]     = -C0210 - M22 - M23 - 2*M24 - M34 - M44 + V2 + V3 + V4;
      Upsilon[0][2][0]          = C0210/2 + M22/2 + M24 + M44/2 - V1/2 - V2 - V3/2 - V4;
    
      Upsilon[0][3][Y13Y22TERM] = -C0201 - C0210;
      Upsilon[0][3][Y13Y2TERM]  = 2*C0201 + 2*C0210;
      Upsilon[0][3][Y13TERM]    = -C0201 - C0210;
      Upsilon[0][3][Y12Y22TERM] = C0201 - M22/2;
      Upsilon[0][3][Y12Y2TERM]  = -4*C0201 - 2*C0210 + M22 + M23 + M24;
      Upsilon[0][3][Y12TERM]    = 3*C0201 + 2*C0210 - M22/2 - M23 - M24 - M33/2 - M34 - M44/2;
      Upsilon[0][3][Y1Y2TERM]   = 2*C0201 - M12/2 - M22 - M23/2 - M24;
      Upsilon[0][3][Y1TERM]     = -3*C0201 - C0210 + M12/2 + M13/2 + M14/2 + M22 + (3*M23)/2 + 2*M24 + M33/2 + (3*M34)/2 + M44;
      Upsilon[0][3][0]          = C0201 - M11/6 - M12/2 - M13/3 - M14/2 - M22/2 - M23/2 - M24 - M33/6 - M34/2 - M44/2;
    
      Upsilon[0][4][Y13Y22TERM] = C0201 + C0210;
      Upsilon[0][4][Y13Y2TERM]  = -2*C0201 - 2*C0210;
      Upsilon[0][4][Y13TERM]    = C0201 + C0210;
      Upsilon[0][4][Y12Y22TERM] = -C0201 - C0210/2;
      Upsilon[0][4][Y12Y2TERM]  = 4*C0201 + 3*C0210;
      Upsilon[0][4][Y12TERM]    = -3*C0201 - (5*C0210)/2;
      Upsilon[0][4][Y1Y2TERM]   = -2*C0201 - C0210;
      Upsilon[0][4][Y1TERM]     = 3*C0201 + 2*C0210 + C2001/3 + C2010/3;
      Upsilon[0][4][0]          = -C0201 - C0210/2 - C2001/3 - C2010/4;

      Upsilon[1][0][0]          = C2010/4 + M11/6 + M13/3 + M33/6 + S + V1/2 + V3/2;
    
      Upsilon[1][1][Y1Y2TERM]   = C2001/3 + M14/2 + M34/2 + V4;
      Upsilon[1][1][Y1TERM]     = -C2001/3 - (2*C2010)/3 - M11/2 - M12/2 - M13/2 - M14/2 - M23/2 - M34/2 - V1 - V2 - V4;
      Upsilon[1][1][0]          = C2001/3 + M12/2 + M14/2 + M23/2 + M34/2 - S + V2 + V4;
    
      Upsilon[1][2][Y12Y22TERM] = M44/2;
      Upsilon[1][2][Y12Y2TERM]  = -C2001 - M14 - M24 - M44;
      Upsilon[1][2][Y12TERM]    = C0210/2 + C2001 + C2010/2 + M11/2 + M12 + M14 + M22/2 + M24 + M44/2;
      Upsilon[1][2][Y1Y2TERM]   = M24 + M44 - V4;
      Upsilon[1][2][Y1TERM]     = -C0210 - C2001 - M12 - M14 - M22 - 2*M24 - M44 + V1 + V2 + V4;
      Upsilon[1][2][0]          = C0210/2 + M22/2 + M24 + M44/2 - V1/2 - V2 - V3/2 - V4;
    
      Upsilon[1][3][Y13Y2TERM]  = C0201 + C2001;
      Upsilon[1][3][Y13TERM]    = -C0201 - C2001;
      Upsilon[1][3][Y12Y22TERM] = -M44/2;
      Upsilon[1][3][Y12Y2TERM]  = -2*C0201 + M14 + M24 + M44;
      Upsilon[1][3][Y12TERM]    = 3*C0201 + C2001 - M11/2 - M12 - M14 - M22/2 - M24 - M44/2;
      Upsilon[1][3][Y1Y2TERM]   = C0201 - M14/2 - M24 - M34/2 - M44;
      Upsilon[1][3][Y1TERM]     = -3*C0201 + M11/2 + (3*M12)/2 + M13/2 + (3*M14)/2 + M22 + M23/2 + 2*M24 + M34/2 + M44;
      Upsilon[1][3][0]          = C0201 - M11/6 - M12/2 - M13/3 - M14/2 - M22/2 - M23/2 - M24 - M33/6 - M34/2 - M44/2;
    
      Upsilon[1][4][Y13Y2TERM]  = -C0201 - C2001;
      Upsilon[1][4][Y13TERM]    = C0201 + C2001;
      Upsilon[1][4][Y12Y2TERM]  = 2*C0201 + C2001;
      Upsilon[1][4][Y12TERM]    = -3*C0201 - C0210/2 - 2*C2001 - C2010/2;
      Upsilon[1][4][Y1Y2TERM]   = -C0201 - C2001/3;
      Upsilon[1][4][Y1TERM]     = 3*C0201 + C0210 + (4*C2001)/3 + (2*C2010)/3;
      Upsilon[1][4][0]          = -C0201 - C0210/2 - C2001/3 - C2010/4;

      Upsilon[2][0][0]          = C2010/4 + M11/6 + M13/3 + M33/6 + S + V1/2 + V3/2;
     
      Upsilon[2][1][Y1Y2TERM]   = -C2001/3 - C2010/3 - M13/2 - M14/2 - M33/2 - M34/2 - V3 - V4;
      Upsilon[2][1][Y1TERM]     = -M12/2 - M23/2 - V2;
      Upsilon[2][1][0]          = C2001/3 + M12/2 + M14/2 + M23/2 + M34/2 - S + V2 + V4;
     
      Upsilon[2][2][Y12Y22TERM] = M33/2 + M34 + M44/2;
      Upsilon[2][2][Y12Y2TERM]  = M23 + M24;
      Upsilon[2][2][Y12TERM]    = C0210/2 + M22/2;
      Upsilon[2][2][Y1Y2TERM]   = -M23 - M24 - M34 - M44 + V3 + V4;
      Upsilon[2][2][Y1TERM]     = -C0210 - M22 - M24 + V2;
      Upsilon[2][2][0]          = C0210/2 + M22/2 + M24 + M44/2 - V1/2 - V2 - V3/2 - V4;
     
      Upsilon[2][3][Y13Y2TERM]  = -C0201 - C0210;
      Upsilon[2][3][Y12Y22TERM] = -M33/2 - M34 - M44/2;
      Upsilon[2][3][Y12Y2TERM]  = 2*C0201 + 2*C0210 - M23 - M24;
      Upsilon[2][3][Y12TERM]    = C0201 - M22/2;
      Upsilon[2][3][Y1Y2TERM]   = -C0201 - C0210 + M13/2 + M14/2 + M23 + M24 + M33/2 + (3*M34)/2
                                  +M44;
      Upsilon[2][3][Y1TERM]     = -2*C0201 + M12/2 + M22 + M23/2 + M24;
      Upsilon[2][3][0]          = C0201 - M11/6 - M12/2 - M13/3 - M14/2 - M22/2 - M23/2 - M24
                                  -M33/6 - M34/2 - M44/2;
    
      Upsilon[2][4][Y13Y2TERM]  = C0201 + C0210;
      Upsilon[2][4][Y12Y2TERM]  = -2*C0201 - 2*C0210;
      Upsilon[2][4][Y12TERM]    = -C0201 - C0210/2;
      Upsilon[2][4][Y1Y2TERM]   = C0201 + C0210 + C2001/3 + C2010/3;
      Upsilon[2][4][Y1TERM]     = 2*C0201 + C0210;
      Upsilon[2][4][0]          = -C0201 - C0210/2 - C2001/3 - C2010/4;



      Upsilon[3][0][0]          = C2010/4 + M11/6 + M13/3 + M33/6 + S + V1/2 + V3/2;
     
      Upsilon[3][1][Y1Y2TERM]   = (-2*C2010)/3 - M11/2 - M12/2 - M13/2 - M23/2 - V1 - V2;
      Upsilon[3][1][Y1TERM]     = -C2001/3 - M14/2 - M34/2 - V4;
      Upsilon[3][1][0]          = C2001/3 + M12/2 + M14/2 + M23/2 + M34/2 - S + V2 + V4;
     
      Upsilon[3][2][Y12Y22TERM] = C0210/2 + C2010/2 + M11/2 + M12 + M22/2;
      Upsilon[3][2][Y12Y2TERM]  = C2001 + M14 + M24;
      Upsilon[3][2][Y12TERM]    = M44/2;
      Upsilon[3][2][Y1Y2TERM]   = -C0210 - C2001 - M12 - M14 - M22 - M24 + V1 + V2;
      Upsilon[3][2][Y1TERM]     = -M24 - M44 + V4;
      Upsilon[3][2][0]          = C0210/2 + M22/2 + M24 + M44/2 - V1/2 - V2 - V3/2 - V4;
    
      Upsilon[3][3][Y13Y22TERM] = -C0201 - C2001;
      Upsilon[3][3][Y12Y22TERM] = C0201 + C2001 - M11/2 - M12 - M22/2;
      Upsilon[3][3][Y12Y2TERM]  = 2*C0201 - M14 - M24;
      Upsilon[3][3][Y12TERM]    = -M44/2;
      Upsilon[3][3][Y1Y2TERM]   = -2*C0201 + M11/2 + (3*M12)/2 + M13/2 + M14 + M22 + M23/2 + M24;
      Upsilon[3][3][Y1TERM]     = -C0201 + M14/2 + M24 + M34/2 + M44;
      Upsilon[3][3][0]          = C0201 - M11/6 - M12/2 - M13/3 - M14/2 - M22/2 - M23/2 - M24
                                   -M33/6 - M34/2 - M44/2;
    
      Upsilon[3][4][Y13Y22TERM] = C0201 + C2001;
      Upsilon[3][4][Y12Y22TERM] = -C0201 - C0210/2 - C2001 - C2010/2;
      Upsilon[3][4][Y12Y2TERM]  = -2*C0201 - C2001;
      Upsilon[3][4][Y1Y2TERM]   = 2*C0201 + C0210 + C2001 + (2*C2010)/3;
      Upsilon[3][4][Y1TERM]     = C0201 + C2001/3;
      Upsilon[3][4][0]          = -C0201 - C0210/2 - C2001/3 - C2010/4;


      Upsilon[4][0][0]          = C2010/4 + M11/6 + M13/3 + M33/6 + S + V1/2 + V3/2;
    
      Upsilon[4][1][Y1Y2TERM]   = -C2010/3 - M13/2 - M33/2 - V3;
      Upsilon[4][1][Y1TERM]     = -C2001/3 - M14/2 - M34/2 - V4;
      Upsilon[4][1][0]          = C2001/3 + M12/2 + M14/2 + M23/2 + M34/2 - S + V2 + V4;
    
      Upsilon[4][2][Y12Y22TERM] = M33/2;
      Upsilon[4][2][Y12Y2TERM]  = M34;
      Upsilon[4][2][Y12TERM]    = M44/2;
      Upsilon[4][2][Y1Y2TERM]   = -M23 - M34 + V3;
      Upsilon[4][2][Y1TERM]     = -M24 - M44 + V4;
      Upsilon[4][2][0]          = C0210/2 + M22/2 + M24 + M44/2 - V1/2 - V2 - V3/2 - V4;
    
      Upsilon[4][3][Y12Y22TERM] = -M33/2;
      Upsilon[4][3][Y12Y2TERM]  = -M34;
      Upsilon[4][3][Y12TERM]    = -M44/2;
      Upsilon[4][3][Y1Y2TERM]   = -C0210 + M13/2 + M23 + M33/2 + M34;
      Upsilon[4][3][Y1TERM]     = -C0201 + M14/2 + M24 + M34/2 + M44;
      Upsilon[4][3][0]          = C0201 - M11/6 - M12/2 - M13/3 - M14/2 - M22/2 - M23/2 - M24 - M33/6 - M34/2 - M44/2;
    
      Upsilon[4][4][Y1Y2TERM]   = C0210 + C2010/3;
      Upsilon[4][4][Y1TERM]     = C0201 + C2001/3;
      Upsilon[4][4][0]          = -C0201 - C0210/2 - C2001/3 - C2010/4;


      Upsilon[5][0][0]          = C2010/4 + M11/6 + M13/3 + M33/6 + S + V1/2 + V3/2;
    
      Upsilon[5][1][Y1Y2TERM]   = (-2*C2010)/3 - M11/2 - M13/2 - V1;
      Upsilon[5][1][Y1TERM]     = -M12/2 - M23/2 - V2;
      Upsilon[5][1][0]          = C2001/3 + M12/2 + M14/2 + M23/2 + M34/2 - S + V2 + V4;
    
      Upsilon[5][2][Y12Y22TERM] = C2010/2 + M11/2;
      Upsilon[5][2][Y12Y2TERM]  = M12;
      Upsilon[5][2][Y12TERM]    = C0210/2 + M22/2;
      Upsilon[5][2][Y1Y2TERM]   = -C2001 - M12 - M14 + V1;
      Upsilon[5][2][Y1TERM]     = -C0210 - M22 - M24 + V2;
      Upsilon[5][2][0]          = C0210/2 + M22/2 + M24 + M44/2 - V1/2 - V2 - V3/2 - V4;
    
      Upsilon[5][3][Y12Y22TERM] = C2001 - M11/2;
      Upsilon[5][3][Y12Y2TERM]  = -M12;
      Upsilon[5][3][Y12TERM]    = C0201 - M22/2;
      Upsilon[5][3][Y1Y2TERM]   = M11/2 + M12 + M13/2 + M14;
      Upsilon[5][3][Y1TERM]     = -2*C0201 + M12/2 + M22 + M23/2 + M24;
      Upsilon[5][3][0]          = C0201 - M11/6 - M12/2 - M13/3 - M14/2 - M22/2 - M23/2 - M24 - M33/6 - M34/2 - M44/2;
    
      Upsilon[5][4][Y12Y22TERM] = -C2001 - C2010/2;
      Upsilon[5][4][Y12TERM]    = -C0201 - C0210/2;
      Upsilon[5][4][Y1Y2TERM]   = C2001 + (2*C2010)/3;
      Upsilon[5][4][Y1TERM]     = 2*C0201 + C0210;
      Upsilon[5][4][0]          = -C0201 - C0210/2 - C2001/3 - C2010/4;

   }
  else //if (WhichCase == TD_COMMONVERTEX)
   {
     // d=1
     Upsilon[0][0][0] = S;

     Upsilon[0][1][0]           = V1; 
     Upsilon[0][1][Y1TERM]      = V2; 
     Upsilon[0][1][Y2TERM]      = V3; 
     Upsilon[0][1][Y2Y3TERM]    = V4; 

     Upsilon[0][2][0]           = M11/2.0;
     Upsilon[0][2][Y1TERM]      = M12;
     Upsilon[0][2][Y12TERM]     = M22/2.0;
     Upsilon[0][2][Y2TERM]      = M13;
     Upsilon[0][2][Y1Y2TERM]    = M23;
     Upsilon[0][2][Y22TERM]     = M33/2.0;
     Upsilon[0][2][Y2Y3TERM]    = M14;
     Upsilon[0][2][Y1Y2Y3TERM]  = M24;
     Upsilon[0][2][Y22Y3TERM]   = M34;
     Upsilon[0][2][Y22Y32TERM]  = M44/2.0;

     Upsilon[0][3][Y2TERM]      = C2010;
     Upsilon[0][3][Y12Y2TERM]   = C0210;
     Upsilon[0][3][Y2Y3TERM]    = C2001;
     Upsilon[0][3][Y12Y2Y3TERM] = C0201;

     // d=2
     Upsilon[1][0][0]           = S;

     Upsilon[1][1][0]           = V3;
     Upsilon[1][1][Y1TERM]      = V4;
     Upsilon[1][1][Y2TERM]      = V1;
     Upsilon[1][1][Y2Y3TERM]    = V2;

     Upsilon[1][2][0]           = M33/2.0;
     Upsilon[1][2][Y1TERM]      = M34;
     Upsilon[1][2][Y12TERM]     = M44/2.0;
     Upsilon[1][2][Y2TERM]      = M13;
     Upsilon[1][2][Y1Y2TERM]    = M14;
     Upsilon[1][2][Y22TERM]     = M11/2.0;
     Upsilon[1][2][Y2Y3TERM]    = M23;
     Upsilon[1][2][Y1Y2Y3TERM]  = M24;
     Upsilon[1][2][Y22Y3TERM]   = M12;
     Upsilon[1][2][Y22Y32TERM]  = M22/2.0;

     Upsilon[1][3][Y22TERM]      = C2010;
     Upsilon[1][3][Y1Y22TERM]    = C2001;
     Upsilon[1][3][Y22Y32TERM]   = C0210;
     Upsilon[1][3][Y1Y22Y32TERM] = C0201;
   };

} // CMVStoUpsilon

/***************************************************************/
/***************************************************************/
/***************************************************************/
void ComputeGeometricParameters(TaylorDuffyArgStruct *Args, 
                                TDWorkspace *TDW)
{ 
  double *V1       = Args->V1;
  double *V2       = Args->V2;
  double *V3       = Args->V3;
  double *V2P      = Args->V2P;
  double *V3P      = Args->V3P;
  double *Q        = Args->Q;
  double *QP       = Args->QP;
  double *nHat     = Args->nHat;
   
  /***************************************************************/
  /* the manual claims that certain parameters are not referenced*/
  /* in certain cases, which means the users might pass NULL for */
  /* those parameters, which could cause core dumps unless we do */
  /* the following                                               */
  /***************************************************************/
  if (Args->WhichCase==TD_COMMONTRIANGLE)
   { V2P=V2; V3P=V3; } 
  else if (Args->WhichCase==TD_COMMONEDGE)
   V2P=V2;

  double A[3], B[3], AP[3], BP[3], L[3];
  VecSub(V2,V1,A);
  VecSub(V3,V2,B);
  VecSub(V2P,V1,AP);
  VecSub(V3P,V2P,BP);

  TDW->A2    = VecDot(A,A);
  TDW->AdB   = VecDot(A,B);
  TDW->B2    = VecDot(B,B);

  if ( Args->WhichCase == TD_COMMONEDGE )
   { 
     VecSub(BP,B,L);
     TDW->BP2   = VecDot(BP,BP);
     TDW->AdBP  = VecDot(A,BP);
     TDW->L2    = VecDot(L,L);
     TDW->AdBP  = VecDot(A,BP);
     TDW->AdL   = VecDot(A,L);
     TDW->BPdL  = VecDot(BP,L);
   }
  else if ( Args->WhichCase == TD_COMMONVERTEX)
   {
     TDW->AP2   = VecDot(AP,AP);
     TDW->BP2   = VecDot(BP,BP);
     TDW->AdAP  = VecDot(A,AP);
     TDW->AdBP  = VecDot(A,BP);
     TDW->BdAP  = VecDot(B,AP);
     TDW->BdBP  = VecDot(B,BP);
     TDW->APdBP = VecDot(AP,BP);
   };

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  int NumRegions=0;
  switch(Args->WhichCase)
   { case TD_COMMONVERTEX:   NumRegions=2;  break;
     case TD_COMMONEDGE:     NumRegions=6;  break;
     case TD_COMMONTRIANGLE: NumRegions=3;  break;
   };
   
  /***************************************************************/
  /* compute the C, M, V, S parameters and convert to Upsilon    */
  /* vectors for each P polynomial we need.                      */
  /***************************************************************/
  double D[3], DP[3], DeltaD[3], AT[3], BT[3], DT[3], Scratch[3];
  if ( Q && QP) 
   { VecSub(V1,Q,D);
     VecSub(V1,QP,DP);
     VecSub(QP, Q, DeltaD);
   };
  if(nHat)
   { VecCross(nHat, A, AT);
     VecCross(nHat, B, BT);
     VecCross(nHat, D, DT);
   };

  for(int np=0; np<NUMPS; np++)
   { 
     if ( TDW->NeedP[np] == false )
      continue;
   
     double C[4], M[16], V[4], S;
     memset(C, 0,  4*sizeof(double));
     memset(M, 0, 16*sizeof(double));
     memset(V, 0,  4*sizeof(double));
     S=0.0;
     switch(np)
      { 
        case TD_UNITY: 
          S=1; 
          break;
   
        case TD_RNORMAL:
          V[0] =  VecDot(nHat, A);
          V[1] =  VecDot(nHat, B);
          V[2] = -VecDot(nHat, AP);
          V[3] = -VecDot(nHat, BP);
          break;
   
        case TD_PMCHWG1:
          M[ 0 + 2*4 ] = VecDot(A,AP); 
          M[ 0 + 3*4 ] = VecDot(A,BP);
          M[ 1 + 2*4 ] = VecDot(B,AP);
          M[ 1 + 3*4 ] = VecDot(B,BP);
          V[0] =  VecDot(A, DP);
          V[1] =  VecDot(B, DP);
          V[2] =  VecDot(D, AP);
          V[3] =  VecDot(D, BP);
          S    =  VecDot(D, DP); // + 4.0/(IK*IK);
          break;
   
        case TD_PMCHWC:
          M[ 0 + 2*4 ] = -VecDot(VecCross(A,AP,Scratch), DeltaD);
          M[ 0 + 3*4 ] = -VecDot(VecCross(A,BP,Scratch), DeltaD);
          M[ 1 + 2*4 ] = -VecDot(VecCross(B,AP,Scratch), DeltaD);
          M[ 1 + 3*4 ] = -VecDot(VecCross(B,BP,Scratch), DeltaD);
          V[0] = -VecDot(A, VecCross(DP,D,Scratch) );
          V[1] = -VecDot(B, VecCross(DP,D,Scratch) );
          V[2] =  VecDot(AP, VecCross(DP,D,Scratch) );
          V[3] =  VecDot(BP, VecCross(DP,D,Scratch) );
          break;
   
        case TD_NMULLERG1:
          M[ 0 + 2*4 ] = VecDot(AT, AP);
          M[ 0 + 3*4 ] = VecDot(AT, BP);
          M[ 1 + 2*4 ] = VecDot(BT, AP);
          M[ 1 + 3*4 ] = VecDot(BT, BP);
          V[0] =  VecDot(AT, DP);
          V[1] =  VecDot(BT, DP);
          V[2] =  VecDot(DT, AP);
          V[3] =  VecDot(DT, BP);
          S = VecDot(DT, DP);
          break;
   
        case TD_NMULLERG2:
          M[ 0 + 0*4 ] = -4.0*VecDot(AT, A);
          M[ 1 + 1*4 ] = -4.0*VecDot(BT, B);
          M[ 0 + 2*4 ] = 2.0*VecDot(AT, AP);
          M[ 0 + 3*4 ] = 2.0*VecDot(AT, BP);
          M[ 1 + 2*4 ] = 2.0*VecDot(BT, AP);
          M[ 1 + 3*4 ] = 2.0*VecDot(BT, BP);
          V[0] =  0.0;
          V[1] =  0.0;
          V[2] =  2.0*VecDot(DT,AP);
          V[3] =  2.0*VecDot(DT,BP);
          break;
   
        case TD_NMULLERC:
          C[0] = VecDot(AT, VecCross(A, AP, Scratch) );
          C[1] = VecDot(AT, VecCross(A, BP, Scratch) );
          C[2] = VecDot(BT, VecCross(B, AP, Scratch) );
          C[3] = VecDot(BT, VecCross(B, AP, Scratch) );
          M[ 0 + 0*4 ] =  VecDot(AT, VecCross(A,  DP, Scratch) );
          M[ 0 + 2*4 ] = -VecDot(AT, VecCross(AP, DP, Scratch) );
          M[ 0 + 3*4 ] = -VecDot(AT, VecCross(BP, DP, Scratch) );
          M[ 1 + 1*4 ] =  VecDot(BT, VecCross(B,  DP, Scratch) );
          M[ 1 + 2*4 ] = -VecDot(BT, VecCross(AP, DP, Scratch) );
          M[ 1 + 3*4 ] = -VecDot(BT, VecCross(BP, DP, Scratch) );
          V[0] =   VecDot(DT, VecCross(A,  DP, Scratch) );
          V[1] =   VecDot(DT, VecCross(B,  DP, Scratch) );
          V[2] =  -VecDot(DT, VecCross(AP, DP, Scratch) );
          V[3] =  -VecDot(DT, VecCross(BP, DP, Scratch) );
          break;

        case TD_EPPFT1:
          V[2] = 2.0*VecDot(nHat, AP);
          V[3] = 2.0*VecDot(nHat, BP);
          S    = 2.0*VecDot(nHat, DP);
          break;

        case TD_EPPFT2:
          V[0] =  4.0*VecDot(nHat, A );
          V[1] =  4.0*VecDot(nHat, B );
          V[2] = -4.0*VecDot(nHat, AP);
          V[3] = -4.0*VecDot(nHat, BP);
          break;

        case TD_EPPFT3:
          M[0 + 2*4] = 2.0*VecDot(nHat, VecCross(A,AP,Scratch) );
          M[0 + 3*4] = 2.0*VecDot(nHat, VecCross(A,BP,Scratch) );
          M[1 + 2*4] = 2.0*VecDot(nHat, VecCross(B,AP,Scratch) );
          M[1 + 3*4] = 2.0*VecDot(nHat, VecCross(B,BP,Scratch) );
          M[2 + 0*4] = 2.0*VecDot(nHat, VecCross(A,AP,Scratch) );
          M[2 + 1*4] = 2.0*VecDot(nHat, VecCross(B,AP,Scratch) );
          M[3 + 0*4] = 2.0*VecDot(nHat, VecCross(A,BP,Scratch) );
          M[3 + 1*4] = 2.0*VecDot(nHat, VecCross(B,BP,Scratch) );
          V[0] =  2.0*VecDot(nHat, VecCross(A, DP,Scratch));
          V[1] =  2.0*VecDot(nHat, VecCross(B, DP,Scratch));
          V[2] = -2.0*VecDot(nHat, VecCross(AP,DP,Scratch));
          V[3] = -2.0*VecDot(nHat, VecCross(BP,DP,Scratch));
          break;

        case TD_EPPFT4:
          M[0 + 2*4] = VecDot(nHat, VecCross(A,AP,Scratch) );
          M[0 + 3*4] = VecDot(nHat, VecCross(A,BP,Scratch) );
          M[1 + 2*4] = VecDot(nHat, VecCross(B,AP,Scratch) );
          M[1 + 3*4] = VecDot(nHat, VecCross(B,BP,Scratch) );
          M[2 + 0*4] = VecDot(nHat, VecCross(A,AP,Scratch) );
          M[2 + 1*4] = VecDot(nHat, VecCross(B,AP,Scratch) );
          M[3 + 0*4] = VecDot(nHat, VecCross(A,BP,Scratch) );
          M[3 + 1*4] = VecDot(nHat, VecCross(B,BP,Scratch) );
          V[0] =  VecDot(nHat, VecCross(A, DP,Scratch));
          V[1] =  VecDot(nHat, VecCross(B, DP,Scratch));
          V[2] =  VecDot(nHat, VecCross(D, AP,Scratch));
          V[3] =  VecDot(nHat, VecCross(D, BP,Scratch));
          S    =  VecDot(nHat, VecCross(D, DP,Scratch));
          break;

        case TD_EPPFT5:
          M[0 + 2*4] = -2.0*VecDot(nHat, VecCross(A,AP,Scratch) );
          M[0 + 3*4] = -2.0*VecDot(nHat, VecCross(A,BP,Scratch) );
          M[1 + 2*4] = -2.0*VecDot(nHat, VecCross(B,AP,Scratch) );
          M[1 + 3*4] = -2.0*VecDot(nHat, VecCross(B,BP,Scratch) );
          M[2 + 0*4] = -2.0*VecDot(nHat, VecCross(A,AP,Scratch) );
          M[2 + 1*4] = -2.0*VecDot(nHat, VecCross(B,AP,Scratch) );
          M[3 + 0*4] = -2.0*VecDot(nHat, VecCross(A,BP,Scratch) );
          M[3 + 1*4] = -2.0*VecDot(nHat, VecCross(B,BP,Scratch) );
          V[0] =  2.0*VecDot(nHat, VecCross(D, A,Scratch));
          V[1] =  2.0*VecDot(nHat, VecCross(D, B,Scratch));
          V[2] = -2.0*VecDot(nHat, VecCross(D, AP,Scratch));
          V[3] = -2.0*VecDot(nHat, VecCross(D, BP,Scratch));
          break;

        case TD_EPPFT6:
          break;

      };
   
     /***************************************************************/
     /* convert C, M, V, S into Upsilon coefficients, then          */
     /* determine the lower and upper limits of the sum over n,     */
     /* and the number of nonzero monomial coefficients in the      */
     /* Upsilon vector                                              */
     /* The check on the magnitude of Upsilon below should really   */
     /* involve a lengthscale these the Upsilon values are          */
     /* dimensionful quantities; but                                */
     /***************************************************************/
     CMVStoUpsilon(Args->WhichCase, C, M, V, S, TDW->Upsilon[np]);
   
     int nMin=NUMWPOWERS-1;
     int nMax=-1;
     int MaxMonomial=-1;
     for(int n=0; n<NUMWPOWERS; n++)
      for(int nr=0; nr<NumRegions; nr++)
       for(int nm=0; nm<NUMMONOMIALS; nm++)
 //       if ( TDW->Upsilon[np][nr][n][nm] != 0.0 )
        if ( fabs(TDW->Upsilon[np][nr][n][nm]) > 1.0e-8 )
         { if (n<nMin) nMin=n;
           if (n>nMax) nMax=n;
           if (nm>MaxMonomial) MaxMonomial=nm;
         };
     TDW->nMin[np]=nMin;
     TDW->nMax[np]=nMax;
     TDW->MaxMonomial[np]=MaxMonomial;

   }; // for(int np=0; np<NumPS; np++)

} // GetGeometricParameters routine

} // namespace scuff
