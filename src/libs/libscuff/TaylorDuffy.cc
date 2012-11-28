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
 *                  evaluating panel-panel integrals over pairs
 *                  of panels with common vertices
 *
 * homer reid       10/2012
 *
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

/********************************************************************/
/* nitty-gritty subroutines, implemented at the bottom of this file */
/********************************************************************/
static void GetAlphaBetaGamma2(TMWorkspace *TMW, const double *y, 
                               double *A, double *B, double *G2);

static void GetX(TMWorkspace *TMW, const double *yVector, double *X);

static void GetScriptP(TMWorkspace *TMW, int WhichP, const double *yVector, 
                       int *nMin, int *nMax, cdouble P[7][5][2]);

static void GetScriptJL(int WhichK, cdouble KParam,
                        double Alpha, double Beta, double Gamma,
                        int nMin, int nMax, cdouble JVector[7], cdouble LVector[7]);

static void GetScriptK(int WhichK, cdouble KParam, double X,
                       int nMin, int nMax, cdouble KVector[7]);

/***************************************************************/
/***************************************************************/
/***************************************************************/
void TaylorDuffySum(unsigned ndim, const double *yVector, void *parms, 
                    unsigned nfun, double *f)
{
  (void) ndim;
  (void) nfun;

  TMWorkspace *TMW = (TMWorkspace *)parms;

  int WhichCase       = TMW->WhichCase;
  int TwiceIntegrable = TMW->TwiceIntegrable;
  int NumPKs          = TMW->NumPKs;
  int *PIndex         = TMW->PIndex;
  int *KIndex         = TMW->KIndex;
  cdouble *KParam     = TMW->KParam;

  TMW->nCalls++;

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  int NumRegions, nOffset;
  double y, Jacobian;
  if (WhichCase==TM_COMMONTRIANGLE)
   { NumRegions=3;
     nOffset=1;
     y = (TwiceIntegrable ? 0.0 : yVector[0]);
     Jacobian=1.0;
   }
  else if (WhichCase==TM_COMMONEDGE)
   { NumRegions=6;
     nOffset=2;
     y = yVector[1];
     Jacobian=yVector[0];
   }
  else // (WhichCase==TM_COMMONVERTEX)
   { NumRegions=2;
     nOffset=3;
     y = yVector[2];
     Jacobian=yVector[1];
   };

  /*--------------------------------------------------------------*/
  /*- prefetch values of the X function (once-integrable case) or */
  /*- the Alpha, Beta, Gamma coefficients (twice-integrable case) */
  /*- for all subregions.                                         */
  /*--------------------------------------------------------------*/
  double X[6], A[6], B[6], G2[6];
  if (TwiceIntegrable)
   GetAlphaBetaGamma2(TMW, yVector, A, B, G2);
  else
   GetX(TMW, yVector, X);

  /*--------------------------------------------------------------*/
  /*- loop over all Ps/Ks                                        -*/
  /*--------------------------------------------------------------*/
  cdouble P[6][5][2], J[6][7], L[6][7], K[6][7];
  int nMin, nMax;
  cdouble *Sum=(cdouble *)f;
  for(int npk=0; npk<NumPKs; npk++)
   { 
     TMW->CurrentKParam=TMW->KParam[npk]; // this is poor programming style
 
     // get the \mathcal{P} coefficients for this P polynomial
     GetScriptP(TMW, PIndex[npk], yVector, &nMin, &nMax, P);
  
     // get the first or second integrals (the \mathcal{K} or
     // \mathcal{J,L} functions) for this kernel. (Note we only need 
     // to recompute if the kernel has changed since the last loop iteration).
     // 20121031 whoops! the values of nMin and nMax may have changed
     // even though the kernel itself hasn't changed.
     if ( 1 ) //npk==0 || KIndex[npk]!=KIndex[npk-1] || KParam[npk]!=KParam[npk-1] )
      { if (TwiceIntegrable)
         for(int d=0; d<NumRegions; d++)
          GetScriptJL( KIndex[npk], KParam[npk], A[d], B[d], G2[d], 
                       nMin+nOffset, nMax+nOffset, J[d], L[d]);
        else // once integrable
         for(int d=0; d<NumRegions; d++)
          GetScriptK( KIndex[npk], KParam[npk], X[d], nMin+nOffset, nMax+nOffset, K[d]);
      };

     // add contributions of all subregions and all n-values 
     Sum[npk]=0.0;
     if (TwiceIntegrable)
      for(int n=nMin; n<=nMax; n++)
       for(int d=0; d<NumRegions; d++)
        Sum[npk] += P[d][n][0]*J[d][n+nOffset] + P[d][n][1]*L[d][n+nOffset];
     else // once integrable
      for(int n=nMin; n<=nMax; n++)
       for(int d=0; d<NumRegions; d++)
        Sum[npk] += (P[d][n][0] + y*P[d][n][1]) * K[d][n+nOffset];

     Sum[npk] *= Jacobian/(4.0*M_PI);

   };

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void InitTaylorDuffyArgs(TaylorDuffyArgStruct *Args)
{
  Args->Q=0;
  Args->QP=0;  
  Args->AbsTol=0.0;     // DEFABSTOL;
  Args->RelTol=1.0e-10; // DEFRELTOL;
  Args->MaxEval=10000;  // DEFRELTOL;
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

  double *V1       = Args->V1;
  double *V2       = Args->V2;
  double *V3       = Args->V3;
  double *V2P      = Args->V2P;
  double *V3P      = Args->V3P;
  double *Q        = Args->Q;
  double *QP       = Args->QP;

  double AbsTol    = Args->AbsTol;
  double RelTol    = Args->RelTol;
  double MaxEval   = Args->MaxEval;
   
  /***************************************************************/
  /* the manual claims that certain parameters are not referenced*/
  /* in certain cases, which means the users might pass NULL for */
  /* those parameters, which could cause core dumps unless we do */
  /* the ollowing                                                */
  /***************************************************************/
  if (WhichCase==TM_COMMONTRIANGLE)
   { V2P=V2; V3P=V3; } 
  else if (WhichCase==TM_COMMONEDGE)
   V2P=V2;

  /***************************************************************/
  /* initialize TMW structure to pass data to integrand routines */
  /***************************************************************/
  TMWorkspace MyTMW, *TMW=&MyTMW;
  TMW->WhichCase = WhichCase;
  TMW->NumPKs    = NumPKs;
  TMW->PIndex    = PIndex;
  TMW->KIndex    = KIndex;
  TMW->KParam    = KParam;

  /***************************************************************/
  /* compute geometric parameters                                */
  /***************************************************************/
  double A[3], AP[3], B[3], BP[3], L[3], D[3], DP[3], QmQP[3], QxQP[3], TV[3];

  VecSub(V2,V1,A);
  VecSub(V3,V2,B);
  VecSub(V2P,V1,AP);
  VecSub(V3P,V2P,BP);
  VecSub(BP,B,L);

  TMW->A2    = VecDot(A,A);
  TMW->B2    = VecDot(B,B);
  TMW->AP2   = VecDot(AP,AP);
  TMW->BP2   = VecDot(BP,BP);
  TMW->L2    = VecDot(L,L);
  TMW->AdB   = VecDot(A,B);
  TMW->AdAP  = VecDot(A,AP);
  TMW->AdBP  = VecDot(A,BP);
  TMW->AdL   = VecDot(A,L);
  TMW->BdAP  = VecDot(B,AP);
  TMW->BdBP  = VecDot(B,BP);
  TMW->APdBP = VecDot(AP,BP);
  TMW->BPdL  = VecDot(BP,L);

  int npk;
  int NeedDot=0;
  for(NeedDot=0, npk=0; NeedDot==0 && npk<NumPKs; npk++)
   if ( PIndex[npk]==TM_DOT || PIndex[npk]==TM_DOTPLUS )
    NeedDot=1;
  if (NeedDot)
   {
     VecSub(V1,Q,D);
     VecSub(V1,QP,DP);
     TMW->AdD   = VecDot(A,D);
     TMW->AdDP  = VecDot(A,DP);
     TMW->BdDP  = VecDot(B,DP);
     TMW->APdD  = VecDot(AP,D);
     TMW->BPdD  = VecDot(BP,D);
     TMW->DdDP  = VecDot(D,DP);
   };

  int NeedCross=0;
  for(NeedCross=0, npk=0; NeedCross==0 && npk<NumPKs; npk++)
   if ( PIndex[npk]==TM_CROSS )
    NeedCross=1;
  if (NeedCross)
   { 
     VecSub(Q,QP,QmQP);
     VecCross(Q,QP,QxQP);

     TMW->AdQxQP     = VecDot(A,  QxQP );
     TMW->APdQxQP    = VecDot(AP, QxQP );
     TMW->BdQxQP     = VecDot(B,  QxQP );
     TMW->BPdQxQP    = VecDot(BP, QxQP );
     TMW->LdQxQP     = VecDot(L,  QxQP );
     TMW->V1xAdQmQP  = VecDot( VecCross(V1,A,TV),  QmQP );
     TMW->V1xAPdQmQP = VecDot( VecCross(V1,AP,TV), QmQP );
     TMW->V1xBdQmQP  = VecDot( VecCross(V1,B,TV),  QmQP );
     TMW->V1xBPdQmQP = VecDot( VecCross(V1,BP,TV), QmQP );
     TMW->AxAPdQmQP  = VecDot( VecCross(A,AP,TV),  QmQP );
     TMW->AxBdQmQP   = VecDot( VecCross(A,B,TV),   QmQP );
     TMW->AxBPdQmQP  = VecDot( VecCross(A,BP,TV),  QmQP );
     TMW->BxAPdQmQP  = VecDot( VecCross(B,AP,TV),  QmQP );
     TMW->BxBPdQmQP  = VecDot( VecCross(B,BP,TV),  QmQP );

   };

  /***************************************************************/
  /* assume we are twice integrable and check for otherwise      */
  /***************************************************************/
  int TwiceIntegrable=1;
  if (Args->ForceOnceIntegrable) TwiceIntegrable=0;
  for(npk=0; TwiceIntegrable==1 && npk<NumPKs; npk++)
   if (KIndex[npk]==TM_HELMHOLTZ || KIndex[npk]==TM_GRADHELMHOLTZ)
    TwiceIntegrable=0;
  TMW->TwiceIntegrable=TwiceIntegrable;

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
  TMW->nCalls=0;
  int IntegralDimension = 4 - WhichCase - TwiceIntegrable;

  if (IntegralDimension==0)
   TaylorDuffySum(0, 0, (void *)TMW, 0, dResult);
  else
   adapt_integrate(fDim, TaylorDuffySum, (void *)TMW, IntegralDimension, 
                   Lower, Upper, MaxEval, AbsTol, RelTol, dResult, dError);

  Args->nCalls = TMW->nCalls;

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
static void GetAlphaBetaGamma2(TMWorkspace *TMW, const double *yVector,
                               double *AVector, double *BVector, double *G2Vector)
{
  if (TMW->WhichCase==TM_COMMONTRIANGLE) 
   { 
     double A2  = TMW->A2; 
     double AdB = TMW->AdB;
     double B2  = TMW->B2; 

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
  else if (TMW->WhichCase==TM_COMMONEDGE) 
   { 
     double y1 = yVector[0], y12=y1*y1;

     double A2   = TMW->A2;
     double BP2  = TMW->BP2;
     double L2   = TMW->L2;
     double AdBP = TMW->AdBP;
     double AdL  = TMW->AdL;
     double BPdL = TMW->BPdL;

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
  else // (WhichCase==TM_COMMONVERTEX) 
   {
     double y1 = yVector[0];
     double y2 = yVector[1], y22=y2*y2;

     double A2    = TMW->A2;
     double B2    = TMW->B2;
     double AP2   = TMW->AP2;
     double BP2   = TMW->BP2;
     double AdB   = TMW->AdB;
     double AdAP  = TMW->AdAP;
     double AdBP  = TMW->AdBP;
     double BdAP  = TMW->BdAP;
     double BdBP  = TMW->BdBP;
     double APdBP = TMW->APdBP;

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
static void GetX(TMWorkspace *TMW, const double *yVector, double *X)
{
  if (TMW->WhichCase==TM_COMMONTRIANGLE)
   { double y=yVector[0], u1, u2;
     u1=1.0; u2=y;     X[0]=sqrt( TMW->A2*u1*u1 + 2.0*TMW->AdB*u1*u2 + TMW->B2*u2*u2);
     u1=y;   u2=(y-1); X[1]=sqrt( TMW->A2*u1*u1 + 2.0*TMW->AdB*u1*u2 + TMW->B2*u2*u2);
     u1=y;   u2=1.0;   X[2]=sqrt( TMW->A2*u1*u1 + 2.0*TMW->AdB*u1*u2 + TMW->B2*u2*u2);
   }
  else if (TMW->WhichCase==TM_COMMONEDGE)
   { 
     double y1=yVector[0], y2=yVector[1], u1, u2, xi2;
     u1=-y1;    u2=-y1*y2;       xi2=(1.0-y1+y1*y2);
     X[0]=sqrt(       u1*u1*TMW->A2 + u2*u2*TMW->BP2 + xi2*xi2*TMW->L2 
                + 2.0*(u1*(u2*TMW->AdBP + xi2*TMW->AdL) + u2*xi2*TMW->BPdL) );

     u1=y1;     u2=y1*y2;        xi2=(1.0-y1);
     X[1]=sqrt(       u1*u1*TMW->A2 + u2*u2*TMW->BP2 + xi2*xi2*TMW->L2 
                + 2.0*(u1*(u2*TMW->AdBP + xi2*TMW->AdL) + u2*xi2*TMW->BPdL) );

     u1=-y1*y2; u2=y1*(1.0-y2);  xi2=(1.0-y1);
     X[2]=sqrt(       u1*u1*TMW->A2 + u2*u2*TMW->BP2 + xi2*xi2*TMW->L2 
                + 2.0*(u1*(u2*TMW->AdBP + xi2*TMW->AdL) + u2*xi2*TMW->BPdL) );

     u1=y1*y2;  u2=-y1*(1.0-y2); xi2=(1.0-y1*y2);
     X[3]=sqrt(       u1*u1*TMW->A2 + u2*u2*TMW->BP2 + xi2*xi2*TMW->L2 
                + 2.0*(u1*(u2*TMW->AdBP + xi2*TMW->AdL) + u2*xi2*TMW->BPdL) );

     u1=-y1*y2; u2=-y1;          xi2=1.0;
     X[4]=sqrt(       u1*u1*TMW->A2 + u2*u2*TMW->BP2 + xi2*xi2*TMW->L2 
                + 2.0*(u1*(u2*TMW->AdBP + xi2*TMW->AdL) + u2*xi2*TMW->BPdL) );

     u1=y1*y2;  u2=y1;           xi2=(1.0-y1);
     X[5]=sqrt(       u1*u1*TMW->A2 + u2*u2*TMW->BP2 + xi2*xi2*TMW->L2 
                + 2.0*(u1*(u2*TMW->AdBP + xi2*TMW->AdL) + u2*xi2*TMW->BPdL) );

   }
  else // (TMW->WhichCase==TM_COMMONVERTEX)
   { double y1=yVector[0], y2=yVector[1], y3=yVector[2], xi1, xi2, eta1, eta2;

     xi1=1.0; xi2=y1; eta1=y2; eta2=y2*y3;
     X[0] = sqrt(     xi1*xi1*TMW->A2    +   xi2*xi2*TMW->B2 
                    + eta1*eta1*TMW->AP2   + eta2*eta2*TMW->BP2
                    + 2.0*xi1*(xi2*TMW->AdB - eta1*TMW->AdAP - eta2*TMW->AdBP)
                    - 2.0*xi2*(eta1*TMW->BdAP + eta2*TMW->BdBP) 
                    + 2.0*eta1*eta2*TMW->APdBP
                 );

     xi1=y2; xi2=y2*y3; eta1=1.0; eta2=y1;
     X[1] = sqrt(     xi1*xi1*TMW->A2    +   xi2*xi2*TMW->B2 
                    + eta1*eta1*TMW->AP2   + eta2*eta2*TMW->BP2
                    + 2.0*xi1*(xi2*TMW->AdB - eta1*TMW->AdAP - eta2*TMW->AdBP)
                    - 2.0*xi2*(eta1*TMW->BdAP + eta2*TMW->BdBP) 
                    + 2.0*eta1*eta2*TMW->APdBP
                 );

   };
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void ScriptP_One(const double *yVec, TMWorkspace *TMW,
                 int *nMin, int *nMax, cdouble P[7][5][2])
{ 
  (void) yVec; // unused

  switch(TMW->WhichCase)
   { 
     case TM_COMMONTRIANGLE: 
       *nMin=0;
       *nMax=2; 
       P[0][0][0]=P[1][0][0]=P[2][0][0]=1.0;
       P[0][0][1]=P[1][0][1]=P[3][0][1]=0.0;

       P[0][1][0]=P[1][1][0]=P[2][1][0]=-2.0;
       P[0][1][1]=P[1][1][1]=P[2][1][1]=0.0;

       P[0][2][0]=P[1][2][0]=P[2][2][0]=1.0;
       P[0][2][1]=P[1][2][1]=P[2][2][1]=0.0;
       return; 

     case TM_COMMONEDGE: 
       *nMin=0;
       *nMax=1; 
       P[0][0][0]=P[1][0][0]=P[2][0][0]=P[3][0][0]=P[4][0][0]=P[5][0][0]=1.0;
       P[0][0][1]=P[1][0][1]=P[2][0][1]=P[3][0][1]=P[4][0][1]=P[5][0][1]=0.0;

       P[0][1][0]=P[1][1][0]=P[2][1][0]=P[3][1][0]=P[4][1][0]=P[5][1][0]=-1.0;
       P[0][1][1]=P[1][1][1]=P[2][1][1]=P[3][1][1]=P[4][1][1]=P[5][1][1]=0.0;
       return;

     case TM_COMMONVERTEX: 
       *nMin=0;
       *nMax=0; 
       P[0][0][0]=P[1][0][0]=1.0;
       return;
   };

} 

/***************************************************************/
/***************************************************************/
/***************************************************************/
void ScriptP_Dot(const double *yVec, TMWorkspace *TMW,
                 int *nMin, int *nMax, cdouble P[7][5][2])
{
  double y1, y1_2, y2;
  double A2, B2, AdAP, AdB, AdBP, AdD, BdD, AdDP, APdD, BdAP, BdBP, BdDP, BPdD, DdDP;

  switch(TMW->WhichCase)
   { 
     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     case TM_COMMONTRIANGLE: 
       *nMin=0;
       *nMax=4; 

       A2=TMW->A2;
       B2=TMW->B2;
       AdB=TMW->AdB;
       AdD=TMW->APdD;
       AdDP=TMW->AdDP;
       BdD=TMW->BPdD;
       BdDP=TMW->BdDP;
       DdDP=TMW->DdDP;

       P[0][0][0] = (B2+3*A2+6*DdDP+2*BdDP+2*BdD+4*AdDP+4*AdD+3*AdB)/6.0;
       P[0][0][1] = 0.0;

       P[0][1][0] = (-4*B2-8*A2-12*DdDP-6*BdDP-6*BdD-9*AdDP -9*AdD-10*AdB)/6.0;
       P[0][1][1] = (2*B2+3*BdDP+3*BdD+4*AdB)/6.0;

       P[0][2][0] = B2+A2+DdDP+BdDP+BdD+AdDP+AdD+2*AdB;
       P[0][2][1] = (-B2-BdDP-BdD-2*AdB);

       P[0][3][0] = (-4*B2-2*BdDP-2*BdD-AdDP-AdD-6*AdB)/6.0;
       P[0][3][1] = (6*B2+3*BdDP+3*BdD+12*AdB)/6.0;

       P[0][4][0] = -(-B2+A2-AdB)/6.0;
       P[0][4][1] = -(2*B2+4*AdB)/6.0;
       
       P[1][0][0] = (B2+3*A2+6*DdDP+2*BdDP+2*BdD+4*AdDP+4*AdD+3*AdB)/6.0;
       P[1][0][1] = 0.0;
 
       P[1][1][0] = -(+2*B2+4*A2+12*DdDP+3*BdDP+3*BdD+6*AdDP+6*AdD+4*AdB)/6.0;
       P[1][1][1] = -((2*B2+4*A2+3*BdDP+3*BdD+3*AdDP+3*AdD+6*AdB))/6.0;

       P[1][2][0] = DdDP;
       P[1][2][1] = (B2+A2+BdDP+BdD+AdDP+AdD+2*AdB);

       P[1][3][0] = -(-2*B2-BdDP-BdD-2*AdDP-2*AdD)/6.0;
       P[1][3][1] = -((6*B2+3*BdDP+3*BdD+3*AdDP+3*AdD+6*AdB))/6.0;

       P[1][4][0] = (-B2+A2+AdB)/6;
       P[1][4][1] = (2*B2-2*A2)/6;
       
       P[2][0][0] = (B2+3*A2+6*DdDP+2*BdDP+2*BdD+4*AdDP+4*AdD+3*AdB)/6;
       P[2][0][1] = 0.0;

       P[2][1][0] = -( +2*B2+4*A2+12*DdDP+3*BdDP+3*BdD+6*AdDP+6*AdD+4*AdB )/6.0;
       P[2][1][1] = -(4*A2+3*AdDP+3*AdD+2*AdB)/6.0;

       P[2][2][0] = DdDP;
       P[2][2][1] = (A2+AdDP+AdD);

       P[2][3][0] = -(-2*B2-BdDP-BdD-2*AdDP-2*AdD)/6.0;
       P[2][3][1] = -((3*AdDP+3*AdD-6*AdB))/6.0;

       P[2][4][0] = -(B2-A2-AdB)/6.0;
       P[2][4][1] = -(2*A2+4*AdB)/6.0;

       return;

     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     case TM_COMMONEDGE: 
       *nMin=0;
       *nMax=3; 

       y1=yVec[0];
       y1_2=y1*y1;

       A2=TMW->A2;
       AdB=TMW->AdB;
       AdBP=TMW->AdBP;
       AdD=TMW->AdD;
       AdDP=TMW->AdDP;
       BdBP=TMW->BdBP;
       BdDP=TMW->BdDP;
       BPdD=TMW->BPdD;
       DdDP=TMW->DdDP;

       P[0][0][0] = (2*A2+6*DdDP+3*AdDP+3*AdD)/6;
       P[0][0][1] = 0.0;

       P[0][1][0] = ((-A2-2*BdDP-2*BPdD-2*AdD-AdBP-AdB)*y1
                      -2*DdDP+2*BdDP+2*BPdD+AdBP+AdB)/2;
       P[0][1][1] = (2*BdDP+AdB)*y1/2.0;

       P[0][2][0] = -((-2*BdBP-2*AdB)*y1_2+(-2*BdDP+4*BdBP-2*BPdD-2*AdD+2*AdB)*y1
                      +2*BdDP -2*BdBP+2*BPdD+AdDP+AdD)/2;
       P[0][2][1] = -((2*BdBP+2*AdB)*y1_2+(2*BdDP-2*BdBP)*y1)/2;

       P[0][3][0] = ((-6*BdBP-6*AdB)*y1_2+(3*A2+12*BdBP+3*AdBP+9*AdB)*y1
                     -2*A2-6*BdBP -3*AdBP-3*AdB)/6;
       P[0][3][1] = (((6*BdBP+6*AdB)*y1_2+(-6*BdBP-3*AdB)*y1))/6.0;
        
       P[1][0][0] = (2*A2+6*DdDP+3*AdDP+3*AdD)/6;
       P[1][0][1] = 0.0;

       P[1][1][0] = ((-A2-2*BdDP-2*BPdD-2*AdDP-AdBP-AdB)*y1
                      -2*DdDP+2*BdDP+2*BPdD+AdBP+AdB)/2;
       P[1][1][1] = (2.0*BPdD+AdBP)*y1/2.0;

       P[1][2][0] = -((-2*BdBP-2*AdBP)*y1_2+(-2*BdDP+4*BdBP-2*BPdD-2*AdDP
                       +2*AdBP)*y1+2*BdDP-2*BdBP+2*BPdD+AdDP+AdD)/2;
       P[1][2][1] = -((2*BdBP+2*AdBP)*y1_2+(-2*BdBP+2*BPdD)*y1)/2;

       P[1][3][0] = ((-6*BdBP-6*AdBP)*y1_2+(3*A2+12*BdBP+9*AdBP+3*AdB)*y1
                     -2*A2-6*BdBP-3*AdBP-3*AdB)/6;
       P[1][3][1] = (((6*BdBP+6*AdBP)*y1_2+(-6*BdBP-3*AdBP)*y1))/6.0;
        
        
       P[2][0][0] = (2*A2+6*DdDP+3*AdDP+3*AdD)/6;
       P[2][0][1] = 0.0;

       P[2][1][0] = -((2*BdDP+AdB)*y1+2*DdDP -2*BdDP-2*BPdD -AdBP-AdB)/2;
       P[2][1][1] = -((A2+2*BPdD+2*AdD+AdBP)*y1)/2;

       P[2][2][0] = ((2*BdDP-2*BdBP)*y1-2*BdDP+2*BdBP-2*BPdD-AdDP-AdD)/2;
       P[2][2][1] = ((2*BdBP+2*AdB)*y1_2+(-2*BdBP+2*BPdD+2*AdD-2*AdB)*y1)/2.0;

       P[2][3][0] = -((-6*BdBP-3*AdB)*y1+2*A2+6*BdBP+3*AdBP+3*AdB)/6;
       P[2][3][1] = -((6*BdBP+6*AdB)*y1_2+(-3*A2-6*BdBP-3*AdBP-6*AdB)*y1)/6.0;
        
       P[3][0][0] = (2*A2+6*DdDP+3*AdDP+3*AdD)/6;
       P[3][0][1] = 0.0;

       P[3][1][0] = -((2*BPdD+AdBP)*y1+2*DdDP-2*BdDP-2*BPdD-AdBP-AdB)/2;
       P[3][1][1] = -((A2+2*BdDP+2*AdDP+AdB)*y1)/2.0;

       P[3][2][0] = (+(-2*BdBP+2*BPdD)*y1-2*BdDP+2*BdBP-2*BPdD-AdDP-AdD)/2;
       P[3][2][1] = ((2*BdBP+2*AdBP)*y1_2+(2*BdDP-2*BdBP+2*AdDP-2*AdBP)*y1)/2.0;

       P[3][3][0] = -((-6*BdBP-3*AdBP)*y1+2*A2+6*BdBP+3*AdBP+3*AdB)/6;
       P[3][3][1] = -((6*BdBP+6*AdBP)*y1_2+(-3*A2-6*BdBP-6*AdBP-3*AdB)*y1)/6.0;
        
        
       P[4][0][0] = (2*A2+6*DdDP+3*AdDP+3*AdD)/6;
       P[4][0][1] = 0.0;

       P[4][1][0] = -((2*BPdD+AdBP)*y1+2*DdDP-2*BdDP-2*BPdD -AdBP-AdB)/2;
       P[4][1][1] = -(A2+2*AdD)*y1/2;

       P[4][2][0] = ((-2*BdBP+2*BPdD)*y1-2*BdDP+2*BdBP -2*BPdD-AdDP-AdD)/2; 
       P[4][2][1] = (2*AdD-2*AdB)*y1/2;

       P[4][3][0] = ((6*BdBP+3*AdBP)*y1-2*A2-6*BdBP-3*AdBP -3*AdB)/6;
       P[4][3][1] = (3*A2+6*AdB)*y1/6.0;
       
       
       P[5][0][0] = (2*A2+6*DdDP+3*AdDP+3*AdD)/6;
       P[5][0][1] = 0.0;

       P[5][1][0] = -((2*BdDP+AdB)*y1+2*DdDP-2*BdDP-2*BPdD-AdBP-AdB)/2;
       P[5][1][1] = -(A2+2*AdDP)*y1/2;

       P[5][2][0] = ((2*BdDP-2*BdBP)*y1-2*BdDP+2*BdBP -2*BPdD -AdDP-AdD)/2;
       P[5][2][1] = (2*AdDP-2*AdBP)*y1/2.0;

       P[5][3][0] = ((6*BdBP+3*AdB)*y1-2*A2-6*BdBP-3*AdBP -3*AdB)/6;
       P[5][3][1] = (3*A2+6*AdBP)*y1/6.0;
       
       return;
        
     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     case TM_COMMONVERTEX: 
       *nMin=0;
       *nMax=2;

       y1=yVec[0];
       y2=yVec[1];

       A2=TMW->A2;
       AdAP=TMW->AdAP;
       AdB=TMW->AdB;
       AdBP=TMW->AdBP;
       AdD=TMW->AdD;
       AdDP=TMW->AdDP;
       BdAP=TMW->BdAP;
       BdBP=TMW->BdBP;
       BdDP=TMW->BdDP;
       BPdD=TMW->BPdD;
       APdD=TMW->APdD;
       BPdD=TMW->BPdD;
       DdDP=TMW->DdDP;

       P[0][0][0] = DdDP;
       P[0][0][1] = 0.0;

       P[0][1][0] = APdD*y2+BdDP*y1+AdDP;
       P[0][1][1] = BPdD*y2;

       P[0][2][0] = (BdAP*y1+AdAP)*y2;
       P[0][2][1] = (BdBP*y1+AdBP)*y2;

       P[1][0][0] = DdDP;
       P[1][0][1] = 0.0;

       P[1][1][0] = AdDP*y2+BPdD*y1+APdD;
       P[1][1][1] = BdDP*y2;

       P[1][2][0] = (AdBP*y1+AdAP)*y2;
       P[1][2][1] = (BdBP*y1+BdAP)*y2;
              
       return;

   }; // switch(WhichCase)
}

/***************************************************************/
/* this is equivalent to 'Dot' - 4/(K*K) * 'One'               */
/* where K is the parameter that enters into the kernela       */
/* e^{i*K*r}/(4*pi*r)                                          */
/***************************************************************/
void ScriptP_DotPlus(const double *yVec, TMWorkspace *TMW,
                     int *nMin, int *nMax, cdouble P[7][5][2])
{

   /***************************************************************/
   /* get the 'Dot' contributions *********************************/
   /***************************************************************/
   ScriptP_Dot(yVec, TMW, nMin, nMax, P);

   /***************************************************************/
   /* add in the 'One' contributions ******************************/
   /***************************************************************/
   cdouble Factor = -4.0/(TMW->CurrentKParam*TMW->CurrentKParam);
   switch(TMW->WhichCase)
     { 
       case TM_COMMONTRIANGLE: 
         P[0][0][0] += Factor;
         P[1][0][0] += Factor;
         P[2][0][0] += Factor;

         P[0][1][0] += -2.0*Factor;
         P[1][1][0] += -2.0*Factor;
         P[2][1][0] += -2.0*Factor;

         P[0][2][0] += Factor;
         P[1][2][0] += Factor;
         P[2][2][0] += Factor;
         break;

       case TM_COMMONEDGE: 
         P[0][0][0] += Factor;
         P[1][0][0] += Factor;
         P[2][0][0] += Factor;
         P[3][0][0] += Factor;
         P[4][0][0] += Factor;
         P[5][0][0] += Factor;

         P[0][1][0] -= Factor;
         P[1][1][0] -= Factor;
         P[2][1][0] -= Factor;
         P[3][1][0] -= Factor;
         P[4][1][0] -= Factor;
         P[5][1][0] -= Factor;
         break;

       case TM_COMMONVERTEX: 
         P[0][0][0] += Factor;
         P[1][0][0] += Factor;
         break;

     }; //switch(WhichCase)

}
               
/***************************************************************/
/***************************************************************/
/***************************************************************/
void ScriptP_Cross(const double *yVec, TMWorkspace *TMW,
                   int *nMin, int *nMax, cdouble P[7][5][2])
{
  double y1, y1_2, y2;

  double AdQxQP, APdQxQP, BdQxQP, BPdQxQP, LdQxQP; 
  double V1xAdQmQP, V1xBdQmQP, V1xAPdQmQP, V1xBPdQmQP; 
  double AxAPdQmQP, BxAPdQmQP, AxBdQmQP, AxBPdQmQP, BxBPdQmQP;

  switch(TMW->WhichCase)
   { 
     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     case TM_COMMONTRIANGLE: 
       *nMin=0; 
       *nMax=-1;
       return;

     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     case TM_COMMONEDGE: 
       *nMin=1; 
       *nMax=3;

       y1=yVec[0];
       y1_2=y1*y1;

       AdQxQP     = TMW->AdQxQP;
       BPdQxQP    = TMW->BPdQxQP;
       LdQxQP     = TMW->LdQxQP;
       V1xAdQmQP  = TMW->V1xAdQmQP;
       V1xBdQmQP  = TMW->V1xBdQmQP;
       V1xBPdQmQP = TMW->V1xBPdQmQP;
       AxBdQmQP   = TMW->AxBdQmQP;
       AxBPdQmQP  = TMW->AxBPdQmQP;
       BxBPdQmQP  = TMW->BxBPdQmQP;

       P[0][1][0] = -((-2*V1xBdQmQP+2*V1xBPdQmQP+2*V1xAdQmQP-2*LdQxQP-AxBdQmQP
                      +AxBPdQmQP-2*AdQxQP)*y1+2*V1xBdQmQP-2*V1xBPdQmQP
                      +2*LdQxQP+AxBdQmQP-AxBPdQmQP)/2;

       P[0][1][1] = -((2*V1xBdQmQP+2*LdQxQP-2*BPdQxQP+AxBdQmQP)*y1)/2.0;

       P[0][2][0] = +(BxBPdQmQP-AxBdQmQP)*y1_2
                    +(-V1xBdQmQP+V1xBPdQmQP+V1xAdQmQP-LdQxQP-2*BxBPdQmQP+AxBdQmQP-AdQxQP)*y1
                    +V1xBdQmQP-V1xBPdQmQP+LdQxQP+BxBPdQmQP;

       P[0][2][1] = (-BxBPdQmQP+AxBdQmQP)*y1_2+(V1xBdQmQP+LdQxQP+BxBPdQmQP-BPdQxQP)*y1;

       P[0][3][0] = ((-2*BxBPdQmQP+2*AxBdQmQP)*y1_2+
                     (4*BxBPdQmQP-3*AxBdQmQP+AxBPdQmQP)*y1
                     -2*BxBPdQmQP+AxBdQmQP-AxBPdQmQP)/2;

       P[0][3][1] = ((2*BxBPdQmQP-2*AxBdQmQP)*y1_2+(-2*BxBPdQmQP+AxBdQmQP)*y1)/2.0;
       
       P[1][1][0] = (+(2*V1xBdQmQP-2*V1xBPdQmQP+2*V1xAdQmQP+2*LdQxQP+AxBdQmQP-AxBPdQmQP
                      -2*AdQxQP)*y1
                     -2*V1xBdQmQP+2*V1xBPdQmQP-2*LdQxQP-AxBdQmQP+AxBPdQmQP)/2;

       P[1][1][1] = (2*V1xBPdQmQP-2*BPdQxQP+AxBPdQmQP)*y1/2.0;

       P[1][2][0] = +(BxBPdQmQP+AxBPdQmQP)*y1_2
                    +(-V1xBdQmQP+V1xBPdQmQP-V1xAdQmQP-LdQxQP-2*BxBPdQmQP-AxBPdQmQP+AdQxQP)*y1
                    +V1xBdQmQP-V1xBPdQmQP+LdQxQP+BxBPdQmQP;

       P[1][2][1] = ((-BxBPdQmQP-AxBPdQmQP)*y1_2+(-V1xBPdQmQP+BxBPdQmQP+BPdQxQP)*y1);


       P[1][3][0] = ((-2*BxBPdQmQP-2*AxBPdQmQP)*y1_2+(4*BxBPdQmQP-AxBdQmQP+3*AxBPdQmQP)*y1
                      -2*BxBPdQmQP+AxBdQmQP-AxBPdQmQP)/2;
       P[1][3][1] = ((2*BxBPdQmQP+2*AxBPdQmQP)*y1_2+(-2*BxBPdQmQP-AxBPdQmQP)*y1)/2.0;
       
       P[2][1][0] = -((-2*V1xBdQmQP-2*LdQxQP+2*BPdQxQP-AxBdQmQP)*y1+2*V1xBdQmQP
                       -2*V1xBPdQmQP+2*LdQxQP+AxBdQmQP-AxBPdQmQP)/2;

       P[2][1][1] = -(2*V1xBPdQmQP+2*V1xAdQmQP-2*BPdQxQP+AxBPdQmQP-2*AdQxQP)*y1/2.0;


       P[2][2][0] = +(-V1xBdQmQP-LdQxQP-BxBPdQmQP+BPdQxQP)*y1
                    +V1xBdQmQP-V1xBPdQmQP+LdQxQP +BxBPdQmQP;

       P[2][2][1] = (BxBPdQmQP-AxBdQmQP)*y1_2
                     +(V1xBPdQmQP+V1xAdQmQP-BxBPdQmQP-BPdQxQP +AxBdQmQP-AdQxQP)*y1;

       P[2][3][0] = -((-2*BxBPdQmQP+AxBdQmQP)*y1+2*BxBPdQmQP-AxBdQmQP+AxBPdQmQP)/2;

       P[2][3][1] = -(((2*BxBPdQmQP-2*AxBdQmQP)*y1_2
                    +(-2*BxBPdQmQP+2*AxBdQmQP-AxBPdQmQP)*y1))/2;

       
       P[3][1][0] = (+(-2*V1xBPdQmQP+2*BPdQxQP-AxBPdQmQP)*y1-2*V1xBdQmQP+2*V1xBPdQmQP
                   -2*LdQxQP-AxBdQmQP+AxBPdQmQP)/2;

       P[3][1][1] = ((2*V1xBdQmQP+2*V1xAdQmQP+2*LdQxQP-2*BPdQxQP+AxBdQmQP-2*AdQxQP)*y1)/2.0;


       P[3][2][0] = +(V1xBPdQmQP-BxBPdQmQP-BPdQxQP)*y1
                    +V1xBdQmQP-V1xBPdQmQP+LdQxQP+BxBPdQmQP;

       P[3][2][1] = (BxBPdQmQP+AxBPdQmQP)*y1_2
                   +(-V1xBdQmQP-V1xAdQmQP-LdQxQP-BxBPdQmQP+BPdQxQP-AxBPdQmQP+AdQxQP)*y1;

       P[3][3][0] = -((-2*BxBPdQmQP-AxBPdQmQP)*y1+2*BxBPdQmQP-AxBdQmQP+AxBPdQmQP)/2;

       P[3][3][1] = -( (  (2*BxBPdQmQP+2*AxBPdQmQP)*y1_2
                         +(-2*BxBPdQmQP+AxBdQmQP-2*AxBPdQmQP)*y1))/2.0;

       P[4][1][0] = -((2*V1xBPdQmQP-2*BPdQxQP+AxBPdQmQP)*y1
                      +2*V1xBdQmQP-2*V1xBPdQmQP+2*LdQxQP+AxBdQmQP-AxBPdQmQP)/2;
       P[4][1][1] = -(2*V1xAdQmQP-2*AdQxQP)*y1/2.0;

       P[4][2][0] = +(V1xBPdQmQP-BxBPdQmQP-BPdQxQP)*y1
                     +V1xBdQmQP-V1xBPdQmQP+LdQxQP+BxBPdQmQP;
       P[4][2][1] = (V1xAdQmQP+AxBdQmQP-AdQxQP)*y1;

       P[4][3][0] = -((-2*BxBPdQmQP-AxBPdQmQP)*y1+2*BxBPdQmQP-AxBdQmQP+AxBPdQmQP)/2;
       P[4][3][1] = -AxBdQmQP*y1;
       
       P[5][1][0] = ( (2*V1xBdQmQP+2*LdQxQP-2*BPdQxQP +AxBdQmQP)*y1
                      -2*V1xBdQmQP+2*V1xBPdQmQP-2*LdQxQP -AxBdQmQP+AxBPdQmQP)/2;
       P[5][1][1] = (2*V1xAdQmQP-2*AdQxQP)*y1/2.0;

       P[5][2][0] = +(-V1xBdQmQP-LdQxQP-BxBPdQmQP+BPdQxQP)*y1
                    +V1xBdQmQP-V1xBPdQmQP+LdQxQP+BxBPdQmQP;

       P[5][2][1] = (-V1xAdQmQP-AxBPdQmQP+AdQxQP)*y1;

       P[5][3][0] = ((2*BxBPdQmQP-AxBdQmQP)*y1-2*BxBPdQmQP+AxBdQmQP -AxBPdQmQP)/2;

       P[5][3][1] = AxBPdQmQP*y1;

       return;

     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     case TM_COMMONVERTEX: 

       *nMin=1; 
       *nMax=2;

       y1=yVec[0];
       y2=yVec[1];

       AdQxQP     = TMW->AdQxQP;
       APdQxQP    = TMW->APdQxQP;
       BdQxQP     = TMW->BdQxQP;
       BPdQxQP    = TMW->BPdQxQP;
       V1xAdQmQP  = TMW->V1xAdQmQP;
       V1xAPdQmQP = TMW->V1xAPdQmQP;
       V1xBdQmQP  = TMW->V1xBdQmQP;
       V1xBPdQmQP = TMW->V1xBPdQmQP;
       AxAPdQmQP  = TMW->AxAPdQmQP;
       AxBPdQmQP  = TMW->AxBPdQmQP;
       BxAPdQmQP  = TMW->BxAPdQmQP;
       BxBPdQmQP  = TMW->BxBPdQmQP;

       P[0][1][0] = +(V1xAPdQmQP-APdQxQP)*y2
                    +(-V1xBdQmQP+BdQxQP)*y1-V1xAdQmQP+AdQxQP;

       P[0][1][1] = (V1xBPdQmQP-BPdQxQP)*y2;

       P[0][2][0] = (BxAPdQmQP*y1+AxAPdQmQP)*y2;
       P[0][2][1] = (BxBPdQmQP*y1+AxBPdQmQP)*y2;

       P[1][1][0] = +(-V1xAdQmQP+AdQxQP)*y2 +(V1xBPdQmQP-BPdQxQP)*y1+V1xAPdQmQP-APdQxQP;
       P[1][1][1] = (-V1xBdQmQP+BdQxQP)*y2;

       P[1][2][0] = (AxBPdQmQP*y1+AxAPdQmQP)*y2;
       P[1][2][1] = (BxBPdQmQP*y1+BxAPdQmQP)*y2;

       return;

   };
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GetScriptP(TMWorkspace *TMW, int WhichP, const double *yVector, 
                int *nMin, int *nMax, cdouble P[7][5][2])
{
   switch(WhichP)
    { case TM_ONE:     ScriptP_One(yVector, TMW, nMin, nMax, P); break;
      case TM_DOT:     ScriptP_Dot(yVector, TMW, nMin, nMax, P); break;
      case TM_DOTPLUS: ScriptP_DotPlus(yVector, TMW, nMin, nMax, P); break;
      case TM_CROSS:   ScriptP_Cross(yVector, TMW, nMin, nMax, P); break;
    };
  
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
static void GetQFPIntegral(double P, double Q2, int p, double *IntQFP, double *IntyQFP)
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

  if (WhichK==TM_RP)
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
  else if (WhichK==TM_HIGHK_HELMHOLTZ)
   { 
     for(int n=nMin; n<=nMax; n++)
      { cdouble ikAlphaMN=pow(-II*KParam*Alpha, -(double)n);
        GetQFPIntegral(Beta, Gamma2, -(n+1), &IntQFP, &IntyQFP);
        JVector[n] = FactorialTable[n-1]*ikAlphaMN*IntQFP/Alpha;
        LVector[n] = FactorialTable[n-1]*ikAlphaMN*IntyQFP/Alpha;
      };
   }
  else if (WhichK==TM_HIGHK_GRADHELMHOLTZ)
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
cdouble ExpRelV2P0(int n, cdouble Z)
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
static void GetScriptK(int WhichK, cdouble KParam, double X,
                       int nMin, int nMax, cdouble KVector[7])
{
  if (WhichK==TM_RP)
   { 
     double P=real(KParam);
     double XP=pow(X,P);
     for(int n=nMin; n<=nMax; n++)
      KVector[n] = XP / ( 1.0 + P + (double)n );
   }
  else if (WhichK==TM_HELMHOLTZ)
   { cdouble IK = II*KParam, IKX = IK*X, eIKX = exp(IKX);
     for(int n=nMin; n<=nMax; n++)
      KVector[n] = eIKX * ExpRelV2P0(n,-IKX) / (n*X);
   }
  else if (WhichK==TM_GRADHELMHOLTZ)
   { cdouble IK = II*KParam, IKX = IK*X, eIKX = exp(IKX);
     cdouble ExpRelTable[10]; 
     if ( (nMin-2) < 0 || (nMax-1) >=10) ErrExit("%s:%i: internal error",__FILE__,__LINE__);
     // precompute ExpRel for all values we need; for better efficiency,
     // maybe implement a routine that computes these values all at once?

     for(int n=nMin-2; n<=nMax-1; n++)
      ExpRelTable[n] = ExpRelV2P0(n,-IKX);

     for(int n=nMin; n<=nMax; n++)
      KVector[n] = eIKX * ( IK*ExpRelTable[n-1] / ((double)n-1.0)
                              -ExpRelTable[n-2] / (((double)n-2.0)*X)
                          ) / (X*X);
   };
  
}

} // namespace scuff
