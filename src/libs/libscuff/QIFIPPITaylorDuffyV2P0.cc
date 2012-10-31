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
 * QIFIPPITaylorDuffy.cc:  a customized version of my master TaylorDuffy
 *                         routine for the special case in which we 
 *                         want to compute all 33 of the Q-independent
 *                         FIPPIs for a given pair of panels. The  
 *                         computation could be done using the existing
 *                         TaylorDuffy routine, but this version is
 *                         specially tweaked for speed.
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

namespace scuff {

#define ABSTOL 0.0
#define RELTOL 1.0e-5
#define MAXFEVALS 10000

#define NUMFUNCS 33

int QIFIPPITaylorDuffyV2P0Calls;

/***************************************************************/
/***************************************************************/
/***************************************************************/
typedef struct TDWorkspaceFIPPI
 { 
   /* geometric data on triangles */
   double A[3], B[3], AP[3], BP[3]; 
   double V0xAP[3], V0xBP[3], AxV0[3], BxV0[3];
   double AxAP[3], AxBP[3], BxAP[3], BxBP[3];

   double A2, B2, AP2, BP2, L2;
   double AdB, AdAP, AdBP, AdL;
   double BdAP, BdBP, APdBP, BPdL;

   double PCross[6][6][4][2], PXXEE[9][6][4][2];

   int WhichCase;
   int nCalls;

 } TDWorkspaceFIPPI;


/********************************************************************/
/* nitty-gritty subroutines, implemented at the bottom of this file */
/********************************************************************/
void GetAlphaBetaGamma_FIPPI(TDWorkspaceFIPPI *W, const double *y, 
                             double *A, double *B, double *G);

void GetScriptP_CommonEdge(TDWorkspaceFIPPI *W, double y1,
                           double PCross[6][6][4][2], double PXXEE[9][6][4][2]);

void GetScriptP_CommonVertex(TDWorkspaceFIPPI *W, double y1, double y2, 
                             double PCross[6][6][4][2], double PXXEE[9][6][4][2]);

void GetQFPIntegrals_FIPPI(double P, double Q,
                           double IntQFP[4], double IntyQFP[4]);

/***************************************************************/
/***************************************************************/
/***************************************************************/
void TaylorDuffySum_FIPPI(unsigned ndim, const double *yVector, void *parms, 
                          unsigned nfun, double *f)
{
  (void) ndim; // unused 
  (void) nfun; // unused 

  TDWorkspaceFIPPI *W = (TDWorkspaceFIPPI *)parms;
  int WhichCase = W->WhichCase;
  W->nCalls++;

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  int NumRegions, nOffset;
  int nMinCross, nMaxCross, nMinXXEE, nMaxXXEE;
  double Jacobian;
  if (WhichCase==TM_COMMONEDGE)
   { NumRegions=6;
     nOffset=2;
     Jacobian=yVector[0];
     nMinCross=1;
     nMaxCross=3;
     nMinXXEE=0;
     nMaxXXEE=3;
   }
  else // (WhichCase==TM_COMMONVERTEX)
   { NumRegions=2;
     nOffset=3;
     Jacobian=yVector[1];
     nMinCross=1;
     nMaxCross=2;
     nMinXXEE=0;
     nMaxXXEE=2;
   };

  /*--------------------------------------------------------------*/
  /*- prefetch values of the the Alpha, Beta, Gamma coefficients  */
  /*- for all subregions.                                         */
  /*--------------------------------------------------------------*/
  double A[6], B[6], G2[6];
  GetAlphaBetaGamma_FIPPI(W, yVector, A, B, G2);

  /*--------------------------------------------------------------*/
  /*- precompute the \mathcal{J} and \mathcal{L} functions for    */
  /*- all regions and all kernels. there are four kernels:        */
  /*- r^{-3}, r^{-1}, r^{1}, r^{2}.                               */
  /*-                                                             */
  /*- How it works:                                               */
  /*-                                                             */
  /*-  For p=-3 (kernel r^{-3})                                   */
  /*-   JrM3[d][n] = \mathcal{J}_n(A_d, B_d, G_d)                 */
  /*-   LrM3[d][n] = \mathcal{L}_n(A_d, B_d, G_d)                 */
  /*-                                                             */
  /*-  For the kernel r^{-1}:                                     */
  /*-   Jrp[0][d][n] = \mathcal{J}_n(A_d, B_d, G_d)               */
  /*-   Lrp[0][d][n] = \mathcal{L}_n(A_d, B_d, G_d)               */
  /*-                                                             */
  /*-  For the kernel r^{1}:                                      */
  /*-   Jrp[1][d][n] = \mathcal{J}_n(A_d, B_d, G_d)               */
  /*-   Lrp[1][d][n] = \mathcal{L}_n(A_d, B_d, G_d)               */
  /*-                                                             */
  /*-  For the kernel r^{2}:                                      */
  /*-   Jrp[2][d][n] = \mathcal{J}_n(A_d, B_d, G_d)               */
  /*-   Lrp[2][d][n] = \mathcal{L}_n(A_d, B_d, G_d)               */
  /*--------------------------------------------------------------*/
  double JrM3[6][6], LrM3[6][6], Jrp[3][6][6], Lrp[3][6][6];
  double A2, A3;
  double IntQFP[4], IntyQFP[4];
  for(int d=0; d<NumRegions; d++)
   { 
     GetQFPIntegrals_FIPPI(B[d], G2[d], IntQFP, IntyQFP);
     A2=A[d]*A[d];
     A3=A2*A[d];
   
     // p=-3
     for(int n=nMinCross+nOffset; n<=nMaxCross+nOffset; n++)
      { JrM3[d][n] = IntQFP[0] / (A3*(1.0 + n - 3.0));
        LrM3[d][n] = IntyQFP[0] / (A3*(1.0 + n - 3.0));
      };

     for(int n=nMinXXEE+nOffset; n<=nMaxXXEE+nOffset; n++)
      { 
        // p=-1
        Jrp[0][d][n] = IntQFP[1] / (A[d]*(1.0 + n - 1.0));
        Lrp[0][d][n] = IntyQFP[1] / (A[d]*(1.0 + n - 1.0));

        // p=1
        Jrp[1][d][n] = A[d]*IntQFP[2] / (1.0 + n + 1.0);
        Lrp[1][d][n] = A[d]*IntyQFP[2] / (1.0 + n + 1.0);

        // p=2
        Jrp[2][d][n] = A2*IntQFP[3] / (1.0 + n + 2.0);
        Lrp[2][d][n] = A2*IntyQFP[3] / (1.0 + n + 2.0);
      };

   };

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  if (WhichCase==TM_COMMONEDGE)
   GetScriptP_CommonEdge(W, yVector[0], W->PCross, W->PXXEE);
  else
   GetScriptP_CommonVertex(W, yVector[0], yVector[1], W->PCross, W->PXXEE);

  /*--------------------------------------------------------------*/
  /*- first 6 slots in the output buffer -------------------------*/
  /*--------------------------------------------------------------*/
  memset(f,0,NUMFUNCS*sizeof(double));
  int nSum=0;
  for(int np=0; np<6; np++, nSum++)
   for(int d=0; d<NumRegions; d++)
    for(int n=nMinCross; n<=nMaxCross; n++)
     f[nSum] += Jacobian*(  W->PCross[np][d][n][0]*JrM3[d][n+nOffset] 
                          + W->PCross[np][d][n][1]*LrM3[d][n+nOffset]    );

  /*--------------------------------------------------------------*/
  /*- next 27 slots ----------------------------------------------*/
  /*--------------------------------------------------------------*/
  for(int nk=0; nk<3; nk++)
   for(int np=0; np<9; np++, nSum++)
    for(int d=0; d<NumRegions; d++)
     for(int n=nMinXXEE; n<=nMaxXXEE; n++)
      f[nSum] += Jacobian*(  W->PXXEE[np][d][n][0]*Jrp[nk][d][n+nOffset] 
                           + W->PXXEE[np][d][n][1]*Lrp[nk][d][n+nOffset] );
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void ComputeQIFIPPIData_TaylorDuffyV2P0(double *V1, double *V2, double *V3,
                                        double *V2P, double *V3P,
                                        QIFIPPIData *QIFD)
{
  TDWorkspaceFIPPI MyTDWorkspaceFIPPI, *W=&MyTDWorkspaceFIPPI;
  double Result[NUMFUNCS], Error[NUMFUNCS];
  double A[3], AP[3], B[3], BP[3], L[3];
  static double Lower[2]={0.0, 0.0};
  static double Upper[2]={1.0, 1.0};

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  for(int a=0; a<9; a++)
   for(int b=0; b<6; b++)
    for(int c=0; c<4; c++)
     for(int d=0; d<2; d++)
      W->PXXEE[a][b][c][d]=0.0;

  /***************************************************************/
  /* 1. compute geometric parameters.                            */
  /***************************************************************/
  VecSub(V2,V1,A);
  VecSub(V3,V2,B);
  VecSub(V2P,V1,AP);
  VecSub(V3P,V2P,BP);
  VecSub(BP,B,L);

  W->A2    = VecDot(A,A);
  W->B2    = VecDot(B,B);
  W->AP2   = VecDot(AP,AP);
  W->BP2   = VecDot(BP,BP);
  W->L2    = VecDot(L,L);
  W->AdB   = VecDot(A,B);
  W->AdAP  = VecDot(A,AP);
  W->AdBP  = VecDot(A,BP);
  W->AdL   = VecDot(A,L);
  W->BdAP  = VecDot(B,AP);
  W->BdBP  = VecDot(B,BP);
  W->APdBP = VecDot(AP,BP);
  W->BPdL  = VecDot(BP,L);

  memcpy(W->A,  A, 3*sizeof(double));
  memcpy(W->B,  B, 3*sizeof(double));
  memcpy(W->AP, AP, 3*sizeof(double));
  memcpy(W->BP, BP, 3*sizeof(double));
  VecCross(V1,AP,W->V0xAP);
  VecCross(V1,BP,W->V0xBP);
  VecCross(A,V1,W->AxV0);
  VecCross(B,V1,W->BxV0);
  VecCross(A,AP,W->AxAP);
  VecCross(B,AP,W->BxAP);
  VecCross(A,BP,W->AxBP);
  VecCross(B,BP,W->BxBP);

  /***************************************************************/
  /* 2. evaluate the taylor-duffy integral                       */
  /***************************************************************/
  int IntegralDimension;
  if( VecEqualFloat(V2, V2P) ) // common-edge case 
   { W->WhichCase=TM_COMMONEDGE;
     IntegralDimension=1;
   }
  else
   { W->WhichCase=TM_COMMONVERTEX;
     IntegralDimension=2;
   };
 
  W->nCalls=0;

  adapt_integrate(NUMFUNCS, TaylorDuffySum_FIPPI, (void *)W,
                  IntegralDimension, Lower, Upper,
                  MAXFEVALS, ABSTOL, RELTOL, Result, Error);

  QIFIPPITaylorDuffyV2P0Calls = W->nCalls;

  /***************************************************************/
  /* 3. pack the integrals into the appropriate slots in the     */
  /*    QIFIPPIData structure.                                   */
  /***************************************************************/
  memcpy(QIFD->xMxpRM3,   Result+0, 3*sizeof(double));
  memcpy(QIFD->xXxpRM3,   Result+3, 3*sizeof(double));
  memcpy(QIFD->uvupvpRM1, Result+6, 9*sizeof(double));
  memcpy(QIFD->uvupvpR1,  Result+15, 9*sizeof(double));
  memcpy(QIFD->uvupvpR2,  Result+24, 9*sizeof(double));

}

/***************************************************************/
/* convert a one-variable quadratic expression into a new form:*/
/*                                                             */
/*  Px^2 + 2Qx + R -> A^2 [ (x+B)^2 + G^2 ]                    */
/*                                                             */
/* note: the quantity G2 is G^2.                               */
/*                                                             */
/* 20121029 note: the test case for why we need the second     */
/* clause for assigning values to B and G2 is the following    */
/* common-edge triangle pair:                                  */
/* V1  = -3.320000e+00 0.000000e+00 -4.166667e+01              */
/* V2  = -2.347595e+00 2.347595e+00 -4.166667e+01              */
/* V3  = -3.320000e+00 0.000000e+00 -4.500000e+01              */
/* V3P = -3.320000e+00 0.000000e+00 -4.166667e+01              */
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

void GetAlphaBetaGamma_FIPPI(TDWorkspaceFIPPI *W, const double *yVector,
                             double *AVector, double *BVector, double *G2Vector)
{
  if (W->WhichCase==TM_COMMONEDGE)
   { 
     double y1 = yVector[0], y12=y1*y1;

     double A2   = W->A2;
     double BP2  = W->BP2;
     double L2   = W->L2;
     double AdBP = W->AdBP;
     double AdL  = W->AdL;
     double BPdL = W->BPdL;

     PQRtoABG2( (BP2 - 2.0*BPdL + L2)*y12, 
                y1*(L2 + BPdL*(y1-1.0) - (AdL + L2-AdBP)*y1),
                L2 - 2*AdL*y1 - 2*L2*y1 + A2*y12 + 2*AdL*y12 + L2*y12,
                AVector+0, BVector+0, G2Vector+0);

     PQRtoABG2( BP2*y12,
                (BPdL + AdBP*y1 - BPdL*y1)*y1,
                L2 + 2*AdL*y1 - 2*L2*y1 + (A2 - 2*AdL + L2)*y12,
                AVector+1, BVector+1, G2Vector+1);

     PQRtoABG2( (A2 + 2*AdBP + BP2)*y12, 
                y1*(AdL*(-1.0 + y1) + BPdL*(-1.0 + y1) - (AdBP + BP2)*y1),
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

     double A2    = W->A2;
     double B2    = W->B2;
     double AP2   = W->AP2;
     double BP2   = W->BP2;
     double AdB   = W->AdB;
     double AdAP  = W->AdAP;
     double AdBP  = W->AdBP;
     double BdAP  = W->BdAP;
     double BdBP  = W->BdBP;
     double APdBP = W->APdBP;

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
/* IntQFP[0]  = \int_0^1 dy [(y+P)^2 + Q^2]^( -3/2)            */
/* IntQFP[1]  = \int_0^1 dy [(y+P)^2 + Q^2]^( -1/2)            */
/* IntQFP[2]  = \int_0^1 dy [(y+P)^2 + Q^2]^(  1/2)            */
/* IntQFP[3]  = \int_0^1 dy [(y+P)^2 + Q^2]^(  2/2)            */
/*                                                             */
/* IntyQFP[0] = \int_0^1 dy y [(y+P)^2 + Q^2]^( -3/2)          */
/* IntyQFP[1] = \int_0^1 dy y [(y+P)^2 + Q^2]^( -1/2)          */
/* IntyQFP[2] = \int_0^1 dy y [(y+P)^2 + Q^2]^(  1/2)          */
/* IntyQFP[3] = \int_0^1 dy y [(y+P)^2 + Q^2]^(  2/2)          */
/***************************************************************/
void GetQFPIntegrals_FIPPI(double P, double Q2,
                           double IntQFP[4], double IntyQFP[4])
{
  double P2     = P*P; 
  double PP1    = P+1.0;

  if ( fabs(Q2) < 1.0e-8 )
   {
     double SignP = P>0.0 ? 1.0 : -1.0;

     // p=-3
     IntQFP[0]  = SignP * (P+0.5) / ( P2*PP1*PP1 );
     IntyQFP[0] = -P*IntQFP[0] + SignP/(P*PP1);

     // p=-1
     IntQFP[1]  = SignP*log( 1.0 + 1.0/P );
     IntyQFP[1] = -P*IntQFP[1] + SignP*1.0;

     // p=+1
     IntQFP[2]  = SignP*(P+0.5);
     IntyQFP[2] = -P*IntQFP[2] + SignP*(P*PP1 + 1.0/3.0);

     // p=+2
     IntQFP[3]  = P*PP1 + 1.0/3.0;
     IntyQFP[3] = -P*IntQFP[3] + P*P2 + (3.0/2.0)*P2 + P + 0.25;
   }  
  else
   { 
     double Q4     = Q2*Q2;
     double S2     = P*P + Q2,  S=sqrt(S2), S3=S*S2, S4=S2*S2;
     double T2     = PP1*PP1 + Q2, T=sqrt(T2), T3=T*T2, T4=T2*T2;
     double LogFac = log( (PP1+T) / (P+S) );
   
     // p=-3
     IntQFP[0]  = ( PP1/T - P/S ) / Q2;
     IntyQFP[0] = -P*IntQFP[0] + 1.0/S - 1.0/T;
   
     // p=-1
     IntQFP[1]  = LogFac;
     IntyQFP[1] = -P*IntQFP[1] + (T-S);
   
     // p=+1
     IntQFP[2]  = (T+P*(T-S) + Q2*LogFac) / 2.0;
     IntyQFP[2] = -P*IntQFP[2] + (T3-S3)/3.0;
   
     // p=+2
     IntQFP[3]  = (T2+S2)/2.0 - 1.0/6.0;
     IntyQFP[3] = -P*IntQFP[3] + (T4-S4)/4.0;
   };
}

/***************************************************************/
/* content of this routine computer-generated by maxima        */
/***************************************************************/
void GetScriptP_CommonEdge(TDWorkspaceFIPPI *W, double y1,
                           double PCross[6][6][4][2], double PXXEE[9][6][4][2])
{
  double *A       = W->A; 
  double *B       = W->B; 
  double *AP      = W->AP; 
  double *BP      = W->BP; 
  double *V0xAP   = W->V0xAP;
  double *V0xBP   = W->V0xBP;
  double *AxV0    = W->AxV0;
  double *BxV0    = W->BxV0;
  double *AxAP    = W->AxAP;
  double *AxBP    = W->AxBP;
  double *BxAP    = W->BxAP;
  double *BxBP    = W->BxBP;

  // slots [0, 1, 2] of the PCross vector all have basically the 
  // same expressions but involving the [0,1,2] component of the
  // vectors on the RHS. (and similarly for slots [3, 4, 5].) so
  for(int Mu=0; Mu<3; Mu++)
   { 
     PCross[Mu][0][1][0] = (BP[Mu]-B[Mu]+AP[Mu])*y1-BP[Mu]+B[Mu];
     PCross[Mu][0][1][1] = B[Mu]*y1;
     PCross[Mu][0][2][0] = -((2.0*BP[Mu]-2.0*B[Mu]+2.0*AP[Mu])*y1-2.0*BP[Mu]+2.0*B[Mu]-AP[Mu] +A[Mu])/2.0;
     PCross[Mu][0][2][1] = -B[Mu]*y1;
     PCross[Mu][0][3][0] = 0.0;
     PCross[Mu][0][3][1] = 0.0;

     PCross[Mu][1][1][0] = (BP[Mu]-B[Mu]-A[Mu])*y1-BP[Mu]+B[Mu];
     PCross[Mu][1][1][1] = -BP[Mu]*y1;
     PCross[Mu][1][2][0] = ((-2.0*BP[Mu]+2.0*B[Mu]+2.0*A[Mu])*y1+2.0*BP[Mu]-2.0*B[Mu]+AP[Mu] -A[Mu])/2.0;
     PCross[Mu][1][2][1] = BP[Mu]*y1;
     PCross[Mu][1][3][0] = 0.0;
     PCross[Mu][1][3][1] = 0.0;

     PCross[Mu][2][1][0] = -B[Mu]*y1-BP[Mu]+B[Mu];
     PCross[Mu][2][1][1] = (BP[Mu]+AP[Mu])*y1;
     PCross[Mu][2][2][0] = -(-2.0*B[Mu]*y1-2.0*BP[Mu]+2.0*B[Mu]-AP[Mu]+A[Mu])/2.0;
     PCross[Mu][2][2][1] = -(BP[Mu]+AP[Mu])*y1;
     PCross[Mu][2][3][0] = 0.0;
     PCross[Mu][2][3][1] = 0.0;

     PCross[Mu][3][1][0] = BP[Mu]*y1-BP[Mu]+B[Mu];
     PCross[Mu][3][1][1] = (-B[Mu]-A[Mu])*y1;
     PCross[Mu][3][2][0] = (-2.0*BP[Mu]*y1+2.0*BP[Mu]-2.0*B[Mu]+AP[Mu]-A[Mu])/2.0;
     PCross[Mu][3][2][1] = (B[Mu]+A[Mu])*y1;
     PCross[Mu][3][3][0] = 0.0;
     PCross[Mu][3][3][1] = 0.0;

     PCross[Mu][4][1][0] = BP[Mu]*y1-BP[Mu]+B[Mu];
     PCross[Mu][4][1][1] = AP[Mu]*y1;
     PCross[Mu][4][2][0] = -(2.0*BP[Mu]*y1-2.0*BP[Mu]+2.0*B[Mu]-AP[Mu]+A[Mu])/2.0;
     PCross[Mu][4][2][1] = -AP[Mu]*y1;
     PCross[Mu][4][3][0] = 0.0;
     PCross[Mu][4][3][1] = 0.0;

     PCross[Mu][5][1][0] = -B[Mu]*y1-BP[Mu]+B[Mu];
     PCross[Mu][5][1][1] = -A[Mu]*y1;
     PCross[Mu][5][2][0] = (2.0*B[Mu]*y1+2.0*BP[Mu]-2.0*B[Mu]+AP[Mu]-A[Mu])/2.0;
     PCross[Mu][5][2][1] = A[Mu]*y1;
     PCross[Mu][5][3][0] = 0.0;
     PCross[Mu][5][3][1] = 0.0;

     PCross[Mu+3][0][1][0] = (+(-2.0*V0xBP[Mu]-2.0*V0xAP[Mu]-2.0*BxV0[Mu]-BxAP[Mu]-AxBP[Mu]-AxAP[Mu])*y1 +2.0*V0xBP[Mu]+2.0*BxV0[Mu]+BxAP[Mu]+AxBP[Mu]) /2.0;
     PCross[Mu+3][0][1][1] = (2.0*BxV0[Mu]+BxAP[Mu])*y1/2.0;
     PCross[Mu+3][0][2][0] = -(+(-2.0*BxBP[Mu]-2.0*BxAP[Mu])*y1*y1 +(-2.0*V0xBP[Mu]-2.0*V0xAP[Mu]-2.0*BxV0[Mu]+4*BxBP[Mu]+2.0*BxAP[Mu])*y1 +2.0*V0xBP[Mu]+V0xAP[Mu]+2.0*BxV0[Mu]-2.0*BxBP[Mu]+AxV0[Mu]) /2.0;
     PCross[Mu+3][0][2][1] = -(((2.0*BxBP[Mu]+2.0*BxAP[Mu])*y1*y1+(2.0*BxV0[Mu]-2.0*BxBP[Mu])*y1))/2.0;
     PCross[Mu+3][0][3][0] = ( +(-6.0*BxBP[Mu]-6.0*BxAP[Mu])*y1*y1 +(12.0*BxBP[Mu]+9*BxAP[Mu]+3.0*AxBP[Mu]+3.0*AxAP[Mu])*y1-6.0*BxBP[Mu] -3.0*BxAP[Mu]-3.0*AxBP[Mu]-2.0*AxAP[Mu]) /6.0;
     PCross[Mu+3][0][3][1] = ((6.0*BxBP[Mu]+6.0*BxAP[Mu])*y1*y1+(-6.0*BxBP[Mu]-3.0*BxAP[Mu])*y1) /6.0;
   
     PCross[Mu+3][1][1][0] = (+(-2.0*V0xBP[Mu]-2.0*BxV0[Mu]-BxAP[Mu]-2.0*AxV0[Mu]-AxBP[Mu]-AxAP[Mu])*y1 +2.0*V0xBP[Mu]+2.0*BxV0[Mu]+BxAP[Mu]+AxBP[Mu]) /2.0;
     PCross[Mu+3][1][1][1] = ((2.0*V0xBP[Mu]+AxBP[Mu])*y1)/2.0;
     PCross[Mu+3][1][2][0] = -(+(-2.0*BxBP[Mu]-2.0*AxBP[Mu])*y1*y1 +(-2.0*V0xBP[Mu]-2.0*BxV0[Mu]+4*BxBP[Mu]-2.0*AxV0[Mu]+2.0*AxBP[Mu])*y1 +2.0*V0xBP[Mu]+V0xAP[Mu]+2.0*BxV0[Mu]-2.0*BxBP[Mu]+AxV0[Mu]) /2.0;
     PCross[Mu+3][1][2][1] = -((2.0*BxBP[Mu]+2.0*AxBP[Mu])*y1*y1+(2.0*V0xBP[Mu]-2.0*BxBP[Mu])*y1)/2.0;
     PCross[Mu+3][1][3][0] = (+(-6.0*BxBP[Mu]-6.0*AxBP[Mu])*y1*y1 +(12.0*BxBP[Mu]+3.0*BxAP[Mu]+9*AxBP[Mu]+3.0*AxAP[Mu])*y1-6.0*BxBP[Mu] -3.0*BxAP[Mu]-3.0*AxBP[Mu]-2.0*AxAP[Mu]) /6.0;
     PCross[Mu+3][1][3][1] = ((6.0*BxBP[Mu]+6.0*AxBP[Mu])*y1*y1+(-6.0*BxBP[Mu]-3.0*AxBP[Mu])*y1)/6.0;
   
     PCross[Mu+3][2][1][0] = -(+(2.0*BxV0[Mu]+BxAP[Mu])*y1-2.0*V0xBP[Mu]-2.0*BxV0[Mu]-BxAP[Mu]-AxBP[Mu]) /2.0;
     PCross[Mu+3][2][1][1] = -(2.0*V0xBP[Mu]+2.0*V0xAP[Mu]+AxBP[Mu]+AxAP[Mu])*y1/2.0;
     PCross[Mu+3][2][2][0] = (+(2.0*BxV0[Mu]-2.0*BxBP[Mu])*y1-2.0*V0xBP[Mu]-V0xAP[Mu]-2.0*BxV0[Mu]+2.0*BxBP[Mu] -AxV0[Mu]) /2.0;
     PCross[Mu+3][2][2][1] = ( (2.0*BxBP[Mu]+2.0*BxAP[Mu])*y1*y1 +(2.0*V0xBP[Mu]+2.0*V0xAP[Mu]-2.0*BxBP[Mu]-2.0*BxAP[Mu])*y1)/2.0;
     PCross[Mu+3][2][3][0] = -(+(-6.0*BxBP[Mu]-3.0*BxAP[Mu])*y1+6.0*BxBP[Mu]+3.0*BxAP[Mu]+3.0*AxBP[Mu] +2.0*AxAP[Mu]) /6.0;
     PCross[Mu+3][2][3][1] = -(((6.0*BxBP[Mu]+6.0*BxAP[Mu])*y1*y1 +(-6.0*BxBP[Mu]-6.0*BxAP[Mu]-3.0*AxBP[Mu]-3.0*AxAP[Mu])*y1))/6.0;
   
     PCross[Mu+3][3][1][0] = -(+(2.0*V0xBP[Mu]+AxBP[Mu])*y1-2.0*V0xBP[Mu]-2.0*BxV0[Mu]-BxAP[Mu]-AxBP[Mu]) /2.0;
     PCross[Mu+3][3][1][1] = -(2.0*BxV0[Mu]+BxAP[Mu]+2.0*AxV0[Mu]+AxAP[Mu])*y1/2.0;
     PCross[Mu+3][3][2][0] = (+(2.0*V0xBP[Mu]-2.0*BxBP[Mu])*y1-2.0*V0xBP[Mu]-V0xAP[Mu]-2.0*BxV0[Mu] +2.0*BxBP[Mu]-AxV0[Mu]) /2.0;
     PCross[Mu+3][3][2][1] = ((2.0*BxBP[Mu]+2.0*AxBP[Mu])*y1*y1 +(2.0*BxV0[Mu]-2.0*BxBP[Mu]+2.0*AxV0[Mu]-2.0*AxBP[Mu])*y1)/2.0; 
     PCross[Mu+3][3][3][0] = -(+(-6.0*BxBP[Mu]-3.0*AxBP[Mu])*y1+6.0*BxBP[Mu]+3.0*BxAP[Mu]+3.0*AxBP[Mu] +2.0*AxAP[Mu]) /6.0;
     PCross[Mu+3][3][3][1] = -((6.0*BxBP[Mu]+6.0*AxBP[Mu])*y1*y1 +(-6.0*BxBP[Mu]-3.0*BxAP[Mu]-6.0*AxBP[Mu]-3.0*AxAP[Mu])*y1)/6.0;
   
     PCross[Mu+3][4][1][0] = -(+(2.0*V0xBP[Mu]+AxBP[Mu])*y1-2.0*V0xBP[Mu] -2.0*BxV0[Mu]-BxAP[Mu]-AxBP[Mu]) /2.0;
     PCross[Mu+3][4][1][1] = -((2.0*V0xAP[Mu]+AxAP[Mu])*y1)/2.0;
     PCross[Mu+3][4][2][0] = (+(2.0*V0xBP[Mu]-2.0*BxBP[Mu])*y1 -2.0*V0xBP[Mu]-V0xAP[Mu]-2.0*BxV0[Mu] +2.0*BxBP[Mu]-AxV0[Mu]) /2.0;
     PCross[Mu+3][4][2][1] = (2.0*V0xAP[Mu]-2.0*BxAP[Mu])*y1/2.0;
     PCross[Mu+3][4][3][0] = (+(6.0*BxBP[Mu]+3.0*AxBP[Mu])*y1-6.0*BxBP[Mu] -3.0*BxAP[Mu]-3.0*AxBP[Mu]-2.0*AxAP[Mu]) /6.0;
     PCross[Mu+3][4][3][1] = (6.0*BxAP[Mu]+3.0*AxAP[Mu])*y1/6.0;
   
     PCross[Mu+3][5][1][0] = -(+(2.0*BxV0[Mu]+BxAP[Mu])*y1-2.0*V0xBP[Mu] -2.0*BxV0[Mu]-BxAP[Mu]-AxBP[Mu]) /2.0;
     PCross[Mu+3][5][1][1] = -(2.0*AxV0[Mu]+AxAP[Mu])*y1/2.0;
     PCross[Mu+3][5][2][0] = (+(2.0*BxV0[Mu]-2.0*BxBP[Mu])*y1 -2.0*V0xBP[Mu]-V0xAP[Mu]-2.0*BxV0[Mu] +2.0*BxBP[Mu]-AxV0[Mu]) /2.0;
     PCross[Mu+3][5][2][1] = (2.0*AxV0[Mu]-2.0*AxBP[Mu])*y1/2.0;
     PCross[Mu+3][5][3][0] = (+(6.0*BxBP[Mu]+3.0*BxAP[Mu])*y1-6.0*BxBP[Mu] -3.0*BxAP[Mu]-3.0*AxBP[Mu]-2.0*AxAP[Mu]) /6.0;
     PCross[Mu+3][5][3][1] = (6.0*AxBP[Mu]+3.0*AxAP[Mu])*y1/6.0;
   };

  // note that only the nonzero entries are filled in; 
  // somebody before us needs to make sure that everything
  // else is zeroed out 
  PXXEE[0][0][0][0] = 1.0;
  PXXEE[0][0][1][0] = -1.0;
  PXXEE[0][1][0][0] = 1.0;
  PXXEE[0][1][1][0] = -1.0;
  PXXEE[0][2][0][0] = 1.0;
  PXXEE[0][2][1][0] = -1.0;
  PXXEE[0][3][0][0] = 1.0;
  PXXEE[0][3][1][0] = -1.0;
  PXXEE[0][4][0][0] = 1.0;
  PXXEE[0][4][1][0] = -1.0;
  PXXEE[0][5][0][0] = 1.0;
  PXXEE[0][5][1][0] = -1.0;
  PXXEE[1][0][0][0] = 1.0/2.0;
  PXXEE[1][0][1][0] = -y1;
  PXXEE[1][0][2][0] = (2.0*y1-1.0)/2.0;
  PXXEE[1][1][0][0] = 1.0/2.0;
  PXXEE[1][1][2][0] = -1.0/2.0;
  PXXEE[1][2][0][0] = 1.0/2.0;
  PXXEE[1][2][1][1] = -y1;
  PXXEE[1][2][2][0] = -1.0/2.0;
  PXXEE[1][2][2][1] = y1;
  PXXEE[1][3][0][0] = 1.0/2.0;
  PXXEE[1][3][2][0] = -1.0/2.0;
  PXXEE[1][4][0][0] = 1.0/2.0;
  PXXEE[1][4][1][1] = -y1;
  PXXEE[1][4][2][0] = -1.0/2.0;
  PXXEE[1][4][2][1] = y1;
  PXXEE[1][5][0][0] = 1.0/2.0;
  PXXEE[1][5][2][0] = -1.0/2.0;
  PXXEE[2][0][1][0] = -y1+1.0;
  PXXEE[2][0][2][0] = y1-1.0;
  PXXEE[2][1][1][0] = -y1+1.0;
  PXXEE[2][1][1][1] = y1;
  PXXEE[2][1][2][0] = +y1-1.0;
  PXXEE[2][1][2][1] = -y1;
  PXXEE[2][2][1][0] = 1.0;
  PXXEE[2][2][1][1] = -y1;
  PXXEE[2][2][2][0] = -1.0;
  PXXEE[2][2][2][1] = y1;
  PXXEE[2][3][1][0] = -y1+1.0;
  PXXEE[2][3][2][0] = y1-1.0;
  PXXEE[2][4][1][0] = -y1+1.0;
  PXXEE[2][4][2][0] = y1-1.0;
  PXXEE[2][5][1][0] = 1.0;
  PXXEE[2][5][2][0] = -1.0;
  PXXEE[3][0][0][0] = 1.0/2.0;
  PXXEE[3][0][2][0] = -1.0/2.0;
  PXXEE[3][1][0][0] = 1.0/2.0;
  PXXEE[3][1][1][0] = -y1;
  PXXEE[3][1][2][0] = (2.0*y1-1.0)/2.0;
  PXXEE[3][2][0][0] = 1.0/2.0;
  PXXEE[3][2][2][0] = -1.0/2.0;
  PXXEE[3][3][0][0] = 1.0/2.0;
  PXXEE[3][3][1][1] = -y1;
  PXXEE[3][3][2][0] = -1.0/2.0;
  PXXEE[3][3][2][1] = y1;
  PXXEE[3][4][0][0] = 1.0/2.0;
  PXXEE[3][4][2][0] = -1.0/2.0;
  PXXEE[3][5][0][0] = 1.0/2.0;
  PXXEE[3][5][1][1] = -y1;
  PXXEE[3][5][2][0] = -1.0/2.0;
  PXXEE[3][5][2][1] = y1;
  PXXEE[4][0][0][0] = 1.0/3.0;
  PXXEE[4][0][1][0] = -y1/2.0;
  PXXEE[4][0][3][0] = (3.0*y1-2.0)/6.0;
  PXXEE[4][1][0][0] = 1.0/3.0;
  PXXEE[4][1][1][0] = -y1/2.0;
  PXXEE[4][1][3][0] = (3.0*y1-2.0)/6.0;
  PXXEE[4][2][0][0] = 1.0/3.0;
  PXXEE[4][2][1][1] = -y1/2.0;
  PXXEE[4][2][3][0] = -1.0/3.0;
  PXXEE[4][2][3][1] = y1/2.0;
  PXXEE[4][3][0][0] = 1.0/3.0;
  PXXEE[4][3][1][1] = -y1/2.0;
  PXXEE[4][3][3][0] = -1.0/3.0;
  PXXEE[4][3][3][1] = y1/2.0;
  PXXEE[4][4][0][0] = 1.0/3.0;
  PXXEE[4][4][1][1] = -y1/2.0;
  PXXEE[4][4][3][0] = -1.0/3.0;
  PXXEE[4][4][3][1] = y1/2.0;
  PXXEE[4][5][0][0] = 1.0/3.0;
  PXXEE[4][5][1][1] = -y1/2.0;
  PXXEE[4][5][3][0] = -1.0/3.0;
  PXXEE[4][5][3][1] = y1/2.0;
  PXXEE[5][0][1][0] = -(y1-1.0)/2.0;
  PXXEE[5][0][3][0] = (y1-1.0)/2.0;
  PXXEE[5][1][1][0] = (-y1+1.0)/2.0;
  PXXEE[5][1][1][1] = y1/2.0;
  PXXEE[5][1][2][0] = +y1*y1-y1;
  PXXEE[5][1][2][1] = -y1*y1;
  PXXEE[5][1][3][0] = (-2.0*y1*y1+3.0*y1-1.0)/2.0;
  PXXEE[5][1][3][1] = (2.0*y1*y1-y1)/2.0;
  PXXEE[5][2][1][0] = -(-1.0)/2.0;
  PXXEE[5][2][1][1] = -y1/2.0;
  PXXEE[5][2][3][0] = (-1.0)/2.0;
  PXXEE[5][2][3][1] = y1/2.0;
  PXXEE[5][3][1][0] = -(y1-1.0)/2.0;
  PXXEE[5][3][2][1] = (y1*y1-y1);
  PXXEE[5][3][3][0] = -(-y1+1.0)/2.0;
  PXXEE[5][3][3][1] = -(y1*y1-y1);
  PXXEE[5][4][1][0] = -(y1-1.0)/2.0;
  PXXEE[5][4][3][0] = (y1-1.0)/2.0;
  PXXEE[5][5][1][0] = 1.0/2.0;
  PXXEE[5][5][2][1] = -y1;
  PXXEE[5][5][3][0] = -1.0/2.0;
  PXXEE[5][5][3][1] = y1;
  PXXEE[6][0][1][0] = -y1+1.0;
  PXXEE[6][0][1][1] = y1;
  PXXEE[6][0][2][0] = +y1-1.0;
  PXXEE[6][0][2][1] = -y1;
  PXXEE[6][1][1][0] = -y1+1.0;
  PXXEE[6][1][2][0] = y1-1.0;
  PXXEE[6][2][1][0] = -y1+1.0;
  PXXEE[6][2][2][0] = y1-1.0;
  PXXEE[6][3][1][0] = 1.0;
  PXXEE[6][3][1][1] = -y1;
  PXXEE[6][3][2][0] = -1.0;
  PXXEE[6][3][2][1] = y1;
  PXXEE[6][4][1][0] = 1.0;
  PXXEE[6][4][2][0] = -1.0;
  PXXEE[6][5][1][0] = -y1+1.0;
  PXXEE[6][5][2][0] = y1-1.0;
  PXXEE[7][0][1][0] = (-y1+1.0)/2.0;
  PXXEE[7][0][1][1] = y1/2.0;
  PXXEE[7][0][2][0] = y1*y1-y1;
  PXXEE[7][0][2][1] = -y1*y1;
  PXXEE[7][0][3][0] = (-2.0*y1*y1+3.0*y1-1.0)/2.0;
  PXXEE[7][0][3][1] = (2.0*y1*y1-y1)/2.0;
  PXXEE[7][1][1][0] = -(y1-1.0)/2.0;
  PXXEE[7][1][3][0] = (y1-1.0)/2.0;
  PXXEE[7][2][1][0] = -(y1-1.0)/2.0;
  PXXEE[7][2][2][1] = (y1*y1-y1);
  PXXEE[7][2][3][0] = -(-y1+1.0)/2.0;
  PXXEE[7][2][3][1] = -(y1*y1-y1);
  PXXEE[7][3][1][0] = -(-1.0)/2.0;
  PXXEE[7][3][1][1] = -y1/2.0;
  PXXEE[7][3][3][0] = -1.0/2.0;
  PXXEE[7][3][3][1] = y1/2.0;
  PXXEE[7][4][1][0] = 1.0/2.0;
  PXXEE[7][4][2][1] = -y1;
  PXXEE[7][4][3][0] = -1.0/2.0;
  PXXEE[7][4][3][1] = y1;
  PXXEE[7][5][1][0] = -(y1-1.0)/2.0;
  PXXEE[7][5][3][0] = (y1-1.0)/2.0;
  PXXEE[8][0][2][0] = +y1*y1-2.0*y1+1.0;
  PXXEE[8][0][2][1] = (-y1*y1+y1);
  PXXEE[8][0][3][0] = -y1*y1+2.0*y1-1.0;
  PXXEE[8][0][3][1] = (y1*y1-y1);
  PXXEE[8][1][2][0] = +y1*y1-2.0*y1+1.0;
  PXXEE[8][1][2][1] = (-y1*y1+y1);
  PXXEE[8][1][3][0] = -y1*y1+2.0*y1-1.0;
  PXXEE[8][1][3][1] = (y1*y1-y1);
  PXXEE[8][2][2][0] = -y1+1.0;
  PXXEE[8][2][2][1] = (y1*y1-y1);
  PXXEE[8][2][3][0] = y1-1.0;
  PXXEE[8][2][3][1] = (-y1*y1+y1);
  PXXEE[8][3][2][0] = -y1+1.0;
  PXXEE[8][3][2][1] = (y1*y1-y1);
  PXXEE[8][3][3][0] = +y1-1.0;
  PXXEE[8][3][3][1] = (-y1*y1+y1);
  PXXEE[8][4][2][0] = -y1+1.0;
  PXXEE[8][4][3][0] = y1-1.0;
  PXXEE[8][5][2][0] = -y1+1.0;
  PXXEE[8][5][3][0] = y1-1.0;

}

/***************************************************************/
/* content of this routine computer-generated by maxima        */
/***************************************************************/
void GetScriptP_CommonVertex(TDWorkspaceFIPPI *W, double y1, double y2, 
                             double PCross[6][6][4][2], double PXXEE[9][6][4][2])
{
  double *A       = W->A; 
  double *B       = W->B; 
  double *AP      = W->AP; 
  double *BP      = W->BP; 
  double *V0xAP   = W->V0xAP;
  double *V0xBP   = W->V0xBP;
  double *AxV0    = W->AxV0;
  double *BxV0    = W->BxV0;
  double *AxAP    = W->AxAP;
  double *AxBP    = W->AxBP;
  double *BxAP    = W->BxAP;
  double *BxBP    = W->BxBP;

  for(int Mu=0; Mu<3; Mu++)
   { PCross[Mu][0][1][0] = -AP[Mu]*y2+B[Mu]*y1+A[Mu];
     PCross[Mu][0][1][1] = -BP[Mu]*y2;
     PCross[Mu][0][2][0] = 0.0;
     PCross[Mu][0][2][1] = 0.0;
     PCross[Mu][1][1][0] = A[Mu]*y2-BP[Mu]*y1-AP[Mu];
     PCross[Mu][1][1][1] = B[Mu]*y2;
     PCross[Mu][1][2][0] = 0.0;
     PCross[Mu][1][2][1] = 0.0;

     PCross[Mu+3][0][1][0] = V0xAP[Mu]*y2+BxV0[Mu]*y1+AxV0[Mu];
     PCross[Mu+3][0][1][1] = V0xBP[Mu]*y2;
     PCross[Mu+3][0][2][0] = (BxAP[Mu]*y1+AxAP[Mu])*y2;
     PCross[Mu+3][0][2][1] = (BxBP[Mu]*y1+AxBP[Mu])*y2;
     PCross[Mu+3][1][1][0] = AxV0[Mu]*y2+V0xBP[Mu]*y1+V0xAP[Mu];
     PCross[Mu+3][1][1][1] = BxV0[Mu]*y2;
     PCross[Mu+3][1][2][0] = (AxBP[Mu]*y1+AxAP[Mu])*y2;
     PCross[Mu+3][1][2][1] = (BxBP[Mu]*y1+BxAP[Mu])*y2;
   };


  // note that only the nonzero entries are filled in; 
  // somebody before us needs to make sure that everything
  // else is zeroed out 
  PXXEE[0][0][0][0] = 1.0;
  PXXEE[0][1][0][0] = 1.0;

  PXXEE[1][0][1][0] = y2;
  PXXEE[1][1][1][0] = 1.0;

  PXXEE[2][0][1][1] = y2;
  PXXEE[2][1][1][0] = y1;

  PXXEE[3][0][1][0] = 1.0;
  PXXEE[3][1][1][0] = y2;

  PXXEE[4][0][2][0] = y2;
  PXXEE[4][1][2][0] = y2;

  PXXEE[5][0][2][1] = y2;
  PXXEE[5][1][2][0] = y1*y2;

  PXXEE[6][0][1][0] = y1;
  PXXEE[6][1][1][1] = y2;

  PXXEE[7][0][2][0] = y1*y2;
  PXXEE[7][1][2][1] = y2;

  PXXEE[8][0][2][1] = y1*y2;
  PXXEE[8][1][2][1] = y1*y2;
}

} // namespace scuff
