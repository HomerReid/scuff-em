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
 * HighKTaylorDuffy.cc:  an implementation of the high-k (short-wavelength)
 *                       limit of the Taylor-Duffy method for evaluating 
 *                       singular panel-panel integrals
 *
 * homer reid             9/2012
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

/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
//int iOnly=-1, AlphaOnly=-1;
//if ( (iOnly==-1 || i==iOnly) && (AlphaOnly==-1 || Alpha==AlphaOnly) )
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/

namespace scuff {

/***************************************************************/
/***************************************************************/
/***************************************************************/
int AICheck(const char *FileName, int LineNum, 
            double *I, double *E, double AbsTol, double RelTol, int Length)
{
  int ViolationFound=0;

  for (int n=0; n<Length; n++)
   //if ( (E[n] > AbsTol) && (E[n] > RelTol*fabs(I[n]) ) )
   if ( (E[n] > 0.1*fabs(I[n])) )
    { ViolationFound=1;
      Log("AI(%s:%i): (I,E)[%i]=(%.1e,%.1e)",FileName,LineNum,n,I[n],E[n]);
    };

  return ViolationFound;
}

#define DEFABSTOL 0.0
#define DEFRELTOL 1.0e-6
#define MAXFEVALS 10000
#define INTERVALS (MAXFEVALS/15)

#define II cdouble(0,1)

static double FactorialTable[7]={1.0, 1.0, 2.0, 6.0, 24.0, 120.0, 620.0};

void SipAlpha_One(const double *xVec, TMWorkspace *TMW, int WhichCase,
                  int *AlphaMin, int *AlphaMax, cdouble S[7][5][2]);
void SipAlpha_Dot(const double *xVec, TMWorkspace *TMW, int WhichCase,
                  int *AlphaMin, int *AlphaMax, cdouble S[7][5][2]);
void SipAlpha_DotPlus(const double *xVec, TMWorkspace *TMW, int WhichCase,
                      int *AlphaMin, int *AlphaMax, cdouble S[7][5][2]);
void SipAlpha_Cross(const double *xVec, TMWorkspace *TMW, int WhichCase,
                    int *AlphaMin, int *AlphaMax, cdouble S[7][5][2]);

/***************************************************************/
/* This routine evaluates the \mathcal{J} and \mathcal{K}      */
/* integrals defined in the memo for the given values of       */
/* P and Q and for Alpha=1..6. (The Alpha==0 slot in the       */
/* JVector and KVector arrays is not referenced.)              */
/*                                                             */
/* Definitions:                                                */
/*                                                             */
/* J_\alpha = \int_0^1 dx [(x+P)^2 + Q^2]^(-\Alpha/2)          */
/*                                                             */
/* K_\alpha = \int_0^1 dx x*[(x+P)^2 + Q^2]^(-\Alpha/2)        */
/*                                                             */
/***************************************************************/
static void GetJKVectors(double P, double Q, double *JVector, double *KVector)
{ 
  double P2  = P*P;
  double PP1 = P+1.0, PP12=PP1*PP1;
  double Q2  = Q*Q, Q3=Q2*Q, Q4=Q3*Q;
  double S2  = P*P + Q*Q, S=sqrt(S2), S3=S*S2, S4=S3*S;
  double T2  = PP1*PP1 + Q*Q, T=sqrt(T2), T3=T*T2, T4=T3*T;

  JVector[1] = log( (PP1+T) / (P+S) );
  JVector[2] = ( atan(PP1/Q) - atan(P/Q) ) / Q;
  JVector[3] = ( PP1/T - P/S ) / (Q2);
  JVector[4] = ((Q2-P*PP1)/(S2*T2) + JVector[2]) / (2.0*Q2);
  JVector[5] = (PP1*(2.0*PP12 + 3.0*Q2)*S3 - P*(2.0*P2 + 3.0*Q2)*T3)/(3.0*Q4*S3*T3);
  JVector[6] = (3.0*PP1/T2 + 2.0*Q2*PP1/T4 - 3.0*P/S2 - 2.0*Q2*P/S4 + 3*JVector[2])/(8.0*Q4);

  KVector[1] = T - S    -  P*JVector[1];
  KVector[2] = log(T/S) -  P*JVector[2];
  KVector[3] = 1.0/S - 1.0/T - P*JVector[3];
  KVector[4] = 1.0/(2.0*S2) - 1.0/(2.0*T2) - P*JVector[4];
  KVector[5] = 1.0/(3.0*S3) - 1.0/(3.0*T3) - P*JVector[5];
  KVector[6] = 1.0/(4.0*S4) - 1.0/(4.0*T4) - P*JVector[6];

}

/***************************************************************/
/* convert a one-variable quadratic expression into a new form:*/
/*  Ax^2 + 2Bx + C -> M^2 [ (x+P)^2 + Q^2 ]                    */
/***************************************************************/
void ABCtoMPQ(double A, double B, double C, double *M, double *P, double *Q)
{ 
  *M=sqrt(A);
  *P=B/A;
  *Q=sqrt( C/A - (*P)*(*P) );
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
static void GetMPQ(int WhichCase, TMWorkspace *TMW, const double *xVector,
                   double *MVector, double *PVector, double *QVector)
{
  if (WhichCase==TM_COMMONTRIANGLE) 
   { 
     double A2  = TMW->A2; 
     double AdB = TMW->AdB;
     double B2  = TMW->B2; 

     double C2 = A2 + B2 + 2.0*AdB;
     double CdB = AdB + B2;

     MVector[1] = sqrt(B2);
     MVector[2] = sqrt(C2);
     MVector[3] = sqrt(A2);

     PVector[1] = AdB/B2;
     PVector[2] = -CdB/C2;
     PVector[3] = AdB/A2;

     QVector[1] = sqrt( A2/B2 - (AdB*AdB)/(B2*B2) );
     QVector[2] = sqrt( B2/C2 - (CdB*CdB)/(C2*C2) );
     QVector[3] = sqrt( B2/A2 - (AdB*AdB)/(A2*A2) );
   }
  else if (WhichCase==TM_COMMONEDGE) 
   { 
     double x1 = xVector[0], x12=x1*x1;

     double A2   = TMW->A2;
     double BP2  = TMW->BP2;
     double L2   = TMW->L2;
     double AdBP = TMW->AdBP;
     double AdL  = TMW->AdL;
     double BPdL = TMW->BPdL;

     ABCtoMPQ( (BP2 - 2.0*BPdL + L2)*x12, 
               x1*(L2 + BPdL*(x1-1.0) - (AdL + L2-AdBP)*x1),
               L2 - 2*AdL*x1 - 2*L2*x1 + A2*x12 + 2*AdL*x12 + L2*x12,
               MVector+1, PVector+1, QVector+1);

     ABCtoMPQ( BP2*x12,
               (BPdL + AdBP*x1 - BPdL*x1)*x1,
               L2 + 2*AdL*x1 - 2*L2*x1 + (A2 - 2*AdL + L2)*x12,
               MVector+2, PVector+2, QVector+2);

     ABCtoMPQ( (A2 + 2*AdBP + BP2)*x12, 
               x1*(AdL*(-1 + x1) + BPdL*(-1 + x1) - (AdBP + BP2)*x1),
               L2 + 2*BPdL*x1 - 2*L2*x1 + BP2*x12 - 2*BPdL*x12 + L2*x12,
               MVector+3, PVector+3, QVector+3);

     ABCtoMPQ( (A2 + 2*AdBP - 2*AdL + BP2 - 2*BPdL + L2)*x12,
               x1*(AdL - L2 - (AdBP + BP2)*x1 + BPdL*(1 + x1)),
               L2 - 2*BPdL*x1 + BP2*x12,
               MVector+4, PVector+4, QVector+4);

     ABCtoMPQ( A2*x12, 
               (AdBP*x1-AdL)*x1, 
               L2 - 2*BPdL*x1 + BP2*x12,
               MVector+5, PVector+5, QVector+5);

     ABCtoMPQ( A2*x12,
               x1*(AdL + AdBP*x1 - AdL*x1),
               L2 + 2*BPdL*x1 - 2*L2*x1 + BP2*x12 - 2*BPdL*x12 + L2*x12,
               MVector+6, PVector+6, QVector+6);
   }
  else // (WhichCase==TM_COMMONVERTEX) 
   {
     double x1 = xVector[0];
     double x2 = xVector[1], x22=x2*x2;

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

     ABCtoMPQ( BP2*x22, 
               -AdBP*x2 -BdBP*x1*x2 +APdBP*x22,
               A2 + B2*x1*x1 + AP2*x2*x2 + 2.0*AdB*x1 - 2.0*AdAP*x2 -2.0*BdAP*x1*x2,
               MVector+1, PVector+1, QVector+1);

     ABCtoMPQ( B2*x22,
               AdB*x22 - BdAP*x2 - BdBP*x1*x2,
               A2*x2*x2 + AP2 + BP2*x1*x1 - 2*AdAP*x2 - 2*AdBP*x2*x1 + 2*APdBP*x1,
               MVector+2, PVector+2, QVector+2);

   };

}

/***************************************************************/
/* High-k Taylor-Duffy computation for the common-triangle case*/
/***************************************************************/
cdouble HKTD_CommonTriangle(TMWorkspace *TMW)
{
  cdouble k = TMW->GParam;

  /*--------------------------------------------------------------*/
  /*- prefetch the S_{i,\alpha,p} coefficients -------------------*/
  /*--------------------------------------------------------------*/
  int AlphaMin, AlphaMax;
  cdouble S[7][5][2];
  if (TMW->WhichH == TM_ONE )
   SipAlpha_One(0, TMW, TM_COMMONTRIANGLE, &AlphaMin, &AlphaMax, S);
  else if (TMW->WhichH == TM_DOT )
   SipAlpha_Dot(0, TMW, TM_COMMONTRIANGLE, &AlphaMin, &AlphaMax, S);
  else if (TMW->WhichH == TM_DOTPLUS )
   SipAlpha_DotPlus(0, TMW, TM_COMMONTRIANGLE, &AlphaMin, &AlphaMax, S);

  /*--------------------------------------------------------------*/
  /*- prefetch the M_i, P_i, Q_i coefficients --------------------*/
  /*--------------------------------------------------------------*/
  double MVector[4], PVector[4], QVector[4];
  GetMPQ(TM_COMMONTRIANGLE, TMW, 0, MVector, PVector, QVector);

  /*--------------------------------------------------------------*/
  /*- prefetch the values of the J, K integrals ------------------*/
  /*--------------------------------------------------------------*/
  double JVector[4][7], KVector[4][7];
  for(int i=1; i<4; i++)
   GetJKVectors(PVector[i], QVector[i], JVector[i], KVector[i]);

  /*--------------------------------------------------------------*/
  /*- initialize the constant prefactors for the summation loop  -*/
  /*--------------------------------------------------------------*/
  cdouble MinusikAlphaTable[7];
  double MiAlphaTable[4][7];
  for(int Alpha=0; Alpha<7; Alpha++)
   { MinusikAlphaTable[Alpha] = (Alpha==0) ? 1.0: MinusikAlphaTable[Alpha-1]*(-II*k);
     MiAlphaTable[1][Alpha]   = (Alpha==0) ? 1.0: MiAlphaTable[1][Alpha-1]*MVector[1];
     MiAlphaTable[2][Alpha]   = (Alpha==0) ? 1.0: MiAlphaTable[2][Alpha-1]*MVector[2];
     MiAlphaTable[3][Alpha]   = (Alpha==0) ? 1.0: MiAlphaTable[3][Alpha-1]*MVector[3];
   };


  /*--------------------------------------------------------------*/
  /*- sum terms --------------------------------------------------*/
  /*--------------------------------------------------------------*/
  cdouble Sum=0.0;
  for(int Alpha=AlphaMin; Alpha<=AlphaMax; Alpha++)
   { 
     for(int i=1; i<4; i++)
      Sum += FactorialTable[Alpha]*(  S[i][Alpha][0]*JVector[i][Alpha+2] 
                                     +S[i][Alpha][1]*KVector[i][Alpha+2] 
                                   ) / (MinusikAlphaTable[Alpha+1]*MiAlphaTable[i][Alpha+2]);

   };

  return Sum;
 
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void HKTD_CommonEdgeIntegrand(unsigned ndim, const double *x, void *parms,
                              unsigned nfun, double *f)
{
  (void) ndim; (void) nfun; // unused
  TMWorkspace *TMW = (TMWorkspace *)parms;
  TMW->nCalls++;
  cdouble k = TMW->GParam;

  /*--------------------------------------------------------------*/
  /*- prefetch the S_{i,\alpha,p} coefficients -------------------*/
  /*--------------------------------------------------------------*/
  int AlphaMin, AlphaMax;
  cdouble S[7][5][2];
  if (TMW->WhichH == TM_ONE )
   SipAlpha_One(x, TMW, TM_COMMONEDGE, &AlphaMin, &AlphaMax, S);
  else if ( TMW->WhichH == TM_DOT )
   SipAlpha_Dot(x, TMW, TM_COMMONEDGE, &AlphaMin, &AlphaMax, S);
  else if ( TMW->WhichH == TM_DOTPLUS )
   SipAlpha_DotPlus(x, TMW, TM_COMMONEDGE, &AlphaMin, &AlphaMax, S);
  else if ( TMW->WhichH == TM_CROSS )
   SipAlpha_Cross(x, TMW, TM_COMMONEDGE, &AlphaMin, &AlphaMax, S);

  /*--------------------------------------------------------------*/
  /*- prefetch the M_i, P_i, Q_i coefficients --------------------*/
  /*--------------------------------------------------------------*/
  double MVector[7], PVector[7], QVector[7];
  GetMPQ(TM_COMMONEDGE, TMW, x, MVector, PVector, QVector);

  /*--------------------------------------------------------------*/
  /*- prefetch the values of the J, K integrals ------------------*/
  /*--------------------------------------------------------------*/
  double JVector[7][7], KVector[7][7];
  for(int i=1; i<=6; i++)
   GetJKVectors(PVector[i], QVector[i], JVector[i], KVector[i]);

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  cdouble MinusikAlphaTable[7];
  double MiAlphaTable[7][7];
  for(int Alpha=0; Alpha<7; Alpha++)
   { MinusikAlphaTable[Alpha] = (Alpha==0) ? 1.0: MinusikAlphaTable[Alpha-1]*(-II*k);
     MiAlphaTable[1][Alpha]   = (Alpha==0) ? 1.0: MiAlphaTable[1][Alpha-1]*MVector[1];
     MiAlphaTable[2][Alpha]   = (Alpha==0) ? 1.0: MiAlphaTable[2][Alpha-1]*MVector[2];
     MiAlphaTable[3][Alpha]   = (Alpha==0) ? 1.0: MiAlphaTable[3][Alpha-1]*MVector[3];
     MiAlphaTable[4][Alpha]   = (Alpha==0) ? 1.0: MiAlphaTable[4][Alpha-1]*MVector[4];
     MiAlphaTable[5][Alpha]   = (Alpha==0) ? 1.0: MiAlphaTable[5][Alpha-1]*MVector[5];
     MiAlphaTable[6][Alpha]   = (Alpha==0) ? 1.0: MiAlphaTable[6][Alpha-1]*MVector[6];
   };

  /*--------------------------------------------------------------*/
  /*- sum terms --------------------------------------------------*/
  /*--------------------------------------------------------------*/
  cdouble Sum=0.0;
  if (TMW->WhichG==TM_EIKR_OVER_R)
   { for(int Alpha=AlphaMin; Alpha<=AlphaMax; Alpha++)
      for(int i=1; i<=6; i++)
       Sum += FactorialTable[Alpha+1]
              *(S[i][Alpha][0]*JVector[i][Alpha+3] + S[i][Alpha][1]*KVector[i][Alpha+3])
              / ( MinusikAlphaTable[Alpha+2] * MiAlphaTable[i][Alpha+3]);

   }
  else if (TMW->WhichG==TM_GRADEIKR_OVER_R)
   { for(int Alpha=AlphaMin; Alpha<=AlphaMax; Alpha++)
      for(int i=1; i<=6; i++)
       Sum -= (Alpha+1)*FactorialTable[Alpha-1]
              *(S[i][Alpha][0]*JVector[i][Alpha+3] + S[i][Alpha][1]*KVector[i][Alpha+3])
              / ( MinusikAlphaTable[Alpha] * MiAlphaTable[i][Alpha+3]);

   };

  *((cdouble *)f) = x[0]*Sum;
 
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void HKTD_CommonVertexIntegrand(unsigned ndim, const double *x, void *parms,
                                unsigned nfun, double *f)
{
  (void) ndim; (void) nfun; // unused
  TMWorkspace *TMW = (TMWorkspace *)parms;
  TMW->nCalls++;
  cdouble k = TMW->GParam;

  /*--------------------------------------------------------------*/
  /*- prefetch the S_{i,\alpha,p} coefficients -------------------*/
  /*--------------------------------------------------------------*/
  int AlphaMin, AlphaMax;
  cdouble S[7][5][2];
  if (TMW->WhichH == TM_ONE )
   SipAlpha_One(x, TMW, TM_COMMONVERTEX, &AlphaMin, &AlphaMax, S);
  else if (TMW->WhichH == TM_DOT )
   SipAlpha_Dot(x, TMW, TM_COMMONVERTEX, &AlphaMin, &AlphaMax, S);
  else if (TMW->WhichH == TM_DOTPLUS )
   SipAlpha_DotPlus(x, TMW, TM_COMMONVERTEX, &AlphaMin, &AlphaMax, S);
  else if (TMW->WhichH == TM_CROSS )
   SipAlpha_Cross(x, TMW, TM_COMMONVERTEX, &AlphaMin, &AlphaMax, S);

  /*--------------------------------------------------------------*/
  /*- prefetch the M_i, P_i, Q_i coefficients --------------------*/
  /*--------------------------------------------------------------*/
  double MVector[3], PVector[3], QVector[3];
  GetMPQ(TM_COMMONVERTEX, TMW, x, MVector, PVector, QVector);

  /*--------------------------------------------------------------*/
  /*- prefetch the values of the J, K integrals ------------------*/
  /*--------------------------------------------------------------*/
  double JVector[3][7], KVector[3][7];
  for(int i=1; i<=2; i++)
   GetJKVectors(PVector[i], QVector[i], JVector[i], KVector[i]);

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  cdouble MinusikAlphaTable[7];
  double MiAlphaTable[3][7];
  for(int Alpha=0; Alpha<7; Alpha++)
   { MinusikAlphaTable[Alpha] = (Alpha==0) ? 1.0: MinusikAlphaTable[Alpha-1]*(-II*k);
     MiAlphaTable[1][Alpha]   = (Alpha==0) ? 1.0: MiAlphaTable[1][Alpha-1]*MVector[1];
     MiAlphaTable[2][Alpha]   = (Alpha==0) ? 1.0: MiAlphaTable[2][Alpha-1]*MVector[2];
   };

  /*--------------------------------------------------------------*/
  /*- sum terms --------------------------------------------------*/
  /*--------------------------------------------------------------*/
  cdouble Sum=0.0;
  if (TMW->WhichG==TM_EIKR_OVER_R)
   { for(int Alpha=AlphaMin; Alpha<=AlphaMax; Alpha++)
      for(int i=1; i<=2; i++)
       Sum += FactorialTable[Alpha+2]
              *(S[i][Alpha][0]*JVector[i][Alpha+4] + S[i][Alpha][1]*KVector[i][Alpha+4])
              / ( MinusikAlphaTable[Alpha+3] * MiAlphaTable[i][Alpha+4]);

   }
  else if (TMW->WhichG==TM_GRADEIKR_OVER_R)
   { for(int Alpha=AlphaMin; Alpha<=AlphaMax; Alpha++)
      for(int i=1; i<=2; i++)
       Sum -= (Alpha+2)*FactorialTable[Alpha]
              *(S[i][Alpha][0]*JVector[i][Alpha+4] + S[i][Alpha][1]*KVector[i][Alpha+4])
              / ( MinusikAlphaTable[Alpha+1] * MiAlphaTable[i][Alpha+4]);

   };

  *((cdouble *)f) = x[1]*Sum;

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
cdouble HighKTaylorDuffy(TaylorDuffyArgStruct *Args)
{
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  int WhichCase   = Args->WhichCase;
  int WhichG      = Args->WhichG;
  int WhichH      = Args->WhichH;
  cdouble GParam  = Args->GParam;
  double *V1      = Args->V1;
  double *V2      = Args->V2;
  double *V3      = Args->V3;
  double *V2P     = Args->V2P;
  double *V3P     = Args->V3P;
  double *Q       = Args->Q;
  double *QP      = Args->QP;
  double AbsTol   = Args->AbsTol;
  double RelTol   = Args->RelTol;
   
  TMWorkspace MyTMW, *TMW=&MyTMW;
  TMW->WhichH = WhichH;
  TMW->WhichG = WhichG;
  TMW->GParam = GParam;

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
   /* the manual claims that certain parameters are not referenced    */
   /* in certain cases, which means the users might pass NULL for     */
   /* those parameters, which could cause core dumps unless we do the */
   /* following                                                       */
  if (WhichCase==TM_COMMONTRIANGLE)
   { V2P=V2; V3P=V3; } 
  else if (WhichCase==TM_COMMONEDGE)
   V2P=V2;

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

  if (WhichH==TM_DOT || WhichH==TM_DOTPLUS )
   {
     VecSub(V1,Q,D);
     VecSub(V1,QP,DP);
     TMW->AdD   = VecDot(A,D);
     TMW->AdDP  = VecDot(A,DP);
     TMW->BdDP  = VecDot(B,DP);
     TMW->APdD  = VecDot(AP,D);
     TMW->BPdD  = VecDot(BP,D);
     TMW->DdDP  = VecDot(D,DP);
   }
  else if (WhichH==TM_CROSS)
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
  /* 4. evaluate the function (for the common triangle case) or  */
  /*    the integral over x (common-edge case) or over x1,x2     */
  /*    (common-vertex case)                                     */
  /***************************************************************/
  static double Lower[2]={0.0, 0.0};
  static double Upper[2]={1.0, 1.0};
  cdouble Result, Error;
  TMW->nCalls=0;
  switch(WhichCase)
   { 
     case TM_COMMONTRIANGLE:
       Result = HKTD_CommonTriangle(TMW);
       break;

     case TM_COMMONEDGE:
       adapt_integrate(2, HKTD_CommonEdgeIntegrand, (void *)TMW, 1, Lower, Upper,
                       10000, AbsTol, RelTol, (double *)&Result, (double *)&Error);
       break;

     case TM_COMMONVERTEX:
       adapt_integrate(2, HKTD_CommonVertexIntegrand, (void *)TMW, 2, Lower, Upper,
                       10000, AbsTol, RelTol, (double *)&Result, (double *)&Error);
       break;
   };

/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
if (AICheck(__FILE__, __LINE__, (double *)&Result, (double *)&Error, AbsTol, RelTol, 2))
 { printf("Bawonkatage %i! \n",WhichCase);
   printf("V1 =(%e,%e,%e)\n",Args->V1[0],Args->V1[1],Args->V1[2]);
   printf("V2 =(%e,%e,%e)\n",Args->V2[0],Args->V2[1],Args->V2[2]);
   printf("V2P=(%e,%e,%e)\n",Args->V2P[0],Args->V2P[1],Args->V2P[2]);
   printf("V3 =(%e,%e,%e)\n",Args->V3[0],Args->V3[1],Args->V3[2]);
   printf("V3P=(%e,%e,%e)\n",Args->V3P[0],Args->V3P[1],Args->V3P[2]);
 };
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/

  return Result/(4.0*M_PI);

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void SipAlpha_One(const double *xVec, TMWorkspace *TMW, int WhichCase,
                  int *AlphaMin, int *AlphaMax, cdouble S[7][5][2])
{ 
  (void) xVec; (void) TMW; // unused

  switch(WhichCase)
   { 
     case TM_COMMONTRIANGLE: 
       *AlphaMin=0;
       *AlphaMax=2; 
       S[1][0][0]=S[2][0][0]=S[3][0][0]=1.0;
       S[1][0][1]=S[2][0][1]=S[3][0][1]=0.0;

       S[1][1][0]=S[2][1][0]=S[3][1][0]=-2.0;
       S[1][1][1]=S[2][1][1]=S[3][1][1]=0.0;

       S[1][2][0]=S[2][2][0]=S[3][2][0]=1.0;
       S[1][2][1]=S[2][2][1]=S[3][2][1]=0.0;
       return; 

     case TM_COMMONEDGE: 
       *AlphaMin=0;
       *AlphaMax=1; 
       S[1][0][0]=S[2][0][0]=S[3][0][0]=S[4][0][0]=S[5][0][0]=S[6][0][0]=1.0;
       S[1][0][1]=S[2][0][1]=S[3][0][1]=S[4][0][1]=S[5][0][1]=S[6][0][1]=0.0;

       S[1][1][0]=S[2][1][0]=S[3][1][0]=S[4][1][0]=S[5][1][0]=S[6][1][0]=-1.0;
       S[1][1][1]=S[2][1][1]=S[3][1][1]=S[4][1][1]=S[5][1][1]=S[6][1][1]=0.0;
       return;

     case TM_COMMONVERTEX: 
       *AlphaMin=0;
       *AlphaMax=0; 
       S[1][0][0]=S[2][0][0]=1.0;
       return;
   };

} 

/***************************************************************/
/***************************************************************/
/***************************************************************/
void SipAlpha_Dot(const double *xVec, TMWorkspace *TMW, int WhichCase,
                  int *AlphaMin, int *AlphaMax, cdouble S[7][5][2])
{
  double x1, x1_2, x2;
  double A2, B2, AdAP, AdB, AdBP, AdD, BdD, AdDP, APdD, BdAP, BdBP, BdDP, BPdD, DdDP;

  switch(WhichCase)
   { 
     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     case TM_COMMONTRIANGLE: 
       *AlphaMin=0;
       *AlphaMax=4; 

       A2=TMW->A2;
       B2=TMW->B2;
       AdB=TMW->AdB;
       AdD=TMW->APdD;
       AdDP=TMW->AdDP;
       BdD=TMW->BPdD;
       BdDP=TMW->BdDP;
       DdDP=TMW->DdDP;

       S[1][0][0] = (B2+3*A2+6*DdDP+2*BdDP+2*BdD+4*AdDP+4*AdD+3*AdB)/6.0;
       S[1][0][1] = 0.0;

       S[1][1][0] = (-4*B2-8*A2-12*DdDP-6*BdDP-6*BdD-9*AdDP -9*AdD-10*AdB)/6.0;
       S[1][1][1] = (2*B2+3*BdDP+3*BdD+4*AdB)/6.0;

       S[1][2][0] = B2+A2+DdDP+BdDP+BdD+AdDP+AdD+2*AdB;
       S[1][2][1] = (-B2-BdDP-BdD-2*AdB);

       S[1][3][0] = (-4*B2-2*BdDP-2*BdD-AdDP-AdD-6*AdB)/6.0;
       S[1][3][1] = (6*B2+3*BdDP+3*BdD+12*AdB)/6.0;

       S[1][4][0] = -(-B2+A2-AdB)/6.0;
       S[1][4][1] = -(2*B2+4*AdB)/6.0;
       
       S[2][0][0] = (B2+3*A2+6*DdDP+2*BdDP+2*BdD+4*AdDP+4*AdD+3*AdB)/6.0;
       S[2][0][1] = 0.0;
 
       S[2][1][0] = -(+2*B2+4*A2+12*DdDP+3*BdDP+3*BdD+6*AdDP+6*AdD+4*AdB)/6.0;
       S[2][1][1] = -((2*B2+4*A2+3*BdDP+3*BdD+3*AdDP+3*AdD+6*AdB))/6.0;

       S[2][2][0] = DdDP;
       S[2][2][1] = (B2+A2+BdDP+BdD+AdDP+AdD+2*AdB);

       S[2][3][0] = -(-2*B2-BdDP-BdD-2*AdDP-2*AdD)/6.0;
       S[2][3][1] = -((6*B2+3*BdDP+3*BdD+3*AdDP+3*AdD+6*AdB))/6.0;

       S[2][4][0] = (-B2+A2+AdB)/6;
       S[2][4][1] = (2*B2-2*A2)/6;
       
       S[3][0][0] = (B2+3*A2+6*DdDP+2*BdDP+2*BdD+4*AdDP+4*AdD+3*AdB)/6;
       S[3][0][1] = 0.0;

       S[3][1][0] = -( +2*B2+4*A2+12*DdDP+3*BdDP+3*BdD+6*AdDP+6*AdD+4*AdB )/6.0;
       S[3][1][1] = -(4*A2+3*AdDP+3*AdD+2*AdB)/6.0;

       S[3][2][0] = DdDP;
       S[3][2][1] = (A2+AdDP+AdD);

       S[3][3][0] = -(-2*B2-BdDP-BdD-2*AdDP-2*AdD)/6.0;
       S[3][3][1] = -((3*AdDP+3*AdD-6*AdB))/6.0;

       S[3][4][0] = -(B2-A2-AdB)/6.0;
       S[3][4][1] = -(2*A2+4*AdB)/6.0;

       return;

     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     case TM_COMMONEDGE: 
       *AlphaMin=0;
       *AlphaMax=3; 

       x1=xVec[0];
       x1_2=x1*x1;

       A2=TMW->A2;
       AdB=TMW->AdB;
       AdBP=TMW->AdBP;
       AdD=TMW->AdD;
       AdDP=TMW->AdDP;
       BdBP=TMW->BdBP;
       BdDP=TMW->BdDP;
       BPdD=TMW->BPdD;
       DdDP=TMW->DdDP;

       S[1][0][0] = (2*A2+6*DdDP+3*AdDP+3*AdD)/6;
       S[1][0][1] = 0.0;

       S[1][1][0] = ((-A2-2*BdDP-2*BPdD-2*AdD-AdBP-AdB)*x1
                      -2*DdDP+2*BdDP+2*BPdD+AdBP+AdB)/2;
       S[1][1][1] = (2*BdDP+AdB)*x1/2.0;

       S[1][2][0] = -((-2*BdBP-2*AdB)*x1_2+(-2*BdDP+4*BdBP-2*BPdD-2*AdD+2*AdB)*x1
                      +2*BdDP -2*BdBP+2*BPdD+AdDP+AdD)/2;
       S[1][2][1] = -((2*BdBP+2*AdB)*x1_2+(2*BdDP-2*BdBP)*x1)/2;

       S[1][3][0] = ((-6*BdBP-6*AdB)*x1_2+(3*A2+12*BdBP+3*AdBP+9*AdB)*x1
                     -2*A2-6*BdBP -3*AdBP-3*AdB)/6;
       S[1][3][1] = (((6*BdBP+6*AdB)*x1_2+(-6*BdBP-3*AdB)*x1))/6.0;
        
       S[2][0][0] = (2*A2+6*DdDP+3*AdDP+3*AdD)/6;
       S[2][0][1] = 0.0;

       S[2][1][0] = ((-A2-2*BdDP-2*BPdD-2*AdDP-AdBP-AdB)*x1
                      -2*DdDP+2*BdDP+2*BPdD+AdBP+AdB)/2;
       S[2][1][1] = (2.0*BPdD+AdBP)*x1/2.0;

       S[2][2][0] = -((-2*BdBP-2*AdBP)*x1_2+(-2*BdDP+4*BdBP-2*BPdD-2*AdDP
                       +2*AdBP)*x1+2*BdDP-2*BdBP+2*BPdD+AdDP+AdD)/2;
       S[2][2][1] = -((2*BdBP+2*AdBP)*x1_2+(-2*BdBP+2*BPdD)*x1)/2;

       S[2][3][0] = ((-6*BdBP-6*AdBP)*x1_2+(3*A2+12*BdBP+9*AdBP+3*AdB)*x1
                     -2*A2-6*BdBP-3*AdBP-3*AdB)/6;
       S[2][3][1] = (((6*BdBP+6*AdBP)*x1_2+(-6*BdBP-3*AdBP)*x1))/6.0;
        
        
       S[3][0][0] = (2*A2+6*DdDP+3*AdDP+3*AdD)/6;
       S[3][0][1] = 0.0;

       S[3][1][0] = -((2*BdDP+AdB)*x1+2*DdDP -2*BdDP-2*BPdD -AdBP-AdB)/2;
       S[3][1][1] = -((A2+2*BPdD+2*AdD+AdBP)*x1)/2;

       S[3][2][0] = ((2*BdDP-2*BdBP)*x1-2*BdDP+2*BdBP-2*BPdD-AdDP-AdD)/2;
       S[3][2][1] = ((2*BdBP+2*AdB)*x1_2+(-2*BdBP+2*BPdD+2*AdD-2*AdB)*x1)/2.0;

       S[3][3][0] = -((-6*BdBP-3*AdB)*x1+2*A2+6*BdBP+3*AdBP+3*AdB)/6;
       S[3][3][1] = -((6*BdBP+6*AdB)*x1_2+(-3*A2-6*BdBP-3*AdBP-6*AdB)*x1)/6.0;
        
       S[4][0][0] = (2*A2+6*DdDP+3*AdDP+3*AdD)/6;
       S[4][0][1] = 0.0;

       S[4][1][0] = -((2*BPdD+AdBP)*x1+2*DdDP-2*BdDP-2*BPdD-AdBP-AdB)/2;
       S[4][1][1] = -((A2+2*BdDP+2*AdDP+AdB)*x1)/2.0;

       S[4][2][0] = (+(-2*BdBP+2*BPdD)*x1-2*BdDP+2*BdBP-2*BPdD-AdDP-AdD)/2;
       S[4][2][1] = ((2*BdBP+2*AdBP)*x1_2+(2*BdDP-2*BdBP+2*AdDP-2*AdBP)*x1)/2.0;

       S[4][3][0] = -((-6*BdBP-3*AdBP)*x1+2*A2+6*BdBP+3*AdBP+3*AdB)/6;
       S[4][3][1] = -((6*BdBP+6*AdBP)*x1_2+(-3*A2-6*BdBP-6*AdBP-3*AdB)*x1)/6.0;
        
        
       S[5][0][0] = (2*A2+6*DdDP+3*AdDP+3*AdD)/6;
       S[5][0][1] = 0.0;

       S[5][1][0] = -((2*BPdD+AdBP)*x1+2*DdDP-2*BdDP-2*BPdD -AdBP-AdB)/2;
       S[5][1][1] = -(A2+2*AdD)*x1/2;

       S[5][2][0] = ((-2*BdBP+2*BPdD)*x1-2*BdDP+2*BdBP -2*BPdD-AdDP-AdD)/2; 
       S[5][2][1] = (2*AdD-2*AdB)*x1/2;

       S[5][3][0] = ((6*BdBP+3*AdBP)*x1-2*A2-6*BdBP-3*AdBP -3*AdB)/6;
       S[5][3][1] = (3*A2+6*AdB)*x1/6.0;
       
       
       S[6][0][0] = (2*A2+6*DdDP+3*AdDP+3*AdD)/6;
       S[6][0][1] = 0.0;

       S[6][1][0] = -((2*BdDP+AdB)*x1+2*DdDP-2*BdDP-2*BPdD-AdBP-AdB)/2;
       S[6][1][1] = -(A2+2*AdDP)*x1/2;

       S[6][2][0] = ((2*BdDP-2*BdBP)*x1-2*BdDP+2*BdBP -2*BPdD -AdDP-AdD)/2;
       S[6][2][1] = (2*AdDP-2*AdBP)*x1/2.0;

       S[6][3][0] = ((6*BdBP+3*AdB)*x1-2*A2-6*BdBP-3*AdBP -3*AdB)/6;
       S[6][3][1] = (3*A2+6*AdBP)*x1/6.0;
       
       return;
        
     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     case TM_COMMONVERTEX: 
       *AlphaMin=0;
       *AlphaMax=2;

       x1=xVec[0];
       x2=xVec[1];

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

       S[1][0][0] = DdDP;
       S[1][0][1] = 0.0;

       S[1][1][0] = APdD*x2+BdDP*x1+AdDP;
       S[1][1][1] = BPdD*x2;

       S[1][2][0] = (BdAP*x1+AdAP)*x2;
       S[1][2][1] = (BdBP*x1+AdBP)*x2;

       S[2][0][0] = DdDP;
       S[2][0][1] = 0.0;

       S[2][1][0] = AdDP*x2+BPdD*x1+APdD;
       S[2][1][1] = BdDP*x2;

       S[2][2][0] = (AdBP*x1+AdAP)*x2;
       S[2][2][1] = (BdBP*x1+BdAP)*x2;
              
       return;

   }; // switch(WhichCase)
}

/***************************************************************/
/* this is equivalent to 'Dot' - 4/(K*K) * 'One'               */
/* where K is the parameter that enters into the kernela       */
/* e^{i*K*r}/(4*pi*r)                                          */
/***************************************************************/
void SipAlpha_DotPlus(const double *xVec, TMWorkspace *TMW, int WhichCase,
                      int *AlphaMin, int *AlphaMax, cdouble S[7][5][2])
{

   /***************************************************************/
   /* get the 'Dot' contributions *********************************/
   /***************************************************************/
   SipAlpha_Dot(xVec, TMW, WhichCase, AlphaMin, AlphaMax, S);

   /***************************************************************/
   /* add in the 'One' contributions ******************************/
   /***************************************************************/
   cdouble Factor = -4.0/(TMW->GParam*TMW->GParam);
   switch(WhichCase)
     { 
       case TM_COMMONTRIANGLE: 
         S[1][0][0] += Factor;
         S[2][0][0] += Factor;
         S[3][0][0] += Factor;

         S[1][1][0] += -2.0*Factor;
         S[2][1][0] += -2.0*Factor;
         S[3][1][0] += -2.0*Factor;

         S[1][2][0] += Factor;
         S[2][2][0] += Factor;
         S[3][2][0] += Factor;
         break;

       case TM_COMMONEDGE: 
         S[1][0][0] += Factor;
         S[2][0][0] += Factor;
         S[3][0][0] += Factor;
         S[4][0][0] += Factor;
         S[5][0][0] += Factor;
         S[6][0][0] += Factor;

         S[1][1][0] -= Factor;
         S[2][1][0] -= Factor;
         S[3][1][0] -= Factor;
         S[4][1][0] -= Factor;
         S[5][1][0] -= Factor;
         S[6][1][0] -= Factor;
         break;

       case TM_COMMONVERTEX: 
         S[1][0][0] += Factor;
         S[2][0][0] += Factor;
         break;

     }; //switch(WhichCase)

}
               
/***************************************************************/
/***************************************************************/
/***************************************************************/
void SipAlpha_Cross(const double *xVec, TMWorkspace *TMW, int WhichCase,
                    int *AlphaMin, int *AlphaMax, cdouble S[7][5][2])
{
  double x1, x1_2, x2;

  double AdQxQP, APdQxQP, BdQxQP, BPdQxQP, LdQxQP; 
  double V1xAdQmQP, V1xBdQmQP, V1xAPdQmQP, V1xBPdQmQP; 
  double AxAPdQmQP, BxAPdQmQP, AxBdQmQP, AxBPdQmQP, BxBPdQmQP;

  switch(WhichCase)
   { 
     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     case TM_COMMONTRIANGLE: 
       *AlphaMin=0; 
       *AlphaMax=-1;
       return;

     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     case TM_COMMONEDGE: 
       *AlphaMin=1; 
       *AlphaMax=3;

       x1=xVec[0];
       x1_2=x1*x1;

       AdQxQP     = TMW->AdQxQP;
       BPdQxQP    = TMW->BPdQxQP;
       LdQxQP     = TMW->LdQxQP;
       V1xAdQmQP  = TMW->V1xAdQmQP;
       V1xBdQmQP  = TMW->V1xBdQmQP;
       V1xBPdQmQP = TMW->V1xBPdQmQP;
       AxBdQmQP   = TMW->AxBdQmQP;
       AxBPdQmQP  = TMW->AxBPdQmQP;
       BxBPdQmQP  = TMW->BxBPdQmQP;

       S[1][1][0] = -((-2*V1xBdQmQP+2*V1xBPdQmQP+2*V1xAdQmQP-2*LdQxQP-AxBdQmQP
                      +AxBPdQmQP-2*AdQxQP)*x1+2*V1xBdQmQP-2*V1xBPdQmQP
                      +2*LdQxQP+AxBdQmQP-AxBPdQmQP)/2;

       S[1][1][1] = -((2*V1xBdQmQP+2*LdQxQP-2*BPdQxQP+AxBdQmQP)*x1)/2.0;

       S[1][2][0] = +(BxBPdQmQP-AxBdQmQP)*x1_2
                    +(-V1xBdQmQP+V1xBPdQmQP+V1xAdQmQP-LdQxQP-2*BxBPdQmQP+AxBdQmQP-AdQxQP)*x1
                    +V1xBdQmQP-V1xBPdQmQP+LdQxQP+BxBPdQmQP;

       S[1][2][1] = (-BxBPdQmQP+AxBdQmQP)*x1_2+(V1xBdQmQP+LdQxQP+BxBPdQmQP-BPdQxQP)*x1;

       S[1][3][0] = ((-2*BxBPdQmQP+2*AxBdQmQP)*x1_2+
                     (4*BxBPdQmQP-3*AxBdQmQP+AxBPdQmQP)*x1
                     -2*BxBPdQmQP+AxBdQmQP-AxBPdQmQP)/2;

       S[1][3][1] = ((2*BxBPdQmQP-2*AxBdQmQP)*x1_2+(-2*BxBPdQmQP+AxBdQmQP)*x1)/2.0;
       
       S[2][1][0] = (+(2*V1xBdQmQP-2*V1xBPdQmQP+2*V1xAdQmQP+2*LdQxQP+AxBdQmQP-AxBPdQmQP
                      -2*AdQxQP)*x1
                     -2*V1xBdQmQP+2*V1xBPdQmQP-2*LdQxQP-AxBdQmQP+AxBPdQmQP)/2;

       S[2][1][1] = (2*V1xBPdQmQP-2*BPdQxQP+AxBPdQmQP)*x1/2.0;

       S[2][2][0] = +(BxBPdQmQP+AxBPdQmQP)*x1_2
                    +(-V1xBdQmQP+V1xBPdQmQP-V1xAdQmQP-LdQxQP-2*BxBPdQmQP-AxBPdQmQP+AdQxQP)*x1
                    +V1xBdQmQP-V1xBPdQmQP+LdQxQP+BxBPdQmQP;

       S[2][2][1] = ((-BxBPdQmQP-AxBPdQmQP)*x1_2+(-V1xBPdQmQP+BxBPdQmQP+BPdQxQP)*x1);


       S[2][3][0] = ((-2*BxBPdQmQP-2*AxBPdQmQP)*x1_2+(4*BxBPdQmQP-AxBdQmQP+3*AxBPdQmQP)*x1
                      -2*BxBPdQmQP+AxBdQmQP-AxBPdQmQP)/2;
       S[2][3][1] = ((2*BxBPdQmQP+2*AxBPdQmQP)*x1_2+(-2*BxBPdQmQP-AxBPdQmQP)*x1)/2.0;
       
       S[3][1][0] = -((-2*V1xBdQmQP-2*LdQxQP+2*BPdQxQP-AxBdQmQP)*x1+2*V1xBdQmQP
                       -2*V1xBPdQmQP+2*LdQxQP+AxBdQmQP-AxBPdQmQP)/2;

       S[3][1][1] = -(2*V1xBPdQmQP+2*V1xAdQmQP-2*BPdQxQP+AxBPdQmQP-2*AdQxQP)*x1/2.0;


       S[3][2][0] = +(-V1xBdQmQP-LdQxQP-BxBPdQmQP+BPdQxQP)*x1
                    +V1xBdQmQP-V1xBPdQmQP+LdQxQP +BxBPdQmQP;

       S[3][2][1] = (BxBPdQmQP-AxBdQmQP)*x1_2
                     +(V1xBPdQmQP+V1xAdQmQP-BxBPdQmQP-BPdQxQP +AxBdQmQP-AdQxQP)*x1;

       S[3][3][0] = -((-2*BxBPdQmQP+AxBdQmQP)*x1+2*BxBPdQmQP-AxBdQmQP+AxBPdQmQP)/2;

       S[3][3][1] = -(((2*BxBPdQmQP-2*AxBdQmQP)*x1_2
                    +(-2*BxBPdQmQP+2*AxBdQmQP-AxBPdQmQP)*x1))/2;

       
       S[4][1][0] = (+(-2*V1xBPdQmQP+2*BPdQxQP-AxBPdQmQP)*x1-2*V1xBdQmQP+2*V1xBPdQmQP
                   -2*LdQxQP-AxBdQmQP+AxBPdQmQP)/2;

       S[4][1][1] = ((2*V1xBdQmQP+2*V1xAdQmQP+2*LdQxQP-2*BPdQxQP+AxBdQmQP-2*AdQxQP)*x1)/2.0;


       S[4][2][0] = +(V1xBPdQmQP-BxBPdQmQP-BPdQxQP)*x1
                    +V1xBdQmQP-V1xBPdQmQP+LdQxQP+BxBPdQmQP;

       S[4][2][1] = (BxBPdQmQP+AxBPdQmQP)*x1_2
                   +(-V1xBdQmQP-V1xAdQmQP-LdQxQP-BxBPdQmQP+BPdQxQP-AxBPdQmQP+AdQxQP)*x1;

       S[4][3][0] = -((-2*BxBPdQmQP-AxBPdQmQP)*x1+2*BxBPdQmQP-AxBdQmQP+AxBPdQmQP)/2;

       S[4][3][1] = -( (  (2*BxBPdQmQP+2*AxBPdQmQP)*x1_2
                         +(-2*BxBPdQmQP+AxBdQmQP-2*AxBPdQmQP)*x1))/2.0;

       S[5][1][0] = -((2*V1xBPdQmQP-2*BPdQxQP+AxBPdQmQP)*x1
                      +2*V1xBdQmQP-2*V1xBPdQmQP+2*LdQxQP+AxBdQmQP-AxBPdQmQP)/2;
       S[5][1][1] = -(2*V1xAdQmQP-2*AdQxQP)*x1/2.0;

       S[5][2][0] = +(V1xBPdQmQP-BxBPdQmQP-BPdQxQP)*x1
                     +V1xBdQmQP-V1xBPdQmQP+LdQxQP+BxBPdQmQP;
       S[5][2][1] = (V1xAdQmQP+AxBdQmQP-AdQxQP)*x1;

       S[5][3][0] = -((-2*BxBPdQmQP-AxBPdQmQP)*x1+2*BxBPdQmQP-AxBdQmQP+AxBPdQmQP)/2;
       S[5][3][1] = -AxBdQmQP*x1;
       
       S[6][1][0] = ( (2*V1xBdQmQP+2*LdQxQP-2*BPdQxQP +AxBdQmQP)*x1
                      -2*V1xBdQmQP+2*V1xBPdQmQP-2*LdQxQP -AxBdQmQP+AxBPdQmQP)/2;
       S[6][1][1] = (2*V1xAdQmQP-2*AdQxQP)*x1/2.0;

       S[6][2][0] = +(-V1xBdQmQP-LdQxQP-BxBPdQmQP+BPdQxQP)*x1
                    +V1xBdQmQP-V1xBPdQmQP+LdQxQP+BxBPdQmQP;

       S[6][2][1] = (-V1xAdQmQP-AxBPdQmQP+AdQxQP)*x1;

       S[6][3][0] = ((2*BxBPdQmQP-AxBdQmQP)*x1-2*BxBPdQmQP+AxBdQmQP -AxBPdQmQP)/2;

       S[6][3][1] = AxBPdQmQP*x1;

       return;

     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     case TM_COMMONVERTEX: 

       *AlphaMin=1; 
       *AlphaMax=2;

       x1=xVec[0];
       x2=xVec[1];

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

       S[1][1][0] = +(V1xAPdQmQP-APdQxQP)*x2
                    +(-V1xBdQmQP+BdQxQP)*x1-V1xAdQmQP+AdQxQP;

       S[1][1][1] = (V1xBPdQmQP-BPdQxQP)*x2;

       S[1][2][0] = (BxAPdQmQP*x1+AxAPdQmQP)*x2;
       S[1][2][1] = (BxBPdQmQP*x1+AxBPdQmQP)*x2;

       S[2][1][0] = +(-V1xAdQmQP+AdQxQP)*x2 +(V1xBPdQmQP-BPdQxQP)*x1+V1xAPdQmQP-APdQxQP;
       S[2][1][1] = (-V1xBdQmQP+BdQxQP)*x2;

       S[2][2][0] = (AxBPdQmQP*x1+AxAPdQmQP)*x2;
       S[2][2][1] = (BxBPdQmQP*x1+BxAPdQmQP)*x2;

       return;

   };
}

} // namespace scuff
