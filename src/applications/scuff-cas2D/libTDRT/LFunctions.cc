/*
 * LFunctions.cc -- code for computing LFunctions between TDRT basis functions
 * 
 * homer reid    -- 11/2008 -- 10/2010
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "libTDRT.h"
#include "libhrutil.h"

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif
#ifdef USE_PTHREAD
#  include <pthread.h>
#endif

/***************************************************************/
/* Compute the contributions to the L-functions from a         */
/* single pair of line segments, and add these contributions   */
/* to L, dLdX, and dLdY.                                       */
/* Note that dLdX, dLdY are the derivatives of L with respect  */
/* to X and Y displacements of the FIRST basis function (and   */
/* thus the negatives of the derivatives wrt displacements of  */
/* the second basis function).                                 */
/***************************************************************/
void AddSegmentSegmentIntegrals(double *A, double *B, double Sigma,
                                double *C, double *D, double Tau,
                                double Kappa, double q, 
                                StaticSSIDataRecord *SSSIDR,
                                LFBuffer *L, LFBuffer *dLdX, LFBuffer *dLdY)
{ 
  double *Xs, *Xe, *Xsp, *Xep;
  StaticSSIDataRecord SSSIDRBuffer;

  double K2, AlphaKQ2, AlphaKQ;
  double uHat[2], upHat[2], R0[2];
  double LL[2], LLP[2];
  double l, lp, SinTheta, CosTheta;

  IPQRData MyIPQRData, *IPQRD=&MyIPQRData;

  double dI_00_0, dI_01_0, dI_10_0, dI_11_0;
  double dI_11_1, dI_21_1, dI_12_1;

  double T[6], dTdX[6], dTdY[6];

  int NeedDerivatives;

  /*****************************************************************/
  /* do some preliminary setup *************************************/
  /*****************************************************************/
  K2=Kappa*Kappa;
  AlphaKQ2=Kappa*Kappa+q*q;
  AlphaKQ=sqrt(AlphaKQ2);

  VecSub(B,A,uHat);
  l=VecNormalize(uHat);

  VecSub(D,C,upHat);
  lp=VecNormalize(upHat);

  SinTheta = (uHat[0]*upHat[1] - uHat[1]*upHat[0]);
  CosTheta = (uHat[0]*upHat[0] + uHat[1]*upHat[1]);

  NeedDerivatives = (dLdX==0 && dLdY==0) ? 0 : 1;

  /*****************************************************************/
  /* figure out whether we are traversing the line segments in the */
  /* forward or backward directions and set endpoints accordingly  */
  /*****************************************************************/
  if (Sigma==1.0)
   { Xs=A; Xe=B; }
  else
   { Xs=B; Xe=A; }

  if (Tau==1.0)
   { Xsp=C; Xep=D; }
  else
   { Xsp=D; Xep=C; }

  /***************************************************************/
  /* evaluate the actual 2D integrals over u,uprime.             */
  /***************************************************************/
  if ( AlphaKQ > 10.0 )
   { 
     if ( Xs==Xsp && Xe==Xep )
      uupIntegralSameSegment(Xs, Xe, AlphaKQ, 0, IPQRD);
     else if ( Xs==Xep && Xe==Xsp )
      uupIntegralSameSegment(Xs, Xe, AlphaKQ, 1, IPQRD);
     else if ( Xs==Xsp || Xs==Xep || Xe==Xsp || Xe==Xep ) 
      uupIntegralDuffy(Xs, Xe, Xsp, Xep, AlphaKQ, IPQRD);
     else
      uupIntegralCubature(Xs, Xe, Xsp, Xep, AlphaKQ, 7, NeedDerivatives, 0, IPQRD);
   }
  else
   { 
     /***************************************************************/
     /* if we need to desingularize the integrals, but the caller   */
     /* failed to provide us with a set of static segment-segment   */
     /* integrals for this pair of line segments, then we need to   */
     /* compute those now on-the-fly.                               */
     /***************************************************************/
     if ( SSSIDR==0 && Nearby(Xs, Xe, Xsp, Xep) )
      SSSIDR=ComputeStaticSSIData(Xs, Xe, Xsp, Xep, NeedDerivatives, &SSSIDRBuffer);
     uupIntegralCubature(Xs, Xe, Xsp, Xep, AlphaKQ, 7, NeedDerivatives, SSSIDR, IPQRD);
   };

  /***************************************************************/
  /* accumulate the contributions of the segment-segment         */
  /* integrals to the L-functions.                               */
  /* note that LPPTimes, LPZNabla, and LZPNabla contain factors  */
  /* of i which are not included here (they are inserted         */
  /* at the matrix-assembly stage.)                              */
  /***************************************************************/
  T[0] = IPQRD->I_00_0;
  T[1] = IPQRD->I_10_0;
  T[2] = IPQRD->I_01_0;
  T[3] = IPQRD->I_11_0;
  T[4] = IPQRD->I_X_1;
  T[5] = IPQRD->I_Y_1;

  L->LPPBullet += CosTheta*T[3];
  L->LPPNabla  += Sigma*Tau*T[0] / (l*lp);
  L->LPPTimes  += q*SinTheta*T[3];

  L->LPZBullet += 0.0;
  L->LPZNabla  += Sigma*q*T[2] / l;
  L->LPZTimes  += uHat[1]*T[4] - uHat[0]*T[5];

  L->LZPBullet += 0.0;
  L->LZPNabla  += Tau*q*T[1] / lp;
  L->LZPTimes  += upHat[1]*T[4] - upHat[0]*T[5];

  L->LZZBullet += T[3];
  L->LZZNabla  += q*q*T[3];
  L->LZZTimes  += 0.0;

  VecSub(Xs,Xsp,R0);
  VecSub(Xe,Xs,LL);
  VecSub(Xep,Xsp,LLP);

  if (dLdX)
   { 
     dI_00_0 = -R0[0]*IPQRD->I_00_1 - LL[0]*IPQRD->I_10_1 + LLP[0]*IPQRD->I_01_1;
     dI_10_0 = -R0[0]*IPQRD->I_10_1 - LL[0]*IPQRD->I_20_1 + LLP[0]*IPQRD->I_11_1;
     dI_01_0 = -R0[0]*IPQRD->I_01_1 - LL[0]*IPQRD->I_11_1 + LLP[0]*IPQRD->I_02_1;
     dI_11_0 = -R0[0]*IPQRD->I_11_1 - LL[0]*IPQRD->I_21_1 + LLP[0]*IPQRD->I_12_1;
     dI_11_1 = -R0[0]*IPQRD->I_11_2 - LL[0]*IPQRD->I_21_2 + LLP[0]*IPQRD->I_12_2;
     dI_21_1 = -R0[0]*IPQRD->I_21_2 - LL[0]*IPQRD->I_31_2 + LLP[0]*IPQRD->I_22_2;
     dI_12_1 = -R0[0]*IPQRD->I_12_2 - LL[0]*IPQRD->I_22_2 + LLP[0]*IPQRD->I_13_2;

     dTdX[0] = dI_00_0;
     dTdX[1] = dI_10_0;
     dTdX[2] = dI_01_0;
     dTdX[3] = dI_11_0;
     dTdX[4] = IPQRD->I_11_1 + R0[0]*dI_11_1 + LL[0]*dI_21_1 - LLP[0]*dI_12_1;
     dTdX[5] = R0[1]*dI_11_1 + LL[1]*dI_21_1 - LLP[1]*dI_12_1;

     dLdX->LPPBullet += CosTheta*dTdX[3];
     dLdX->LPPNabla  += Sigma*Tau*dTdX[0] / (l*lp);
     dLdX->LPPTimes  += q*SinTheta*dTdX[3];
   
     dLdX->LPZBullet += 0.0;
     dLdX->LPZNabla  += Sigma*q*dTdX[2] / l;
     dLdX->LPZTimes  += uHat[1]*dTdX[4] - uHat[0]*dTdX[5];
   
     dLdX->LZPBullet += 0.0;
     dLdX->LZPNabla  += Tau*q*dTdX[1] / lp;
     dLdX->LZPTimes  += upHat[1]*dTdX[4] - upHat[0]*dTdX[5];
   
     dLdX->LZZBullet += dTdX[3];
     dLdX->LZZNabla  += q*q*dTdX[3];
     dLdX->LZZTimes  += 0.0;

   }; // if (dLdX) ...

  if (dLdY)
   { 
     dI_00_0 = -R0[1]*IPQRD->I_00_1 - LL[1]*IPQRD->I_10_1 + LLP[1]*IPQRD->I_01_1;
     dI_10_0 = -R0[1]*IPQRD->I_10_1 - LL[1]*IPQRD->I_20_1 + LLP[1]*IPQRD->I_11_1;
     dI_01_0 = -R0[1]*IPQRD->I_01_1 - LL[1]*IPQRD->I_11_1 + LLP[1]*IPQRD->I_02_1;
     dI_11_0 = -R0[1]*IPQRD->I_11_1 - LL[1]*IPQRD->I_21_1 + LLP[1]*IPQRD->I_12_1;
     dI_11_1 = -R0[1]*IPQRD->I_11_2 - LL[1]*IPQRD->I_21_2 + LLP[1]*IPQRD->I_12_2;
     dI_21_1 = -R0[1]*IPQRD->I_21_2 - LL[1]*IPQRD->I_31_2 + LLP[1]*IPQRD->I_22_2;
     dI_12_1 = -R0[1]*IPQRD->I_12_2 - LL[1]*IPQRD->I_22_2 + LLP[1]*IPQRD->I_13_2;

     dTdY[0] = dI_00_0;
     dTdY[1] = dI_10_0;
     dTdY[2] = dI_01_0;
     dTdY[3] = dI_11_0;
     dTdY[4] = R0[0]*dI_11_1 + LL[0]*dI_21_1 - LLP[0]*dI_12_1;
     dTdY[5] = IPQRD->I_11_1 + R0[1]*dI_11_1 + LL[1]*dI_21_1 - LLP[1]*dI_12_1;

     dLdY->LPPBullet += CosTheta*dTdY[3];
     dLdY->LPPNabla  += Sigma*Tau*dTdY[0] / (l*lp);
     dLdY->LPPTimes  += q*SinTheta*dTdY[3];
   
     dLdY->LPZBullet += 0.0;
     dLdY->LPZNabla  += Sigma*q*dTdY[2] / l;
     dLdY->LPZTimes  += uHat[1]*dTdY[4] - uHat[0]*dTdY[5];
   
     dLdY->LZPBullet += 0.0;
     dLdY->LZPNabla  += Tau*q*dTdY[1] / lp;
     dLdY->LZPTimes  += upHat[1]*dTdY[4] - upHat[0]*dTdY[5];
   
     dLdY->LZZBullet += dTdY[3];
     dLdY->LZZNabla  += q*q*dTdY[3];
     dLdY->LZZTimes  += 0.0;

   }; // if (dLdY) ...

}
   
/***************************************************************/
/* Compute the L-functions between the TDRT basis functions    */
/* associated with a pair of control points (\rho_a, \rho_b).  */
/*                                                             */
/* INPUTS:                                                     */
/*                                                             */
/*  Oa, niva: control point \rho_a is interior vertex #niva    */
/*            within object Oa.                                */
/*                                                             */
/*  Ob, nivb: control point \rho_b is interior vertex #nivb    */
/*            within object Ob.                                */
/*                                                             */
/*  Kappa, q: wavevectors                                      */
/*                                                             */
/*    SSSIDT: may be NULL, but if not must be the return value */
/*            of a call to CreateStaticSSIDataTable() for the  */
/*            given objects.                                   */
/*                                                             */
/*  L, dLdX, dLdY: pointers to caller-allocated structures of  */
/*                 type LFBuffer.                              */
/*                                                             */
/*                 NOTE: dLdX, dLdY may be NULL, in which      */
/*                 case the corresponding derivatives are      */
/*                 not computed.                               */
/*                                                             */
/* OUTPUTS:                                                    */
/*                                                             */
/*  on return, L is populated with the values of the           */
/*  the L-functions, and dLdX and dLdY (if non-null) are       */
/*  populated with the X and/or Y derivatives of the           */
/*  L-functions.                                               */
/***************************************************************/
void ComputeLFunctions(TDRTObject *Oa, int niva, TDRTObject *Ob, int nivb,
                       double Kappa, double q, StaticSSIDataTable *SSSIDT,
                       LFBuffer *L, LFBuffer *dLdX, LFBuffer *dLdY)
{
  int iCVa, iEV1a, iEV2a;    /* indices of center and end vertices for BF A */
  double *CVa, *EV1a, *EV2a; /* center and end vertices for BF A */

  int iCVb, iEV1b, iEV2b;    /* indices of center and end vertices for BF B */ 
  double *CVb, *EV1b, *EV2b; /* center and end vertices for BF B */

  StaticSSIDataRecord SSSIDRBuffer, *SSSIDR;

  /***************************************************************/
  /* unpack vertices *********************************************/
  /***************************************************************/
  iCVa=Oa->IVs[niva]            ;  // index of center vertex of BF A 
  iEV1a=Oa->Neighbors[2*niva]   ;  // index of end-vertex 1 of BF A
  iEV2a=Oa->Neighbors[2*niva+1] ;  // index of end-vertex 2 of BF A 

  iCVb=Ob->IVs[nivb]            ;  // index of center vertex of BF B 
  iEV1b=Ob->Neighbors[2*nivb]   ;  // index of end-vertex 1 of BF B
  iEV2b=Ob->Neighbors[2*nivb+1] ;  // index of end-vertex 2 of BF B

  CVa=Oa->Vertices   + 2*iCVa;     // center vertex of BF A
  EV1a=Oa->Vertices  + 2*iEV1a;    // end-vertex 1 of BF A 
  EV2a=Oa->Vertices  + 2*iEV2a;    // end-vertex 2 of BF A

  CVb=Ob->Vertices   + 2*iCVb;     // center vertex of BF B
  EV1b=Ob->Vertices  + 2*iEV1b;    // end-vertex 1 of BF B 
  EV2b=Ob->Vertices  + 2*iEV2b;    // end-vertex 2 of BF B

  /***************************************************************/
  /* add contributions to matrix elements from the four pairs of */
  /* line segments                                               */
  /***************************************************************/
  memset(L,0,sizeof(LFBuffer));
  if (dLdX) memset(dLdX,0,sizeof(LFBuffer));
  if (dLdY) memset(dLdY,0,sizeof(LFBuffer));

  SSSIDR=GetStaticSSIData(SSSIDT, Oa, iEV1a, iCVa, Ob, iEV1b, iCVb, &SSSIDRBuffer); 
  AddSegmentSegmentIntegrals(EV1a, CVa, +1.0, EV1b, CVb, +1.0, Kappa, q, SSSIDR, L, dLdX, dLdY);

  SSSIDR=GetStaticSSIData(SSSIDT, Oa, iEV1a, iCVa, Ob, iEV2b, iCVb, &SSSIDRBuffer); 
  AddSegmentSegmentIntegrals(EV1a, CVa, +1.0, CVb, EV2b, -1.0, Kappa, q, SSSIDR, L, dLdX, dLdY);

  SSSIDR=GetStaticSSIData(SSSIDT, Oa, iEV2a, iCVa, Ob, iEV1b, iCVb, &SSSIDRBuffer); 
  AddSegmentSegmentIntegrals(CVa, EV2a, -1.0, EV1b, CVb, +1.0, Kappa, q, SSSIDR, L, dLdX, dLdY);

  SSSIDR=GetStaticSSIData(SSSIDT, Oa, iEV2a, iCVa, Ob, iEV2b, iCVb, &SSSIDRBuffer);
  AddSegmentSegmentIntegrals(CVa, EV2a, -1.0, CVb, EV2b, -1.0, Kappa, q, SSSIDR, L, dLdX, dLdY);

}
