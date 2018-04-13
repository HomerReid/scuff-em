/*
 * StaticSSIDataRecord.cc -- libTDRT routines for computing static 
 *                        -- segment-segment integrals 
 *
 *                        -- this file contains routines that work
 *                        -- with StaticSSIDataRecords only, and knows
 *                        -- nothing about StaticSSIDataTables; all of 
 *                        -- that stuff is now relegated to a separate file.
 *
 * homer reid -- 6/2009 -- 10/2010
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

#include <stdarg.h>

#include "libhrutil.h"
#include "libSGJC.h"
#include "libTDRT.h"

#include <math.h>

#define INTERVALS 10000
#define ABSTOL 1.0e-10
#define RELTOL 1.0e-8
#define NFUN 30

/***************************************************************/
/* workspace structure passed among these routines             */
/***************************************************************/
typedef struct StaticSSIWorkspace
 { double x, r0, CT, CP, CPP;
   double R0[2], L[2], LP[2];
   int NeedDerivatives;
 } StaticSSIWorkspace;

/***************************************************************/
/* 1D integrand routine that evaluates the portion of the      */
/* static line-line integrals that cannot be done analytically */
/*                                                             */
/* the original integrals are                                  */
/*                                                             */
/* J_{pq}^r = \int_0^1 \int_0^1 du dup u^p up^q S_r(R(u,up))   */
/*                                                             */
/* with S_{1,2,3,4}(R) = 1,log(R),1/R^2,1/R^4                  */
/*                                                             */
/* in each case we can perform one of the two integrations     */
/* analytically, leaving a single integral to be done by       */
/* numerical quadrature; this routine calculates the value of  */
/* the integrand of that integral.                             */
/*                                                             */
/* the packing of the integrals into the output vector obeys   */
/* the following correspondence:                               */
/*                                                             */
/* f[00] = J_{00}^2                                            */
/* f[01] = J_{10}^2                                            */
/* f[02] = J_{01}^2                                            */
/* f[03] = J_{11}^2                                            */
/* f[04] = J_{20}^2                                            */
/* f[05] = J_{21}^2                                            */
/* f[06] = J_{12}^2                                            */
/* f[07] = J_{02}^2                                            */
/* f[08] = J_{31}^2                                            */
/* f[09] = J_{22}^2                                            */
/* f[10] = J_{13}^2                                            */
/*                                                             */
/* f[11] = J_{00}^3                                            */
/* f[12] = J_{10}^3                                            */
/* f[13] = J_{01}^3                                            */
/* f[14] = J_{11}^3                                            */
/* f[15] = J_{20}^3                                            */
/* f[16] = J_{21}^3                                            */
/* f[17] = J_{12}^3                                            */
/* f[18] = J_{02}^3                                            */
/* f[19] = J_{31}^3                                            */
/* f[20] = J_{22}^3                                            */
/* f[21] = J_{13}^3                                            */
/*                                                             */
/* f[22] = J_{11}^4                                            */
/* f[23] = J_{21}^4                                            */
/* f[24] = J_{22}^4                                            */
/* f[25] = J_{12}^4                                            */
/* f[26] = J_{31}^4                                            */
/* f[27] = J_{13}^4                                            */
/*                                                             */
/* f[28] = J_{x}^3                                             */
/* f[29] = J_{y}^3                                             */
/*                                                             */
/* if NeedDerivatives=0, then f[6,7,8,12,13,14] are set to 0   */
/* to speed convergence, since those are the most singular     */
/* integrals and we only need them for derivative calculations.*/
/*                                                             */
/* the parameters u0, b that enter below are defined by the    */
/* expression                                                  */
/*                                                             */
/*  R^2(u,up) = l^2[ (u+u0^2) + b^2 ]                          */
/*                                                             */
/* where                                                       */
/*                                                             */
/*   u0=r0*CP - up*x*CT                                        */
/*   b2=r0^2 + x^2*up*up - 2*x*r0*CP*up - u0*u0                */
/*   r0=|R0|/l, x=lp/l, etc.                                   */
/***************************************************************/
static void StaticSSIntegrand(unsigned ndim, const double *upValues, 
                              void *params, unsigned fdim, double *f)
{ 
  (void)ndim; //unused
  (void)fdim; //unused
  double up;
  double x, r0, CP, CT, CPP;
  double x2, r02, u0, u02, u03, u04, b, b2, b3, b4;
  double u0P1, u0P12, u0P13, logb2Pu02, logb2Pu0P12, atanDelta, logDelta;
  double atanDeltaOverB;
  StaticSSIWorkspace *SSSIW=(StaticSSIWorkspace *)params;
  double A0, A1, A2, A3;
  double B0, B1, B2, B3;
  double C1, C2, C3;

  x=SSSIW->x;
  r0=SSSIW->r0;
  CP=SSSIW->CP;
  CT=SSSIW->CT;
  CPP=SSSIW->CPP;

  x2=x*x;
  r02=r0*r0;
  
  up=upValues[0];

  u0=r0*CP - up*x*CT; 
  u02=u0*u0;
  u03=u02*u0;
  u04=u02*u02;
  u0P1=u0+1.0;
  u0P12=u0P1*u0P1;
  u0P13=u0P12*u0P1;

  b2=r02 + x2*up*up - 2*x*r0*CPP*up - u0*u0;
  b=b2 < 1.0e-12 ? 0.0 : sqrt(b2);
  b3=b2*b;
  b4=b2*b2;
   
  //logb2Pu02=log(b2+u02);
  logb2Pu02 = (b2+u02)<1.0e-8 ? 0.0 : log(b2+u02);
  // logb2Pu0P12 = log(b2+u0P12);
  logb2Pu0P12 = (b2+u0P12)<1.0e-8 ? 0.0 : log(b2+u0P12);

  atanDelta=(b==0.0) ? 0.0 : atan(u0P1/b) - atan(u0/b);
  atanDeltaOverB=(b==0.0 ? 0.0 : atanDelta/b);
  logDelta=logb2Pu0P12-logb2Pu02;

  /***************************************************************/
  /* An = \int_0^1 du u^n Log [ Sqrt[ (u+u0)^2 + b^2 ] ]         */
  /* Bn = \int_0^1 du u^n / [ (u+u0)^2 + b^2 ]                   */
  /* Cn = \int_0^1 du u^n / [ (u+u0)^2 + b^2 ]^2                 */
  /***************************************************************/
  A0 = -1.0 + b*atanDelta + 0.5*u0P1*logb2Pu0P12 - 0.5*u0*logb2Pu02;
  A1 = -0.25 + 0.5*u0 - b*u0*atanDelta + 0.25*logb2Pu0P12 +0.25*(b2-u02)*logDelta;
  A2 = -1.0/9.0 + (b2-u02)/3.0 + u0/6.0 - b*(b2-3.0*u02)*atanDelta/3.0
       -(b2*u0/2.0 - u03/6.0)*logDelta +(logb2Pu0P12/6.0);
  A3 = -1/16.0 + b2*(1.0-6.0*u0)/8.0 + u0*(2.0-3.0*u0+6.0*u02)/24.0
       +u0*b*(b2-u02)*atanDelta -(b4-6.0*b2*u02+u04)*logDelta/8.0
       +logb2Pu0P12/8.0;

  if ( b < 1.0e-4 )
   { 
     B0 = 1.0/(u0*u0P1);
     B1 = -1.0/u0P1 + 0.5*logDelta;
     B2 = (1+2.0*u0)/u0P1 -u0*logDelta;
     B3 = 1.5 - 3.0*u0 - 1.0/u0P1 + 1.5*u02*logDelta;

     C1 = (1.0 + 3.0*u0) / (6.0*u02*u0P13);
     C2 = 1.0/ (3.0*u0*u0P13);
     C3 = -(11.0 + 3.0*u0*(5.0+2.0*u0))/(6.0*u0P13) + 0.5*logDelta;
   }
  else
   { 
     B0 = atanDeltaOverB;
     B1 = -u0*B0 + 0.5*logDelta;
     B2 = 1.0 + (u02-b2)*atanDeltaOverB - u0*logDelta;
     B3 = 0.5 - 2.0*u0 - u0*(u02-3.0*b2)*atanDeltaOverB 
          -0.5*(b2-3.0*u02)*logDelta;

     C1 = (b*u0P1 - u0*(b2+u0P12)*atanDelta) / (2.0*b3*(b2+u0P12));
     C2 = ((b2+u02)*(b2+u0P12)*atanDelta - b*(b2+u02+u0))/(2.0*b3*(b2+u0P12));
     C3 = (-1.0 + u0 + u0P1/((b2 + u0P12))) / (2.0*b2)
         - u0*(3.0*b2+u0*u0)*atanDelta/(2.0*b3)
         +0.5*logDelta;
   };

  // i'm pretty sure the following procedure yields correct answers
  // in all cases
  if (!std::isnormal(B0)) B0=0.0;
  if (!std::isnormal(B1)) B1=0.0;
  if (!std::isnormal(B2)) B2=0.0;
  if (!std::isnormal(B3)) B3=0.0;
  if (!std::isnormal(C1)) C1=0.0;
  if (!std::isnormal(C2)) C2=0.0;
  if (!std::isnormal(C3)) C3=0.0;

  /* log\rho terms */
  f[ 0  ] = A0; 
  f[ 1  ] = A1;
  f[ 2  ] = up*A0;
  f[ 3  ] = up*A1;
  f[ 4  ] = A2;
  f[ 5  ] = up*A2;
  f[ 6  ] = up*up*A1; 
  f[ 7  ] = up*up*A0;
  f[ 8  ] = up*A3;
  f[ 9  ] = up*up*A2;
  f[ 10 ] = up*up*up*A1;

  /* 1/R^2 terms */
  f[ 11 ] = B0; 
  f[ 12 ] = B1;
  f[ 13 ] = up*B0;
  f[ 14 ] = up*B1;
  f[ 15 ] = B2;
  f[ 16 ] = up*B2;
  f[ 17 ] = up*up*B1;
  f[ 18 ] = up*up*B0;
  f[ 19 ] = up*B3;
  f[ 20 ] = up*up*B2;
  f[ 21 ] = up*up*up*B1;

  /* 1/R^4 terms */
  f[ 22 ] = up*C1;
  f[ 23 ] = up*C2;
  f[ 24 ] = up*up*C2;
  f[ 25 ] = up*up*C1;
  f[ 26 ] = up*C3;
  f[ 27 ] = up*up*up*C1;

/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
 //    SSSIDR->J_X_3  = R0[0]*SSSIDR->J_11_3 + L[0]*SSSIDR->J_21_3 - LP[0]*SSSIDR->J_12_3;
 //    SSSIDR->J_Y_3  = R0[1]*SSSIDR->J_11_3 + L[1]*SSSIDR->J_21_3 - LP[1]*SSSIDR->J_12_3;
f[ 28 ] = SSSIW->R0[0]*f[ 14 ] + SSSIW->L[0]*f[ 16 ] - SSSIW->LP[0]*f[ 17 ];
f[ 29 ] = SSSIW->R0[1]*f[ 14 ] + SSSIW->L[1]*f[ 16 ] - SSSIW->LP[1]*f[ 17 ];
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/

  if (SSSIW->NeedDerivatives==0) 
   { f[11] = f[12] = f[13] = 0.0;
     f[22] = f[23] = f[24] = 0.0;
     f[25] = f[26] = f[27] = 0.0;
   };

}

/***************************************************************/
/* Compute static SSIs in the special case in which both u and */
/* up integrals are over the same line segment. (The full      */
/* integrals can be evaluated analytically in this case).      */
/*                                                             */
/* Flip=1 means the line segments are going in opposite        */
/* directions, which we handle by taking up \to (1-up) in the  */  
/* up integral.                                                */
/***************************************************************/
void ComputeStaticSSIData_SameSegment(double *Xs, double *Xe, int Flip,
                                      StaticSSIDataRecord *SSSIDR)
{
  double L=VecDistance(Xs,Xe), LogL=log(L);

  SSSIDR->J_00_1=1.0;
  SSSIDR->J_10_1=0.5;
  SSSIDR->J_01_1=0.5;
  SSSIDR->J_11_1=0.25;

  SSSIDR->J_00_2 = LogL - 1.5;
  SSSIDR->J_10_2 = 0.5*LogL-0.75;
  SSSIDR->J_01_2 = 0.5*LogL-0.75;
  SSSIDR->J_11_2 = 0.25*LogL - 7.0/16.0;

  SSSIDR->J_X_3=SSSIDR->J_X_2=SSSIDR->J_X_1=0.0;
  SSSIDR->J_Y_3=SSSIDR->J_Y_2=SSSIDR->J_Y_1=0.0;

  if (Flip)
   SSSIDR->J_11_2=SSSIDR->J_01_2 - SSSIDR->J_11_2;

  /* and we don't need to evaluate any more of the J functions  */ 
  /* because we never need derivatives in the same-segment case */ 
  
}

/*******************************************************************/
/* compute static SSIs for the line segments (Xs,Xe) and (Xsp,Xep) */
/* and return a pointer to the buffer into which we write them.    */
/*******************************************************************/
StaticSSIDataRecord *ComputeStaticSSIData(double *Xs, double *Xe, double *Xsp, double *Xep,
                                          int NeedDerivatives,
                                          StaticSSIDataRecord *SSSIDR)
{

  if ( Xs==Xsp && Xe==Xep )
   ComputeStaticSSIData_SameSegment(Xs, Xe, 0, SSSIDR);
  else if ( Xs==Xep && Xe==Xsp )
   ComputeStaticSSIData_SameSegment(Xs, Xe, 1, SSSIDR);
  else if ( Xs==Xsp || Xs==Xep || Xe==Xsp || Xe==Xep )
   ComputeStaticSSIData_CommonVertex(Xs, Xe, Xsp, Xep, SSSIDR);
  else
   { 
     double Val[NFUN], Err[NFUN];
     double L[2], LP[2], R0[2];
     double l, L2, L4, LogL, lp, x, x2, r0, r02, CT, CP, CPP;

     /* compute some geometrical parameters. */
     VecSub(Xe,Xs,L);
     VecSub(Xep,Xsp,LP);
     VecSub(Xs,Xsp,R0);
     l=VecNorm(L);
     L2=l*l;
     L4=L2*L2;
     LogL=log(l);
     lp=VecNorm(LP);
     x=lp/l;
     x2=x*x;
     r0=VecNorm(R0)/l;
     r02=r0*r0;
     CT  = (L[0]*LP[0] + L[1]*LP[1]) / (x*L2);
     CP  = r0<1.0e-8 ? 0.0 : (L[0]*R0[0]  + L[1]*R0[1]) / (L2*r0);
     CPP = r0<1.0e-8 ? 0.0 : (LP[0]*R0[0] + LP[1]*R0[1]) / (x*L2*r0);

     StaticSSIWorkspace MySSSIW, *SSSIW=&MySSSIW;
     SSSIW->x=x;
     SSSIW->r0=r0;
     SSSIW->CT=CT;
     SSSIW->CP=CP;
     SSSIW->CPP=CPP;
     SSSIW->NeedDerivatives=NeedDerivatives;
 
     /* these ones are easy */
     SSSIDR->J_00_1=1.0;
     SSSIDR->J_10_1=1.0/2.0;
     SSSIDR->J_01_1=1.0/2.0;
     SSSIDR->J_11_1=1.0/4.0;
     SSSIDR->J_20_1=1.0/3.0;
     SSSIDR->J_21_1=1.0/6.0;
     SSSIDR->J_12_1=1.0/6.0;
     SSSIDR->J_02_1=1.0/3.0;
     SSSIDR->J_31_1=1.0/8.0;
     SSSIDR->J_22_1=1.0/9.0;
     SSSIDR->J_13_1=1.0/8.0;
     SSSIDR->J_X_1=R0[0]*0.25 + (L[0]-LP[0])/6.0;
     SSSIDR->J_Y_1=R0[1]*0.25 + (L[1]-LP[1])/6.0;

/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
memcpy(SSSIW->L,L,2*sizeof(double));
memcpy(SSSIW->LP,LP,2*sizeof(double));
memcpy(SSSIW->R0,R0,2*sizeof(double));
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/

     /* the remaining integrals must be evaluated by quadrature */
     double Lower=0.0, Upper=1.0;
     adapt_integrate(NFUN, StaticSSIntegrand, (void *) SSSIW, 1,
                     &Lower, &Upper, 100000, ABSTOL, RELTOL, Val, Err);

     /* */
     SSSIDR->J_00_2 = Val[0]  + LogL;
     SSSIDR->J_10_2 = Val[1]  + LogL/2.0;
     SSSIDR->J_01_2 = Val[2]  + LogL/2.0;
     SSSIDR->J_11_2 = Val[3]  + LogL/4.0;
     SSSIDR->J_20_2 = Val[4]  + LogL/3.0;
     SSSIDR->J_21_2 = Val[5]  + LogL/6.0;
     SSSIDR->J_12_2 = Val[6]  + LogL/6.0;
     SSSIDR->J_02_2 = Val[7]  + LogL/3.0;
     SSSIDR->J_31_2 = Val[8]  + LogL/8.0;
     SSSIDR->J_22_2 = Val[9]  + LogL/9.0;
     SSSIDR->J_13_2 = Val[10] + LogL/8.0;
     SSSIDR->J_X_2  = R0[0]*SSSIDR->J_11_2 + L[0]*SSSIDR->J_21_2 - LP[0]*SSSIDR->J_12_2;
     SSSIDR->J_Y_2  = R0[1]*SSSIDR->J_11_2 + L[1]*SSSIDR->J_21_2 - LP[1]*SSSIDR->J_12_2;

     SSSIDR->J_00_3 = Val[11]/L2;
     SSSIDR->J_10_3 = Val[12]/L2;
     SSSIDR->J_01_3 = Val[13]/L2;
     SSSIDR->J_11_3 = Val[14]/L2; 
     SSSIDR->J_20_3 = Val[15]/L2;
     SSSIDR->J_21_3 = Val[16]/L2;
     SSSIDR->J_12_3 = Val[17]/L2; 
     SSSIDR->J_02_3 = Val[18]/L2;
     SSSIDR->J_31_3 = Val[19]/L2;
     SSSIDR->J_22_3 = Val[20]/L2;
     SSSIDR->J_13_3 = Val[21]/L2;
 //    SSSIDR->J_X_3  = R0[0]*SSSIDR->J_11_3 + L[0]*SSSIDR->J_21_3 - LP[0]*SSSIDR->J_12_3;
 //    SSSIDR->J_Y_3  = R0[1]*SSSIDR->J_11_3 + L[1]*SSSIDR->J_21_3 - LP[1]*SSSIDR->J_12_3;

SSSIDR->J_X_3  = Val[28]/L2;
SSSIDR->J_Y_3  = Val[29]/L2;

     SSSIDR->J_11_4 = Val[22]/L4;
     SSSIDR->J_21_4 = Val[23]/L4;
     SSSIDR->J_22_4 = Val[24]/L4;
     SSSIDR->J_12_4 = Val[25]/L4; 
     SSSIDR->J_31_4 = Val[26]/L4;
     SSSIDR->J_13_4 = Val[27]/L4;

   };

  return SSSIDR;

}  
