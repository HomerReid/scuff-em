/*
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

#include <gsl/gsl_sf.h>

#include <libhrutil.h>
#include "libTDRT.h"

#define NFUN 15

/*******************************************************************/
/* uhh ... arcane much?? *******************************************/
/*******************************************************************/
#define _00_2 0
#define _10_2 1
#define _01_2 2
#define _20_2 3
#define _11_2 4
#define _02_2 5
#define _21_2 6
#define _12_2 7
#define _10_3 8
#define _01_3 9
#define _20_3 10
#define _11_3 11
#define _02_3 12
#define _21_3 13
#define _12_3 14

/*******************************************************************/
/* An = \int_0^1 du u^n Log [ Sqrt[ (u+u0)^2 + b^2 ] ]             */
/* Bn = \int_0^1 du u^n / [ (u+u0)^2 + b^2 ]^2                     */
/* Cn = \int_0^1 du u^n / [ (u+u0)^2 + b^2 ]^4                     */
/*******************************************************************/
void GetABC(double u0, double b, double A[4], double B[4], double C[4])
{
  double u02, u03, u04, b2, b3, b4;
  double u0P1, u0P12, u0P13, logb2Pu02, logb2Pu0P12, atanDelta, logDelta;
  double atanDeltaOverB;

  u02=u0*u0;
  u03=u02*u0;
  u04=u02*u02;
  u0P1=u0+1.0;
  u0P12=u0P1*u0P1;
  u0P13=u0P12*u0P1;

  b2=b*b;
  b3=b2*b;
  b4=b2*b2;
   
  logb2Pu02 = (b2+u02)<1.0e-8 ? 0.0 : log(b2+u02);
  logb2Pu0P12 = (b2+u0P12)<1.0e-8 ? 0.0 : log(b2+u0P12);

  atanDelta=(b==0.0) ? 0.0 : atan(u0P1/b) - atan(u0/b);
  atanDeltaOverB=(b==0.0 ? 0.0 : atanDelta/b);
  logDelta=logb2Pu0P12-logb2Pu02;

  A[0] = -1.0 + b*atanDelta + 0.5*u0P1*logb2Pu0P12 - 0.5*u0*logb2Pu02;
  A[1] = -0.25 + 0.5*u0 - b*u0*atanDelta + 0.25*logb2Pu0P12 +0.25*(b2-u02)*logDelta;
  A[2] = -1.0/9.0 + (b2-u02)/3.0 + u0/6.0 - b*(b2-3.0*u02)*atanDelta/3.0
         -(b2*u0/2.0 - u03/6.0)*logDelta +(logb2Pu0P12/6.0);
  A[3] = -1/16.0 + b2*(1.0-6.0*u0)/8.0 + u0*(2.0-3.0*u0+6.0*u02)/24.0
         +u0*b*(b2-u02)*atanDelta -(b4-6.0*b2*u02+u04)*logDelta/8.0
         +logb2Pu0P12/8.0;

  if ( b < 1.0e-8 )
   { 
     B[0] = 1.0/(u0*u0P1);
     B[1] = -1.0/u0P1 + 0.5*logDelta;
     B[2] = (1+2.0*u0)/u0P1 -u0*logDelta;
     B[3] = 1.5 - 3.0*u0 - 1.0/u0P1 + 1.5*u02*logDelta;

     C[1] = (1.0 + 3.0*u0) / (6.0*u02*u0P13);
     C[2] = 1.0/ (3.0*u0*u0P13);
     C[3] = -(11.0 + 3.0*u0*(5.0+2.0*u0))/(6.0*u0P13) + 0.5*logDelta;
   }
  else
   { 
     B[0] = atanDeltaOverB;
     B[1] = -u0*B[0] + 0.5*logDelta;
     B[2] = 1.0 + (u02-b2)*atanDeltaOverB - u0*logDelta;
     B[3] = 0.5 - 2.0*u0 - u0*(u02-3.0*b2)*atanDeltaOverB 
            -0.5*(b2-3.0*u02)*logDelta;

     C[1] = (b*u0P1 - u0*(b2+u0P12)*atanDelta) / (2.0*b3*(b2+u0P12));
     C[2] = ((b2+u02)*(b2+u0P12)*atanDelta - b*(b2+u02+u0))/(2.0*b3*(b2+u0P12));
     C[3] = (-1.0 + u0 + u0P1/((b2 + u0P12))) / (2.0*b2)
            - u0*(3.0*b2+u0*u0)*atanDelta/(2.0*b3)
            +0.5*logDelta;
   };

}


/*******************************************************************/
/* this routine uses duffy-transform techniques to evaluate the    */
/* following integrals exactly                                     */
/*                                                                 */
/* I[ 0] = \int_0^1 \int_0^1 du dup log(R)                         */
/* I[ 1] = \int_0^1 \int_0^1 du dup u*log(R)                       */
/* I[ 2] = \int_0^1 \int_0^1 du dup up*log(R)                      */
/* I[ 3] = \int_0^1 \int_0^1 du dup u^2*log(R)                     */
/* I[ 4] = \int_0^1 \int_0^1 du dup u*up*log(R)                    */
/* I[ 5] = \int_0^1 \int_0^1 du dup up^2*log(R)                    */
/* I[ 6] = \int_0^1 \int_0^1 du dup u^2*up*log(R)                  */
/* I[ 7] = \int_0^1 \int_0^1 du dup u*up^2*log(R)                  */
/*                                                                 */
/* I[ 8] = \int_0^1 \int_0^1 du dup u/R^2                          */
/* I[ 9] = \int_0^1 \int_0^1 du dup up/R^2                         */
/* I[10] = \int_0^1 \int_0^1 du dup u^2/R^2                        */
/* I[11] = \int_0^1 \int_0^1 du dup u*up/R^2                       */
/* I[12] = \int_0^1 \int_0^1 du dup up^2/R^2                       */
/* I[13] = \int_0^1 \int_0^1 du dup u^2*up/R^2                     */
/* I[14] = \int_0^1 \int_0^1 du dup u*up^2/R^2                     */
/*                                                                 */
/* where R^2 = l^2*u^2 + lp^2*up^2 - 2*u*up*l*lp*CT                */  
/*                                                                 */
/*******************************************************************/
void GetDuffyIntegrals(double l, double lp, double CT, double I[NFUN])
{
  double x, y, ST;
  double A1[4], B1[4], A2[4], B2[4], C[4];
  double LogL, LogLP, L2, LP2;
  int p, q;
  double dp, dq;

  LogL=log(l);
  LogLP=log(lp);
  L2=l*l;
  LP2=lp*lp;

  x=lp/l;
  y=l/lp;
  ST=1.0-CT*CT;
  if (ST<1.0e-8)  
   ST=0.0;
  else
   ST=sqrt(ST);

  GetABC( -y*CT, y*ST, A1, B1, C);
  GetABC( -x*CT, x*ST, A2, B2, C);

  p=0; q=0; dp=(double)p; dq=(double)q;
  I[_00_2] = ( ( (dp+1.0)*LogLP + (dq+1.0)*LogL - 1.0 ) / ((dp+1.0)*(dq+1.0))
              + A1[q] + A2[p] ) / (dp+dq+2.0);

  p=1; q=0; dp=(double)p; dq=(double)q;
  I[_10_2] = ( ( (dp+1.0)*LogLP + (dq+1.0)*LogL - 1.0 ) / ((dp+1.0)*(dq+1.0))
              + A1[q] + A2[p] ) / (dp+dq+2.0);
  I[_10_3] = ( B1[q]/LP2 + B2[p]/L2 ) / (dp+dq); 

  p=0; q=1; dp=(double)p; dq=(double)q;
  I[_01_2] = ( ( (dp+1.0)*LogLP +(dq+1.0)*LogL - 1.0 ) / ((dp+1.0)*(dq+1.0))
              + A1[q] + A2[p] ) / (dp+dq+2.0);
  I[_01_3] = ( B1[q]/LP2 + B2[p]/L2 ) / (dp+dq); 

  p=2; q=0; dp=(double)p; dq=(double)q;
  I[_20_2] = ( ( (dp+1.0)*LogLP + (dq+1.0)*LogL - 1.0 ) / ((dp+1.0)*(dq+1.0))
              + A1[q] + A2[p] ) / (dp+dq+2.0);
  I[_20_3] = ( B1[q]/LP2 + B2[p]/L2 ) / (dp+dq); 

  p=1; q=1; dp=(double)p; dq=(double)q;
  I[_11_2] = ( ( (dp+1.0)*LogLP + (dq+1.0)*LogL - 1.0 ) / ((dp+1.0)*(dq+1.0))
              + A1[q] + A2[p] ) / (dp+dq+2.0);
  I[_11_3] = ( B1[q]/LP2 + B2[p]/L2 ) / (dp+dq); 

  p=0; q=2; dp=(double)p; dq=(double)q;
  I[_02_2] = ( ( (dp+1.0)*LogLP + (dq+1.0)*LogL - 1.0 ) / ((dp+1.0)*(dq+1.0))
              + A1[q] + A2[p] ) / (dp+dq+2.0);
  I[_02_3] = ( B1[q]/LP2 + B2[p]/L2 ) / (dp+dq); 

  p=2; q=1; dp=(double)p; dq=(double)q;
  I[_21_2] = ( ( (dp+1.0)*LogLP + (dq+1.0)*LogL - 1.0 ) / ((dp+1.0)*(dq+1.0))
              + A1[q] + A2[p] ) / (dp+dq+2.0);
  I[_21_3] = ( B1[q]/LP2 + B2[p]/L2 ) / (dp+dq); 

  p=1; q=2; dp=(double)p; dq=(double)q;
  I[_12_2] = ( ( (dp+1.0)*LogLP + (dq+1.0)*LogL - 1.0 ) / ((dp+1.0)*(dq+1.0))
              + A1[q] + A2[p] ) / (dp+dq+2.0);
  I[_12_3] = ( B1[q]/LP2 + B2[p]/L2 ) / (dp+dq); 

}

/*******************************************************************/
/*******************************************************************/
/*******************************************************************/
void ComputeStaticSSIData_CommonVertex(double *Xs, double *Xe, double *Xsp, double *Xep,
                                       StaticSSIDataRecord *SSSIDR)
{
  double I[NFUN];

  double R0[2], L[2], LP[2], l, lp, llp, x, CT;
  
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  VecSub(Xe,Xs,L);
  VecSub(Xep,Xsp,LP);
  VecSub(Xs,Xsp,R0);
  l=VecNorm(L);
  lp=VecNorm(LP);
  llp=l*lp;
  x=lp/l;
  CT=(L[0]*LP[0] + L[1]*LP[1])/(l*lp);

  /***************************************************************/
  /* this step amounts to switching Xs<-->Xe and/or Xep <--> Xsp */
  /* to ensure that the duffy integrals are always computed in   */
  /* the proper sense (i.e. with the singularity at (u,up)=(0,0) */
  /* instead of one of the other three corners of the unit       */
  /* square)                                                     */
  /***************************************************************/
  if ( Xs==Xsp )
   CT*=1.0;
  else if ( Xs==Xep )
   CT*=-1.0;
  else if ( Xe==Xsp )
   CT*=-1.0;
  else if ( Xe==Xep )
   CT*=1.0;
 
  /***************************************************************/
  /* fill in the easy SSSIs **************************************/
  /***************************************************************/
  SSSIDR->J_00_1=1.0;
  SSSIDR->J_10_1=1.0/2.0;
  SSSIDR->J_01_1=1.0/2.0;
  SSSIDR->J_11_1=1.0/4.0;
  SSSIDR->J_X_1=R0[0]*0.25 + (L[0]-LP[0])/6.0;
  SSSIDR->J_Y_1=R0[1]*0.25 + (L[1]-LP[1])/6.0;

  /* the other _1 quantities are finite but not needed since we  */
  /* do not need derivatives in the common-vertex case           */

  /***************************************************************/
  /* evaluate the (u,t) integrals                                */
  /***************************************************************/
  GetDuffyIntegrals(l, lp, CT, I);

  /***************************************************************/
  /* second pass through all the cases to assemble the integrals */
  /* we need from the integrals we computed                      */
  /***************************************************************/
  if ( Xs==Xsp )
   { 
     /* u up  = u up */
     /* u up (u*L - up*LP) -> L*u^2up  - LP*uup^2 */
     SSSIDR->J_00_2 = I[_00_2];
     SSSIDR->J_10_2 = I[_10_2];
     SSSIDR->J_01_2 = I[_01_2];
     SSSIDR->J_11_2 = I[_11_2];
     SSSIDR->J_X_2  = L[0]*I[_21_2] - LP[0]*I[_12_2];
     SSSIDR->J_Y_2  = L[1]*I[_21_2] - LP[1]*I[_12_2];
     SSSIDR->J_X_3  = L[0]*I[_21_3] - LP[0]*I[_12_3];
     SSSIDR->J_Y_3  = L[1]*I[_21_3] - LP[1]*I[_12_3];
   }
  else if ( Xs==Xep )
   { 
     /* u v (R0 + u*L - v*LP) -> u(1-up) (u*L + up*LP)           */
     /*                       -> L(u^2-u^2up) + LP(u*up - uup^2) */
     SSSIDR->J_00_2 = I[_00_2];
     SSSIDR->J_10_2 = I[_10_2];
     SSSIDR->J_01_2 = I[_00_2] - I[_01_2];
     SSSIDR->J_11_2 = I[_10_2] - I[_11_2];
     SSSIDR->J_X_2  = L[0]*(I[_20_2] - I[_21_2]) + LP[0]*(I[_11_2] - I[_12_2]);
     SSSIDR->J_Y_2  = L[1]*(I[_20_2] - I[_21_2]) + LP[1]*(I[_11_2] - I[_12_2]);
     SSSIDR->J_X_3  = L[0]*(I[_20_3] - I[_21_3]) + LP[0]*(I[_11_3] - I[_12_3]);
     SSSIDR->J_Y_3  = L[1]*(I[_20_3] - I[_21_3]) + LP[1]*(I[_11_3] - I[_12_3]);
   }
  else if ( Xe==Xsp )
   { 
     /* v u (R0 + v*L - up*LP) -> (1-u)up (-u*L - up*LP)          */
     /*                        -> L(u^2up-uup) + LP(u*up^2- up^2) */
     SSSIDR->J_00_2 = I[_00_2];
     SSSIDR->J_10_2 = I[_00_2] - I[_10_2];
     SSSIDR->J_01_2 = I[_01_2];
     SSSIDR->J_11_2 = I[_01_2] - I[_11_2];
     SSSIDR->J_X_2  = L[0]*(I[_21_2] - I[_11_2]) + LP[0]*(I[_12_2] - I[_02_2]);
     SSSIDR->J_Y_2  = L[1]*(I[_21_2] - I[_11_2]) + LP[1]*(I[_12_2] - I[_02_2]);
     SSSIDR->J_X_3  = L[0]*(I[_21_3] - I[_11_3]) + LP[0]*(I[_12_3] - I[_02_3]);
     SSSIDR->J_Y_3  = L[1]*(I[_21_3] - I[_11_3]) + LP[1]*(I[_12_3] - I[_02_3]);
   }
  else if ( Xe==Xep )
   { 
     /* v vp (R0 + v*L - vp*LP) -> (1-u)(1-up)(-u*L + up*LP)        */
     /*                         ->   L(u^2  - u   + uup  - u^2 up)  */
     /*                            +LP(up - uup - up^2 + uup^2)     */

     SSSIDR->J_00_2 = I[_00_2];
     SSSIDR->J_10_2 = I[_00_2] - I[_10_2];
     SSSIDR->J_01_2 = I[_00_2] - I[_01_2];
     SSSIDR->J_11_2 = I[_00_2] - I[_10_2] - I[_01_2] + I[_11_2];
     SSSIDR->J_X_2  =   L[0]*( I[_20_2] - I[_10_2] + I[_11_2] - I[_21_2] )
                      +LP[0]*( I[_01_2] - I[_11_2] - I[_02_2] + I[_12_2] );
     SSSIDR->J_Y_2  =   L[1]*( I[_20_2] - I[_10_2] + I[_11_2] - I[_21_2] )
                      +LP[1]*( I[_01_2] - I[_11_2] - I[_02_2] + I[_12_2] );
     SSSIDR->J_X_3  =   L[0]*( I[_20_3] - I[_10_3] + I[_11_3] - I[_21_3] )
                      +LP[0]*( I[_01_3] - I[_11_3] - I[_02_3] + I[_12_3] );
     SSSIDR->J_Y_3  =   L[1]*( I[_20_3] - I[_10_3] + I[_11_3] - I[_21_3] )
                      +LP[1]*( I[_01_3] - I[_11_3] - I[_02_3] + I[_12_3] );

   }
  else
   {
     fprintf(stderr,"**warning: invalid call to StaticSSIDataRecord_CommonVertex\n");
     memset(SSSIDR,0,sizeof(*SSSIDR));
   };
 
}
