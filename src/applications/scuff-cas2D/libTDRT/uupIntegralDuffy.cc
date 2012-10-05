/*
 * uupIntegralDuffy.cc -- evaluate segment-segment integrals using 
 *                     -- duffy-transform techniques
 *
 * homer reid          -- 10/2010
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

#include <libhrutil.h>
#include <libSGJC.h>
#include "libTDRT.h"

#define NFUN 11
#define ABSTOL 1.0e-8
#define RELTOL 1.0e-8

/*******************************************************************/
/* F[0]  = K0             (p=0, q=0)                                */ 
/* F[1]  = u*K0           (p=1, q=0)                                */
/* F[2]  = up*K0          (p=0, q=1)                                */
/* F[3]  = u*up*K0        (p=1, q=1)                                */
/* F[4]  = u*K1/R         (p=1, q=0)                                */
/* F[5]  = up*K1/R        (p=0, q=1)                                */
/* F[6]  = u^2*K1/R       (p=2, q=0)                                */
/* F[7]  = u*up*K1/R      (p=1, q=1)                                */
/* F[8]  = up^2*K1/R      (p=0, q=2)                                */
/* F[9]  = u^2*up*K1/R    (p=2, q=1)                                */
/* F[10] = u*up^2*K1/R    (p=1, q=2)                                */
/*******************************************************************/
typedef struct uupIDData
 { double l, lp, x, CT, Alpha;
 } uupIDData;

static void uupIDIntegrand(unsigned ndim, const double *ut, void *fdata, 
                           unsigned fdim, double *f)

{ 
  uupIDData *uupIDD=(uupIDData *) fdata;
  double l, lp, x, y, CT, Alpha;
  double t, t2, u, u2, u3;
  double X1, X2;

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  u=ut[0]; u2=u*u; u3=u2*u;
  t=ut[1]; t2=t*t;
   
  l=uupIDD->l;
  lp=uupIDD->lp;
  x=uupIDD->x;
  y=1.0/x;
  CT=uupIDD->CT;
  Alpha=uupIDD->Alpha;

  X1=lp*sqrt(t*t - 2.0*t*y*CT + y*y);
  X2=l*sqrt(t*t - 2.0*t*x*CT + x*x);

  if (!std::isnormal(X1)) X1=0.0;
  if (!std::isnormal(X2)) X2=0.0;

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  double KArray[2];
  double AuX1, AuX2;
  double K01, K11oX1;
  double K02, K12oX2;

  AuX1=Alpha*u*X1;
  if ( fabs(X1) < 1.0e-8 || AuX1 > 50.0 )
   { K01 = K11oX1 = 0.0; }
  else
   { HRBesselK(AuX1, 0, KArray);
     K01=KArray[0];
     K11oX1=KArray[1]*Alpha*u;
   };

  AuX2=Alpha*u*X2;
  if ( fabs(X2) < 1.0e-8 || AuX2 > 50.0 )
   { K02 = K12oX2 = 0.0; }
  else
   { 
     HRBesselK(AuX2, 0, KArray);
     K02=KArray[0];
     K12oX2=KArray[1]*Alpha*u;
   };

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  double FuX1, FuX2;

  FuX1 = K01;
  FuX2 = K02;
  f[0] =  u*(   FuX1 +   FuX2);
  f[1] = u2*(   FuX1 + t*FuX2);
  f[2] = u2*( t*FuX1 +   FuX2);
  f[3] = u3*( t*FuX1 + t*FuX2);

  FuX1 = K11oX1;
  FuX2 = K12oX2;
  f[4 ]= u*(   FuX1  +  t*FuX2);
  f[5 ]= u*( t*FuX1  +    FuX2);
  f[6 ]=u2*(   FuX1  + t2*FuX2);
  f[7 ]=u2*( t*FuX1  +  t*FuX2);
  f[8 ]=u2*(t2*FuX1  +    FuX2);
  f[9 ]=u3*( t*FuX1  + t2*FuX2);
  f[10]=u3*(t2*FuX1  +  t*FuX2);
  
}

/*******************************************************************/
/*******************************************************************/
/*******************************************************************/
void uupIntegralDuffy(double *Xs, double *Xe, double *Xsp, double *Xep,
                      double Alpha, IPQRData *IPQRD)
{
  uupIDData MyuupIDData, *uupIDD=&MyuupIDData;
  double I[NFUN], E[NFUN];
  double Lower[2]={0.0, 0.0};
  double Upper[2]={1.0, 1.0};

  double L[2], LP[2], l, lp, llp, x, CT;
  
  VecSub(Xe,Xs,L);
  VecSub(Xep,Xsp,LP);
  l=VecNorm(L);
  lp=VecNorm(LP);
  llp=l*lp;
  x=lp/l;
  CT=(L[0]*LP[0] + L[1]*LP[1])/(l*lp);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  if ( (Xs==Xsp && Xe==Xep) )
   CT*=1.0;
  else if ( Xs==Xep && Xe==Xsp )
   CT*=-1.0;
  else if ( Xs==Xsp )
   CT*=1.0;
  else if ( Xs==Xep )
   CT*=-1.0;
  else if ( Xe==Xsp )
   CT*=-1.0;
  else if ( Xe==Xep )
   CT*=1.0;

  /***************************************************************/
  /* now fill in the structure and evaluate the (u,t) integral   */
  /***************************************************************/
  uupIDD->l=l;
  uupIDD->lp=lp;
  uupIDD->x=x;
  uupIDD->CT=CT;
  uupIDD->Alpha=Alpha;
  adapt_integrate(NFUN, uupIDIntegrand, (void *)uupIDD, 2, 
                  Lower, Upper, 100000, ABSTOL, RELTOL, I, E);

  I[0]*=llp;
  I[1]*=llp;
  I[2]*=llp;
  I[3]*=llp;
  I[4]*=llp*Alpha;
  I[5]*=llp*Alpha;
  I[6]*=llp*Alpha;
  I[7]*=llp*Alpha;
  I[8]*=llp*Alpha;
  I[9]*=llp*Alpha;
  I[10]*=llp*Alpha;

  /***************************************************************/
  /* second pass through all the cases to assemble the integrals */
  /* we need from the integrals we computed                      */
  /***************************************************************/
  if ( Xs==Xsp && Xe==Xep )
   { IPQRD->I_00_0 = I[0];
     IPQRD->I_10_0 = I[1];
     IPQRD->I_01_0 = I[2];
     IPQRD->I_11_0 = I[3];
     IPQRD->I_X_1  = 0.0;
     IPQRD->I_Y_1  = 0.0;
   }
  else if ( Xs==Xep && Xe==Xsp )
   { IPQRD->I_00_0 = I[0];
     IPQRD->I_10_0 = I[1];
     IPQRD->I_01_0 = I[0] - I[2];
     IPQRD->I_11_0 = I[1] - I[3];
     IPQRD->I_X_1  = L[0]*(I[6] - I[9] - I[7] + I[10]);
     IPQRD->I_Y_1  = L[1]*(I[6] - I[9] - I[7] + I[10]);
   }
  else if ( Xs==Xsp )
   { IPQRD->I_00_0 = I[0];
     IPQRD->I_10_0 = I[1];
     IPQRD->I_01_0 = I[2];
     IPQRD->I_11_0 = I[3];
     IPQRD->I_X_1  = L[0]*I[9] - LP[0]*I[10];
     IPQRD->I_Y_1  = L[1]*I[9] - LP[1]*I[10];
   }
  else if ( Xs==Xep )
   { IPQRD->I_00_0 = I[0];
     IPQRD->I_10_0 = I[1];
     IPQRD->I_01_0 = I[0] - I[2];
     IPQRD->I_11_0 = I[1] - I[3];
     IPQRD->I_X_1  = L[0]*(I[6]-I[9]) + LP[0]*(I[7]-I[10]);
     IPQRD->I_Y_1  = L[1]*(I[6]-I[9]) + LP[1]*(I[7]-I[10]);
   }
  else if ( Xe==Xsp )
   { IPQRD->I_00_0 = I[0];
     IPQRD->I_10_0 = I[0] - I[1];
     IPQRD->I_01_0 = I[2];
     IPQRD->I_11_0 = I[2] - I[3];
     IPQRD->I_X_1  = L[0]*(I[9]-I[7]) + LP[0]*(I[10]-I[8]);
     IPQRD->I_Y_1  = L[1]*(I[9]-I[7]) + LP[1]*(I[10]-I[8]);
   }
  else if ( Xe==Xep )
   { IPQRD->I_00_0 = I[0];
     IPQRD->I_10_0 = I[0] - I[1];
     IPQRD->I_01_0 = I[0] - I[2];
     IPQRD->I_11_0 = I[0] - I[1] - I[2] + I[3];
     IPQRD->I_X_1  =  L[0] *(-I[4]+I[6]+I[7]-I[9]) 
                     +LP[0]*(I[5]-I[7]-I[8]+I[10]);
     IPQRD->I_Y_1  =  L[1] *(-I[4]+I[6]+I[7]-I[9]) 
                     +LP[1]*(I[5]-I[7]-I[8]+I[10]);
   }
  else
   {
     fprintf(stderr,"**warning: invalid call to uupIntegralDuffy\n");
     memset(IPQRD,0,sizeof(*IPQRD));
   };
 
}
