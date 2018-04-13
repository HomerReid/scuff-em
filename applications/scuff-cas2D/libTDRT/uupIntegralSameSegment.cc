/*
 * uupIntegralSameSegment.cc -- evaluate segment-segment integrals in 
 *                           -- the common-segment case 
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

#define NFUN 2
#define INTERVALS 100000
#define ABSTOL 1.0e-8
#define RELTOL 1.0e-8

#define TMAX 25.0
#define LN2MGAMMA 0.115931515658412 

/*******************************************************************/
/* G_1(T) = \int_0^T  K_0( t ) dt           ************************/
/* G_2(T) = \int_0^T  t^2 K_0( t ) dt       ************************/
/*******************************************************************/
void GetG1G2(double T, double *G1, double *G2)
{ 

  if (T<0.1)
   { 
     /* use the taylor series for small T */
     double LogT, C1, C2, C3, C4, C5, T3, T5;
     LogT = (T==0.0 ? 0.0 : log(T));
     C1 = 1.0 + LN2MGAMMA - LogT;
     C2 = (1.0 + 3.0*C1) / 36.0;
     C3 = (7.0 + 10.0*C1) / (3200.0);
     C4 = (3.0*C1 - 2.0)/9.0;
     C5 = (5.0*C1 + 1.0)/100.0;
     T3=T*T*T;
     T5=T3*T*T;
     *G1= C1*T + C2*T3 + C3*T5;
     *G2= C4*T3 + C5*T5;
   }
  else if ( T > TMAX)
   { 
     /* use the T=infinity values for large T */
     *G1=*G2=M_PI_2;
   }
  else 
   { /* use interpolation tables for intermediate T */
     double G1G2[2];
     TDRTGeometry::G1G2Interp->Evaluate(T, G1G2);
     *G1 = G1G2[0];
     *G2 = G1G2[1];
   };
 
}
 
/*******************************************************************/
/* F[0]  = K0             (p=0, q=0)                               */ 
/* F[1]  = u*up*K0        (p=1, q=1)                               */
/*******************************************************************/
void uupISSIntegrand(unsigned ndim, const double *x, void *params,
                     unsigned fdim, double *f)
{ 
  (void) ndim; //unused 
  (void) fdim; //unused 


  double Alphal, AL2;
  double s;
  double G1s, G2s;

  Alphal=*((double *)params);
  AL2=Alphal*Alphal;

  s=x[0];
  GetG1G2(Alphal*s, &G1s, &G2s);

  f[0]=2.0*G1s;
  f[1]=(1.0-s+0.5*s*s)*G1s - 0.5*G2s/AL2;

}


/*******************************************************************/
/*******************************************************************/
/*******************************************************************/
void uupIntegralSameSegment(double *Xs, double *Xe, double Alpha, int Flip, 
                            IPQRData *IPQRD)
{
  double I[3], E[3];
  double l, Alphal;
  
  l=VecDistance(Xe,Xs);
  Alphal=Alpha*l;

  double Lower=0.0, Upper=1.0;
  adapt_integrate(NFUN, uupISSIntegrand, (void *)&Alphal, 1,
                  &Lower, &Upper, 0, ABSTOL, RELTOL, I, E);

  I[0]*=l/Alpha;
  I[1]*=l/Alpha;

  if (Flip)
   {
     IPQRD->I_00_0=I[0];
     IPQRD->I_10_0=0.5*I[0];
     IPQRD->I_01_0=0.5*I[0];
     IPQRD->I_11_0=0.5*I[0] - I[1];
   }
  else
   { 
     IPQRD->I_00_0=I[0];
     IPQRD->I_10_0=0.5*I[0];
     IPQRD->I_01_0=0.5*I[0];
     IPQRD->I_11_0=I[1];
   };

  IPQRD->I_X_1=0.0; 
  IPQRD->I_Y_1=0.0; 

}   
