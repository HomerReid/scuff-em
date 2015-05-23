/*
 * MieSurfaceCurrents.cc
 *
 * Homer Reid    -- 1/2013
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <libhrutil.h>
#include <libSpherical.h>

#define II cdouble(0.0,1.0)
#define ZVAC 376.73031346177

/***************************************************************/
/***************************************************************/
/***************************************************************/
cdouble ipow(int n)
{ switch(n%4)
   { case 1:  return +1.0*II;
     case 2:  return -1.0;
     case 3:  return -1.0*II;
     default: return +1.0;
   };
  return 1.0;
}

/***************************************************************/
/* x = sphere radius * vacuum wavenumber                       */
/***************************************************************/
void GetRadialFunctionsAndCoefficients(cdouble x, cdouble Eps, cdouble Mu, int nMax, 
                                       cdouble *PsiX, cdouble *PsiPrimeX, 
                                       cdouble *PsiMX, cdouble *PsiPrimeMX, 
                                       cdouble *XiX, cdouble *XiPrimeX,
                                       cdouble *a, cdouble *b)
{ 
  cdouble m    = sqrt(Eps*Mu);

  /*--------------------------------------------------------------*/
  /*- fetch values of spherical bessel functions j and h ---------*/
  /*--------------------------------------------------------------*/
  cdouble jX[nMax+1], jMX[nMax+1], hX[nMax+1];
  AmosBessel('j',   x, 0, nMax+1, false,  jX);
  AmosBessel('j', m*x, 0, nMax+1, false, jMX);
  AmosBessel('o',   x, 0, nMax+1, false,  hX);

  /*--------------------------------------------------------------*/
  /*- construct the Psi and Xi functions and their derivatives    */
  /*- using d/dx (xz_l(x)) = xz_{l-1}(x) - lz_l(x)                */
  /*--------------------------------------------------------------*/
  PsiX[0]       = x*jX[0];
  PsiMX[0]      = m*x*jMX[0];
  XiX[0]        = x*hX[0];
  PsiPrimeX[0]  = cos(x);
  PsiPrimeMX[0] = cos(m*x);
  XiPrimeX[0]   = -II*exp(II*x);
  for(int n=1; n<=nMax; n++)
   { 
     double dn = (double)n;

     PsiX[n]  = x*jX[n];
     PsiMX[n] = m*x*jMX[n];
     XiX[n]   = x*hX[n];

     PsiPrimeX[n]   = PsiX[n-1]  - dn*jX[n];
     PsiPrimeMX[n]  = PsiMX[n-1] - dn*jMX[n];
     XiPrimeX[n]    = XiX[n-1]   - dn*hX[n];

     a[n] = (m*PsiMX[n]*PsiPrimeX[n] - PsiX[n]*PsiPrimeMX[n])
           /(m*PsiMX[n]*XiPrimeX[n]  - XiX[n]*PsiPrimeMX[n]);

     b[n] = (  PsiMX[n]*PsiPrimeX[n] - m*PsiX[n]*PsiPrimeMX[n])
           /(  PsiMX[n]*XiPrimeX[n]  - m*XiX[n]*PsiPrimeMX[n]);
   };
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GetTauPi(double Theta, int nMax, double *Tau, double *Pi)
{
  double Mu = cos(Theta);

  Pi[0]  = 0.0;
  Pi[1]  = 1.0;
  Tau[0] = 0.0;
  Tau[1] = Mu;

  for(int n=2; n<=nMax; n++)
   { 
     double dn=(double )n;
     Pi[n]  = (2.0*dn-1.0)*Mu*Pi[n-1]/(dn-1.0) - dn*Pi[n-2]/(dn-1.0);
     Tau[n] = dn*Mu*Pi[n] - (dn+1.0)*Pi[n-1];
   };
  
}

/***************************************************************/
/* K[0,1,2] = (r, \theta, \phi) components of K_n              */
/* N[0,1,2] = (r, \theta, \phi) components of N_n              */
/***************************************************************/
void GetMieSurfaceCurrents(double Theta, double Phi, 
                           cdouble Omega, cdouble Epsilon, cdouble Mu, 
                           int nMax, cdouble K[3], cdouble N[3])
{
  cdouble PsiX[nMax+1], PsiPrimeX[nMax+1];
  cdouble PsiMX[nMax+1], PsiPrimeMX[nMax+1];
  cdouble XiX[nMax+1], XiPrimeX[nMax+1];
  cdouble a[nMax+1], b[nMax+1];

  GetRadialFunctionsAndCoefficients(Omega, Epsilon, Mu, nMax,
                                    PsiX, PsiPrimeX, PsiMX, PsiPrimeMX, 
                                    XiX, XiPrimeX, a, b);

  double Tau[nMax+1], Pi[nMax+1];
  GetTauPi(Theta, nMax, Tau, Pi);

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  K[0]=K[1]=K[2]=N[0]=N[1]=N[2]=0.0;
  for(int n=1; n<=nMax; n++)
   { 
     cdouble En = ipow(n)*(2.0*n+1.0)/(n*(n+1.0));

     K[1] += En * (II*(PsiPrimeX[n] - b[n]*XiPrimeX[n])*Pi[n] - (PsiX[n] - a[n]*XiX[n])*Tau[n]);
     K[2] += En * ((PsiX[n] - a[n]*XiX[n])*Pi[n] - II*(PsiPrimeX[n] - b[n]*XiPrimeX[n])*Tau[n]);

     N[1] += En * (II*(PsiPrimeX[n] - a[n]*XiPrimeX[n])*Pi[n] - (PsiX[n]-b[n]*XiX[n])*Tau[n]);
     N[2] += En * ((PsiX[n] - b[n]*XiX[n])*Pi[n] - II*(PsiPrimeX[n]-a[n]*XiPrimeX[n])*Tau[n]);
   };

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  K[1] *= cos(Phi)/(Omega*ZVAC);
  K[2] *= sin(Phi)/(Omega*ZVAC);

  N[1] *= sin(Phi)/Omega;
  N[2] *= -cos(Phi)/Omega;
  
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GetMieCrossSections(double Omega, cdouble Epsilon, cdouble Mu, int nMax, 
                         double *pSigmaScat, double *pSigmaTot, double *pSigmaForce)
{ 
  cdouble PsiX[nMax+2], PsiPrimeX[nMax+2];
  cdouble PsiMX[nMax+2], PsiPrimeMX[nMax+2];
  cdouble XiX[nMax+2], XiPrimeX[nMax+2];
  cdouble a[nMax+2], b[nMax+2];
  GetRadialFunctionsAndCoefficients(Omega, Epsilon, Mu, nMax+1, 
                                    PsiX, PsiPrimeX, PsiMX, PsiPrimeMX, 
                                    XiX, XiPrimeX, a, b);

  double SigmaScat=0.0;
  double SigmaTot=0.0;
  double SigmaForce=0.0;
  for(int n=1; n<=nMax; n++)
   { double dn = (double)n;
     SigmaScat += (2.0*dn+1.0)*( norm(a[n]) + norm(b[n]) ) ; 
     SigmaTot  += (2.0*dn+1.0)*( real(a[n]) + real(b[n]) ) ; 
     SigmaForce += (dn*(dn+2.0)/(dn+1.0)) * real(   a[n]*conj(a[n+1])
                                                  + b[n]*conj(b[n+1])
                                                )
                  +(2.0*dn+1.0)/(dn*(dn+1.0)) * real( a[n]*conj(b[n]) );
   };
  SigmaScat *= 2.0*M_PI/(Omega*Omega);
  SigmaTot *= 2.0*M_PI/(Omega*Omega);
  SigmaForce*= 4.0*M_PI/(Omega*Omega);

  *pSigmaScat  = SigmaScat;
  *pSigmaTot   = SigmaTot;
  *pSigmaForce = SigmaTot - SigmaForce;

}
