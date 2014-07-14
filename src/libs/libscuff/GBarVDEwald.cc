#include <stdlib.h>
#include <math.h>

#include <libhrutil.h>
#include <libMDInterp.h>
#include <libscuff.h>
#include <libscuffInternals.h>

#define ABSTOL 0.0
#define RELTOL 1.0e-8

#define II cdouble(0,1)

#define NSUM 8
#define NFIRSTROUND 1
#define NMAX 10000

#include "Faddeeva.hh"

using namespace Faddeeva;
namespace scuff{

/***************************************************************/
/* compute exp(a)*erfc(b), being careful of overflow/underflow */
/***************************************************************/
static cdouble erfc_s(cdouble a, cdouble b) 
{
  double x = real(b); 
  double y = imag(b);

  cdouble mb2( (y - x) * (x + y), -2*x*y); // -b^2

  if (x >= 0) // erfc(b) = exp(-b^2) erfcx(b) ==> exp(a-b^2) erfcx(b)
   return exp(a + mb2) * Faddeeva::erfcx(b);
  else        // erfc(b) = 2 - exp(-b^2) erfcx(b) ==> 2*exp(a) - exp(a-b^2) erfcx(-b)
   return 2.0*exp(a) - exp(a + mb2) * Faddeeva::erfcx(-b);
}

static bool cisnan(cdouble z) 
 { return ISNAN(real(z)) || ISNAN(imag(z)); }

/***************************************************************/
/* 'EEF' = 'exp*erfc factor'  **********************************/
/*                                                             */
/*   EEF = e^{Q*Z} * erfc[ Q/2E + E*Z]                         */
/*          + e^{-Q*Z} * erfc[ Q/2E - E*Z]                     */
/***************************************************************/
void GetEEF(double z, double E, cdouble Q, cdouble *EEF, cdouble *EEFPrime)
{ 
  cdouble Arg, ExpFac, dExpFac, ErfcFac, dErfcFac;
  cdouble PlusTerm, dPlusTerm, MinusTerm, dMinusTerm;

  // PlusTerm  = exp(  kz*R[2] ) * erfc( 0.5*kz/E  + R[2]*E );
  Arg       = 0.5*Q/E + z*E;
  PlusTerm = erfc_s(Q*z, Arg);
  dPlusTerm = Q*PlusTerm - (M_2_SQRTPI * E) * exp(Q*z - Arg*Arg);

  // MinusTerm  = exp( -kz*R[2] ) * erfc( 0.5*kz/E  - R[2]*E );
  Arg        = 0.5*Q/E - z*E;
  MinusTerm = erfc_s(-Q*z, Arg);
  dMinusTerm = -Q*MinusTerm + (M_2_SQRTPI * E) * exp(-Q*z - Arg*Arg);
  
  *EEF      = PlusTerm + MinusTerm;
  *EEFPrime = dPlusTerm + dMinusTerm;

  if ( !isfinite(*EEF) || cisnan(*EEF) )
   *EEF=0.0;
  if ( !isfinite(*EEFPrime) || cisnan(*EEFPrime) )
   *EEFPrime=0.0;

}

/***************************************************************/
/* add the contribution of a single reciprocal-lattice vector  */
/* G = n1*Gamma1 + n2*Gamma2 to the reciprocal-lattice sum     */
/* that defines GBarDistant.                                   */
/***************************************************************/
void AddGLong(double R[3], cdouble k, double P[2],
              int n1, int n2, double *Gamma[2], int LDim,
              double E, cdouble *GBarVD)
{ 
  double PmG[2];
  cdouble PreFactor, Q, EEF, EEFPrime;
   
  PmG[0] = P[0] - n1*Gamma[0][0] - n2*Gamma[1][0];
  PmG[1] = P[1] - n1*Gamma[1][0] - n2*Gamma[1][1];

  Q = sqrt ( PmG[0]*PmG[0] + PmG[1]*PmG[1] - k*k );

  PreFactor = exp( II * (PmG[0]*R[0] + PmG[1]*R[1]) ) / Q;

  GetEEF(R[2], E, Q, &EEF, &EEFPrime);

  GBarVD[0] += PreFactor * EEF;
  GBarVD[1] += II*PmG[0]*PreFactor*EEF;
  GBarVD[2] += II*PmG[1]*PreFactor*EEF;
  GBarVD[3] += PreFactor*EEFPrime;
  GBarVD[4] += -PmG[0]*PmG[1]*PreFactor*EEF;
  GBarVD[5] += II*PmG[0]*PreFactor*EEFPrime;
  GBarVD[6] += II*PmG[1]*PreFactor*EEFPrime;
  GBarVD[7] += -PmG[0]*PmG[1]*PreFactor*EEFPrime;

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GetGBarDistant(double *R, cdouble k, double *kBloch,
                    double *Gamma[2], int LDim,
                    double E, int *pnCells, cdouble *Sum)
{ 
  /***************************************************************/
  /* start by summing the contributions of a ``first round''     */
  /* of cells near the origin                                    */
  /***************************************************************/
  memset(Sum,0,NSUM*sizeof(cdouble));
  int n2Mult = (LDim==2) ? 1 : 0;
  int nCells=0;
  for (int n1=-NFIRSTROUND; n1<=NFIRSTROUND; n1++)
   for (int n2=-NFIRSTROUND*n2Mult; n2<=NFIRSTROUND*n2Mult; n2++, nCells++)
    AddGLong(R, k, kBloch, n1, n2, Gamma, LDim, E, Sum);

  /***************************************************************/
  /* continue to add contributions of outer cells until converged*/
  /***************************************************************/
  cdouble LastSum[NSUM];
  memcpy(LastSum,Sum,NSUM*sizeof(cdouble));
  int ConvergedIters=0;
  for(int NN=NFIRSTROUND+1; ConvergedIters<3 && NN<=NMAX; NN++)
   {  
     if (LDim==1)
      { AddGLong(R, k, kBloch,  NN, 0, Gamma, LDim, E, Sum);
        AddGLong(R, k, kBloch, -NN, 0, Gamma, LDim, E, Sum);
        nCells+=2;
      }
     else // LDim==2
      { 
        /*--------------------------------------------------------------*/
        /* sum the contributions of the outer perimeter of an           */
        /* NNxNN square of grid cells.                                  */
        /*--------------------------------------------------------------*/
        for(int n=-NN; n<NN; n++)
         { AddGLong(R, k, kBloch,   n,  NN, Gamma, LDim, E, Sum);
           AddGLong(R, k, kBloch,  NN,  -n, Gamma, LDim, E, Sum);
           AddGLong(R, k, kBloch,  -n, -NN, Gamma, LDim, E, Sum);
           AddGLong(R, k, kBloch, -NN,   n, Gamma, LDim, E, Sum);
           nCells+=4;
         };
      };

     /*--------------------------------------------------------------*/
     /* convergence analysis ----------------------------------------*/
     /*--------------------------------------------------------------*/
     double MaxRelDelta=0.0, MaxAbsDelta=0.0;
     for(int ns=0; ns<NSUM; ns++)
      { double Delta=abs(Sum[ns]-LastSum[ns]);
        if ( Delta>MaxAbsDelta )
         MaxAbsDelta=Delta;
        double AbsSum=abs(Sum[ns]);
        if ( AbsSum>0.0 && (Delta > MaxRelDelta*AbsSum) )
         MaxRelDelta=Delta/AbsSum;
      };
     if ( MaxAbsDelta<ABSTOL || MaxRelDelta<RELTOL )
      ConvergedIters++;
     else
      ConvergedIters=0;

     memcpy(LastSum,Sum,NSUM*sizeof(cdouble));

   };

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  double VGamma;
  if (LDim==1)
   VGamma = sqrt(Gamma[0][0]*Gamma[0][0] + Gamma[0][1]*Gamma[0][1]);
  else
   VGamma = Gamma[0][0]*Gamma[1][1] - Gamma[0][1]*Gamma[1][0];

  for(int ns=0; ns<NSUM; ns++)
   Sum[ns] *= VGamma / (16.0*M_PI*M_PI);

  if (pnCells) 
   *pnCells=nCells;

}

/***************************************************************/
/* add the contribution of a single direct lattice vector L    */
/* to the direct-lattice sum, where L = n1*L1 + n2*L2.         */
/*                                                             */
/* note: the summand is:                                       */
/*                                                             */
/*  exp(i kBloch\dot L) * g_1*( g_{2+} g_{3+} + g_{2-}g_{3-} ) */  
/*                                                             */
/* where:                                                      */
/*                                                             */
/*  g_1       = 1/(8*pi*|R|)                                   */
/*  g_{2+}    = exp(i*k*|R|)                                   */
/*  g_{3+}    = erfc(  E*R + i*k/(2*E) )                       */
/*  g_{2-}    = exp(-i*k*|R|)                                  */
/*  g_{3-}    = erfc(  E*R - i*k/(2*E)  )                      */
/*                                                             */
/* with |R| = |r-L|.                                           */ 
/*                                                             */
/* if we define the following combinations:                    */
/*                                                             */
/*  ggPgg = g_{2+} g_{3+}  +  g_{2-} g_{3-}                    */
/*  ggMgg = g_{2+} g_{3+}  -  g_{2-} g_{3-}                    */
/*                                                             */
/* then we have the following derivative relationships:        */
/*                                                             */
/* d/dR (ggPgg) = i*k*ggMgg + g4                               */
/* d/dR (ggMgg) = i*k*ggPgg                                    */
/*                                                             */
/* where g4 = (-4E/sqrt(pi)) * exp( -(E)^2R^2 + k^2/(4(E)^2).  */
/*                                                             */
/***************************************************************/
void AddGShort(double *R, cdouble k, double *kBloch,
               int n1, int n2, double *LBV[2], int LDim,
               double E, cdouble *Sum)
{ 
  cdouble PhaseFactor; 
  double RmL[3], rml2, rml, rml3, rml4, rml5, rml6, rml7;
  cdouble g2p, g3p, g2m, g3m, g4, ggPgg, ggMgg, Term;

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  double L[2];
  if (LDim==1)
   { L[0] = n1*LBV[0][0];
     L[1] = n1*LBV[0][1];
   }
  else // (LDim==2)
   { L[0] = n1*LBV[0][0] + n2*LBV[1][0];
     L[1] = n1*LBV[0][1] + n2*LBV[1][1];
   };
  PhaseFactor=exp( II * (kBloch[0]*L[0] + kBloch[1]*L[1]) ) / (8.0*M_PI);

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  RmL[0] = (R[0]-L[0]);
  RmL[1] = (R[1]-L[1]);
  RmL[2] =  R[2];

  rml2=RmL[0]*RmL[0] + RmL[1]*RmL[1] + RmL[2]*RmL[2];
  rml=sqrt(rml2);
  if ( rml < 1.0e-6 ) 
   return;
  rml3=rml2*rml;
  rml4=rml3*rml;
  rml5=rml4*rml;
  rml6=rml5*rml;
  rml7=rml6*rml;

  double E2 = E*E;
  double E4 = E2*E2;

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  //g2p = exp( II*k*rml );
  //g3p = Faddeeva::erfc( E*rml + II*k/(2.0*E) );
  cdouble g2pTg3p = erfc_s( II*k*rml, E*rml + II*k/(2.0*E) );

  //g2m = exp( -II*k*rml );
  //g3m = Faddeeva::erfc( E*rml - II*k/(2.0*E) );
  cdouble g2mTg3m = erfc_s( -II*k*rml, E*rml - II*k/(2.0*E) );

  //ggPgg = g2p*g3p + g2m*g3m;
  //ggMgg = g2p*g3p - g2m*g3m;
  ggPgg = g2pTg3p + g2mTg3m;
  ggMgg = g2pTg3p - g2mTg3m;

  g4 = -2.0*M_2_SQRTPI*E*exp(-E2*rml2 + k*k/(4.0*E2));

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  Sum[0] += PhaseFactor * ggPgg / rml;

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  Term = -ggPgg/rml3 + (g4 + II*k*ggMgg)/rml2;

  Sum[1] += PhaseFactor * RmL[0] * Term;
  Sum[2] += PhaseFactor * RmL[1] * Term;
  Sum[3] += PhaseFactor * RmL[2] * Term;

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  Term = 3.0*ggPgg/rml5 - 3.0*(g4+II*k*ggMgg)/rml4 
           - k*k*ggPgg/rml3 - 2.0*E2*g4/rml2;

  Sum[4] += PhaseFactor * RmL[0] * RmL[1] * Term;
  Sum[5] += PhaseFactor * RmL[0] * RmL[2] * Term;
  Sum[6] += PhaseFactor * RmL[1] * RmL[2] * Term;

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  Term = -15.0*ggPgg/rml7 + 15.0*(g4+II*k*ggMgg)/rml6 
        + 6.0*k*k*ggPgg/rml5 + 10.0*E2*g4/rml4
         -k*k*(II*k*ggMgg + g4)/rml4 + 4.0*E4*g4/rml2;

  Sum[7] += PhaseFactor * RmL[0] * RmL[1] * RmL[2] * Term;

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GetGBarNearby(double *R, cdouble k, double *kBloch,
                   double *LBV[2], int LDim,
                   double E, bool ExcludeInnerCells,
                   int *pnCells, cdouble *Sum)
{ 
  /***************************************************************/
  /* add the contributions of a 'first round' of grid cells near */
  /* the center                                                  */
  /***************************************************************/
  int nCells=0;
  memset(Sum,0,NSUM*sizeof(cdouble));
  if (LDim==1)
   { 
     for (int n1=-NFIRSTROUND; n1<=NFIRSTROUND; n1++, nCells++)
       if ( !ExcludeInnerCells || abs(n1)>1 )
        AddGShort(R, k, kBloch, n1, 0, LBV, LDim, E, Sum);
   }
  else if (LDim==2)
   { 
     for (int n1=-NFIRSTROUND; n1<=NFIRSTROUND; n1++)
      for (int n2=-NFIRSTROUND; n2<=NFIRSTROUND; n2++, nCells++)
       if ( !ExcludeInnerCells || abs(n1)>1 || abs(n2)>1 )
        AddGShort(R, k, kBloch, n1, n2, LBV, LDim, E, Sum);
   };
         
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  cdouble LastSum[NSUM];
  memcpy(LastSum,Sum,NSUM*sizeof(cdouble));
  int ConvergedIters=0;
  for(int NN=NFIRSTROUND+1; ConvergedIters<3 && NN<=NMAX; NN++)
   {  
     if (LDim==1)
      { AddGShort(R, k, kBloch,  NN, 0, LBV, LDim, E, Sum);
        AddGShort(R, k, kBloch, -NN, 0, LBV, LDim, E, Sum);
        nCells+=2;
      }
     else // LDim==2
      { 
        /*--------------------------------------------------------------*/
        /* sum the contributions of the outer perimeter of an */
        /* NNxNN square of grid cells.                        */
        /*--------------------------------------------------------------*/
        for(int n=-NN; n<NN; n++)
         { AddGShort(R, k, kBloch,   n,  NN, LBV, LDim, E, Sum);
           AddGShort(R, k, kBloch,  NN,  -n, LBV, LDim, E, Sum);
           AddGShort(R, k, kBloch,  -n, -NN, LBV, LDim, E, Sum);
           AddGShort(R, k, kBloch, -NN,   n, LBV, LDim, E, Sum);
           nCells+=4;
         };
      };

     /*--------------------------------------------------------------*/
     /* convergence analysis ----------------------------------------*/
     /*--------------------------------------------------------------*/
     double MaxAbsDelta=0.0, MaxRelDelta=0.0;
     for(int ns=0; ns<NSUM; ns++)
      { double Delta=abs(Sum[ns]-LastSum[ns]);
        if ( Delta>MaxAbsDelta )
         MaxAbsDelta=Delta;
        double AbsSum=abs(Sum[ns]);
        if ( AbsSum>0.0 && (Delta > MaxRelDelta*AbsSum) )
         MaxRelDelta=Delta/AbsSum;
      };
     if ( MaxAbsDelta<ABSTOL || MaxRelDelta<RELTOL )
      ConvergedIters++;
     else
      ConvergedIters=0;

     memcpy(LastSum,Sum,NSUM*sizeof(cdouble));

   };

  if (pnCells) *pnCells=nCells;

}

/***************************************************************/
/* get the contributions of a single real-space lattice cell to*/
/* the full periodic green's function (no ewald decomposition) */
/***************************************************************/
void AddGFull(double R[3], cdouble k, double kBloch[2],
              double Lx, double Ly, cdouble *Sum)
{ 

  double RmL[3];
  double r2, r;
  cdouble PhaseFactor, IKR, Phi, Psi, Zeta, Upsilon;
   
  PhaseFactor=exp( II*(Lx*kBloch[0] + Ly*kBloch[1]) );

  RmL[0]=R[0] - Lx;
  RmL[1]=R[1] - Ly;
  RmL[2]=R[2];

  r2=RmL[0]*RmL[0] + RmL[1]*RmL[1] + RmL[2]*RmL[2];
  r=sqrt(r2);
  if ( r < 1.0e-8 )
   return;
  IKR=II*k*r;
  Phi=exp(IKR)/(4.0*M_PI*r);
  Psi=(IKR-1.0)*Phi/r2;
  Zeta=(3.0 + IKR*(-3.0 + IKR))*Phi/(r2*r2);
  Upsilon=(-15.0 + IKR*(15.0 + IKR*(-6.0 + IKR)))*Phi/(r2*r2*r2);

  Sum[0] += PhaseFactor * Phi;
  Sum[1] += PhaseFactor * RmL[0] * Psi;
  Sum[2] += PhaseFactor * RmL[1] * Psi;
  Sum[3] += PhaseFactor * RmL[2] * Psi;
  Sum[4] += PhaseFactor * RmL[0] * RmL[1] * Zeta;
  Sum[5] += PhaseFactor * RmL[0] * RmL[2] * Zeta;
  Sum[6] += PhaseFactor * RmL[1] * RmL[2] * Zeta;
  Sum[7] += PhaseFactor * RmL[0] * RmL[1] * RmL[2] * Upsilon;
 
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
#define PI32 5.5683279968317078453 // pi^{3/2}
void AddGLongRealSpace(double *R, cdouble k, double *kBloch,
                       int n1, int n2, double *LBV[2], int LDim, 
                       double E, cdouble *Sum)
{
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  double L[2];
  if (LDim==1)
   { L[0] = n1*LBV[0][0];
     L[1] = n1*LBV[0][1];
   }
  else 
   { L[0] = n1*LBV[0][0] + n2*LBV[1][0];
     L[1] = n1*LBV[0][1] + n2*LBV[1][1];
   };

  double RmL[3];
  RmL[0] = (R[0]-L[0]);
  RmL[1] = (R[1]-L[1]);
  RmL[2] =  R[2];

  double rml2=RmL[0]*RmL[0] + RmL[1]*RmL[1] + RmL[2]*RmL[2];
  double rml=sqrt(rml2);
  bool Smallr = ( rml*E < 0.5 ) && ( rml*abs(k) < 1.0 );

  if ( Smallr ) // use small-r expansion derived in memo
   { 
     cdouble k2  = k*k;
     cdouble k3  = k2*k;
     cdouble E2  = E*E;
     cdouble ErfFac = 1.0 + erf(0.5*II*k/E);
     cdouble ExpFac = exp(0.25*k2/E2);
     cdouble C0  = E*ExpFac/(2.0*PI32) + II*k*ErfFac/(4.0*M_PI);
     cdouble C2  = -ExpFac*E*(2.0*E2+k2)/(12.0*PI32) 
                    -II*k3*ErfFac/(24.0*M_PI);
     cdouble C4  = ExpFac*E*(12.0*E2*E2 + 2.0*E2*k2 + k2*k2)/(240.0*PI32)
                    +II*k3*k2*ErfFac/(480.0*M_PI);

     cdouble PhaseFactor=exp(II * (kBloch[0]*L[0] + kBloch[1]*L[1]) );

     Sum[0] += PhaseFactor*(C0 + C2*rml2 + C4*rml2*rml2);
     Sum[1] += PhaseFactor * (2.0*C2*RmL[0] + 4.0*C4*rml2*RmL[0]);
     Sum[2] += PhaseFactor * (2.0*C2*RmL[1] + 4.0*C4*rml2*RmL[1]);
     Sum[3] += PhaseFactor * (2.0*C2*RmL[2] + 4.0*C4*rml2*RmL[2]);
     Sum[4] += PhaseFactor * (8.0*C4*RmL[0]*RmL[1]);
     Sum[5] += PhaseFactor * (8.0*C4*RmL[0]*RmL[2]);
     Sum[6] += PhaseFactor * (8.0*C4*RmL[1]*RmL[2]);
     Sum[7] += 0.0; // third mixed partial vanishes in this approximation
   }
  else // set GLong = GFull - GShort
   { 
     cdouble GFull[NSUM];
     memset(GFull,0,NSUM*sizeof(cdouble));
     AddGFull(R, k, kBloch, L[0], L[1], GFull);

     cdouble GShort[NSUM];
     memset(GShort,0,NSUM*sizeof(cdouble));
     AddGShort(R, k, kBloch, n1, n2, LBV, LDim, E, GShort);

     for(int ns=0; ns<NSUM; ns++)
      { 
        cdouble GLong = GFull[ns] - GShort[ns];
        Sum[ns] += GLong;

        if ( ns==0 && abs(GLong) < 1.0e-6*(abs(GFull[ns]) + abs(GShort[ns])) )
         Warn("loss of precision (r=%e) ( %.8e - %.8e = %.1e ) in GLongRealSpace",
               rml,abs(GFull[ns]),abs(GShort[ns]),abs(GLong));
      };
   };
 
}

/***************************************************************/
/* 'GBar values and derivatives,' computed via Ewald's method  */
/*                                                             */
/* inputs:                                                     */
/*                                                             */
/*  R = 3D coordinate vector                                   */
/*  k = photon wavenumber                                      */
/*  LDim = lattice periodicity dimension                       */
/*        (LDim=1 or 2 for 1D or 2D periodicity)               */
/*                                                             */
/*  kBloch[0,1]  = x,y components of bloch wavevector          */
/*  LBV[0][0..1] = x,y coordinates of lattice basis vector 1   */
/*  LBV[1][0..1] = x,y coordinates of lattice basis vector 2   */
/*                                                             */
/*  (Note kBloch[1] and LBV[1] are only referenced if LDim==2) */
/*                                                             */
/*  E = ewald separation parameter (set to -1.0 for automatic) */
/*                                                             */
/* outputs:                                                    */
/*                                                             */
/*  GBarVD[0] = GBar                                           */
/*  GBarVD[1] = dGBar/dX                                       */
/*  GBarVD[2] = dGBar/dY                                       */
/*  GBarVD[3] = dGBar/dZ                                       */
/*  GBarVD[4] = d^2GBar/dXdY                                   */
/*  GBarVD[5] = d^2GBar/dXdZ                                   */
/*  GBarVD[6] = d^2GBar/dYdZ                                   */
/*  GBarVD[7] = d^3GBar/dXdYdZ                                 */
/*                                                             */
/***************************************************************/
void GBarVDEwald(double *R, cdouble k, double *kBloch,
                 double *LBV[2], int LDim,
                 double E, bool ExcludeInnerCells,
                 cdouble *GBarVD)
{ 
  /*--------------------------------------------------------------*/
  /* the periodic green's function is well-defined at k==0 (i.e.  */
  /* the electrostatic case), but in that case the method used    */
  /* in this file doesn't work; moreover,                         */
  /* in practice the only situation in which this case arises in  */
  /* SCUFF-EM is when we have 'zeroed out' the material properties*/
  /* of a region in order to neglect its contributions to the BEM */
  /* matrix, so we want to return zeros in this case.             */
  /*--------------------------------------------------------------*/
   if (k==0.0)
    { memset(GBarVD, 0, 8*sizeof(cdouble));
      return;
    };
 
  /*--------------------------------------------------------------*/
  /*- get a basis for the reciprocal lattice and choose the      -*/
  /*- separation parameter eta if the user didn't choose it      -*/
  /*--------------------------------------------------------------*/
  double Gamma1[2], Gamma2[2];
  double *Gamma[2];
  Gamma[0]=Gamma1;
  Gamma[1]=Gamma2;
  if (LDim==1)
   { 
     double Length = sqrt( LBV[0][0]*LBV[0][0] + LBV[0][1]*LBV[0][1]) ;
     Gamma[0][0] = (2.0*M_PI/Length) * LBV[0][0];
     Gamma[0][1] = (2.0*M_PI/Length) * LBV[0][1];
     Gamma[1][0] = Gamma[1][1] = 0.0;

     if (E==-1.0)
      E = sqrt(M_PI) / Length;
   }
  else if (LDim==2)
   { 
     double Area= LBV[0][0]*LBV[1][1] - LBV[0][1]*LBV[1][0];
     if (Area==0.0)
      ErrExit("%s:%i: lattice has empty unit cell",__FILE__,__LINE__);
     Gamma[0][0] =  2.0*M_PI*LBV[1][1] / Area;
     Gamma[0][1] = -2.0*M_PI*LBV[0][1] / Area;
     Gamma[1][0] = -2.0*M_PI*LBV[1][0] / Area;
     Gamma[1][1] =  2.0*M_PI*LBV[0][0] / Area;

     if (E==-1.0)
      { double EOpt1 = sqrt(M_PI / Area);
        double G12 = Gamma[0][0]*Gamma[0][0] + Gamma[0][1]*Gamma[0][1];
        double G22 = Gamma[1][0]*Gamma[1][0] + Gamma[1][1]*Gamma[1][1];
        double EOpt2 = sqrt( norm(k) + G12 + G22 );
        E = fmax(EOpt1, EOpt2);
      };
   }
  else
   ErrExit("only 1D or 2D periodicity implemented in GBarVDEwald");

  /***************************************************************/
  /* evaluate 'nearby' and 'distant' sums                        */
  /***************************************************************/
  cdouble GBarNearby[NSUM], GBarDistant[NSUM];
  GetGBarNearby(R, k, kBloch, LBV, LDim,
                E, ExcludeInnerCells, 0, GBarNearby);
  GetGBarDistant(R, k, kBloch, Gamma, LDim, E, 0, GBarDistant);
   
  for(int ns=0; ns<NSUM; ns++)
   GBarVD[ns] = GBarNearby[ns] + GBarDistant[ns];

  /***************************************************************/
  /* subtract off the contributions to the 'distant' sum coming  */
  /* from the inner grid cells in real space                     */
  /***************************************************************/
  if (ExcludeInnerCells)
   { 
     cdouble GLongInner[NSUM];
     memset(GLongInner,0,NSUM*sizeof(cdouble));
     int n2Mult = (LDim==2) ? 1 : 0;
     for(int n1=-1; n1<=1; n1++)
      for(int n2=-1*n2Mult; n2<=1*n2Mult; n2++)
       AddGLongRealSpace(R, k, kBloch, n1, n2, LBV, LDim, E, GLongInner);
     for(int ns=0; ns<NSUM; ns++)
      GBarVD[ns] -= GLongInner[ns];
   };

  /***************************************************************/
  /* detect evaluation points at the origin or lattice-equivalent*/
  /* to the origin so that we can explicitly zero out the        */
  /* appropriate derivatives in this case. NOTE 20140425: This   */
  /* step is perhaps not needed now that we have the improved    */
  /* treatment of evaluation points on lattice sites.            */
  /***************************************************************/
  // convert R to lattice basis, RL = inv(LBV') * R; note LBV is transposed
  double RL[2];
  RL[0] = (Gamma[0][0]*R[0] + Gamma[1][0]*R[1]) / (2.0*M_PI);
  if (LDim==2)
   RL[1] = (Gamma[0][1]*R[0] + Gamma[1][1]*R[1]) / (2.0*M_PI);
  else
   RL[1] = 2.0;

  bool ZeroCoordinate[3]={false, false, false};
  const double tol = 1e-8;
  ZeroCoordinate[0] = fabs(RL[0]) < tol || fabs( (fabs(RL[0]) - 1.0) ) < tol;
  ZeroCoordinate[1] = fabs(RL[1]) < tol || fabs( (fabs(RL[1]) - 1.0) ) < tol;
  ZeroCoordinate[2] = fabs(R[2]) < tol;

  if ( ZeroCoordinate[0] )
   GBarVD[1]=GBarVD[4]=GBarVD[5]=GBarVD[7]=0.0;
  if ( ZeroCoordinate[1] )
   GBarVD[2]=GBarVD[4]=GBarVD[6]=GBarVD[7]=0.0;
  if ( ZeroCoordinate[2] )
   GBarVD[3]=GBarVD[5]=GBarVD[6]=GBarVD[7]=0.0;


} 

/***************************************************************/
/* this is an entry point for GBarVD that has the proper       */
/* prototype for passage to the Interp3D() initialization      */
/* routine; the structure GBarData is defined in libscuffInternals.h.*/
/***************************************************************/
void GBarVDPhi3D(double X1, double X2, double X3, void *UserData, double *PhiVD)
{
  GBarData *GBD = (GBarData *)UserData;

  double R[3];
  R[0]=X1;
  R[1]=X2;
  R[2]=X3;

  cdouble GBarVD[8];
  GBarVDEwald(R, GBD->k, GBD->kBloch, GBD->LBV, GBD->LDim,
              GBD->E, GBD->ExcludeInnerCells, GBarVD);
 
  PhiVD[ 0] = real(GBarVD[0]);
  PhiVD[ 1] = real(GBarVD[1]);
  PhiVD[ 2] = real(GBarVD[2]);
  PhiVD[ 3] = real(GBarVD[3]);
  PhiVD[ 4] = real(GBarVD[4]);
  PhiVD[ 5] = real(GBarVD[5]);
  PhiVD[ 6] = real(GBarVD[6]);
  PhiVD[ 7] = real(GBarVD[7]);
  PhiVD[ 8] = imag(GBarVD[0]);
  PhiVD[ 9] = imag(GBarVD[1]);
  PhiVD[10] = imag(GBarVD[2]);
  PhiVD[11] = imag(GBarVD[3]);
  PhiVD[12] = imag(GBarVD[4]);
  PhiVD[13] = imag(GBarVD[5]);
  PhiVD[14] = imag(GBarVD[6]);
  PhiVD[15] = imag(GBarVD[7]);

} 

} // namespace scuff
