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
/***************************************************************/
/***************************************************************/
static void GetReciprocalBasis(double LBVinv[2][2], 
			       double Gamma1[2], double Gamma2[2])
{
  // compute the reciprocal lattice vectors Gamma1 and Gamma2:
  // (Gamma1 Gamma2) matrix = 2*pi * inverse[transpose[(LBV[0] LBV[1])]]
  //                        == 2*pi * inverse(LBV),
  // where we denote a 2x2 matrix by (column1 column2)
  Gamma1[0] = LBVinv[0][0] * (2*M_PI);
  Gamma1[1] = LBVinv[1][0] * (2*M_PI);
  Gamma2[0] = LBVinv[0][1] * (2*M_PI);
  Gamma2[1] = LBVinv[1][1] * (2*M_PI);
}

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
/* to the Fourier-space sum that defines GBarDistant           */
/***************************************************************/
void AddGLong(double *R, cdouble k, double *P,
              double GammaX, double GammaY,
              double E, cdouble *GBarVD)
{ 
  double PmG[2];
  cdouble PreFactor, Q, EEF, EEFPrime;

  PmG[0] = P[0] - GammaX;
  PmG[1] = P[1] - GammaY;

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
void GetGBarDistant( double *R, cdouble k, double *kBloch,
                     double LBVinv[2][2],
                     double E, int *pnCells, cdouble *Sum)
{ 
  int n1, n2;
  cdouble LastSum[NSUM];
  double MaxRelDelta, MaxAbsDelta;
  double Delta, AbsSum;
  int i, NN, ConvergedIters;
  int nCells=0;

  double Gamma1[2], Gamma2[2], AGamma;

  /*--------------------------------------------------------------*/  
  /*--------------------------------------------------------------*/  
  /*--------------------------------------------------------------*/  
  GetReciprocalBasis(LBVinv, Gamma1, Gamma2);
  AGamma = Gamma1[0]*Gamma2[1] - Gamma2[0]*Gamma1[1]; // det(Gamma1 Gamma2)

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  memset(Sum,0,NSUM*sizeof(cdouble));

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  for (n1=-NFIRSTROUND; n1<=NFIRSTROUND; n1++)
   for (n2=-NFIRSTROUND; n2<=NFIRSTROUND; n2++, nCells++)
    AddGLong(R, k, kBloch,
             n1*Gamma1[0] + n2*Gamma2[0],
             n1*Gamma1[1] + n2*Gamma2[1],
             E, Sum);
         
  /***************************************************************/
  /* each iteration of this loop sums the contributions of the   */
  /* outer perimeter of an NNxNN square of grid cells.           */
  /***************************************************************/
  memcpy(LastSum,Sum,NSUM*sizeof(cdouble));
  ConvergedIters=0;
  for(NN=NFIRSTROUND+1; ConvergedIters<3 && NN<=NMAX; NN++)
   {  
     for(n1=-NN; n1<=NN; n1++)
      for(n2=-NN; n2<=NN; n2++)
       { 
         if ( (abs(n1)<NN) && (abs(n2)<NN) )
          continue;
         nCells++;
         AddGLong(R, k, kBloch,
                  n1*Gamma1[0] + n2*Gamma2[0],
                  n1*Gamma1[1] + n2*Gamma2[1], 
                  E, Sum);
       };

     /*--------------------------------------------------------------*/
     /* convergence analysis ----------------------------------------*/
     /*--------------------------------------------------------------*/
     MaxAbsDelta=MaxRelDelta=0.0;
     for(i=0; i<NSUM; i++)
      { Delta=abs(Sum[i]-LastSum[i]);
        if ( Delta>MaxAbsDelta )
         MaxAbsDelta=Delta;
        AbsSum=abs(Sum[i]);
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
  int ns; 
  for(ns=0; ns<NSUM; ns++)
   Sum[ns] *= AGamma / (16.0*M_PI*M_PI);

  if (pnCells) 
   *pnCells=nCells;

}

/***************************************************************/
/* add the contribution of a single direct lattice vector L    */
/* to the direct-lattice sum.                                  */
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
               double Lx, double Ly, double E, cdouble *Sum)
{ 
  cdouble PhaseFactor; 
  double RmL[3], rml2, rml, rml3, rml4, rml5, rml6, rml7;
  cdouble g2p, g3p, g2m, g3m, g4, ggPgg, ggMgg, Term;

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  PhaseFactor=exp( II * (kBloch[0]*Lx + kBloch[1]*Ly) ) / (8.0*M_PI);

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  RmL[0] = (R[0]-Lx);
  RmL[1] = (R[1]-Ly);
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
                   double *LBV[2], double E, bool ExcludeInnerCells,
                   int *pnCells, cdouble *Sum)
{ 
  int n1, n2;
  cdouble LastSum[NSUM];
  double MaxRelDelta, MaxAbsDelta;
  double Delta, AbsSum;
  int i, NN, ConvergedIters;

  int nCells=0;

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  memset(Sum,0,NSUM*sizeof(cdouble));

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  for (n1=-NFIRSTROUND; n1<=NFIRSTROUND; n1++)
   for (n2=-NFIRSTROUND; n2<=NFIRSTROUND; n2++, nCells++)
    if ( !ExcludeInnerCells || abs(n1)>1 || abs(n2)>1 )
     AddGShort(R, k, kBloch,
               n1*LBV[0][0] + n2*LBV[1][0],
               n1*LBV[0][1] + n2*LBV[1][1],
               E, Sum);
         
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  memcpy(LastSum,Sum,NSUM*sizeof(cdouble));
  ConvergedIters=0;
  for(NN=NFIRSTROUND+1; ConvergedIters<3 && NN<=NMAX; NN++)
   {  
     for(n1=-NN; n1<=NN; n1++)
      for(n2=-NN; n2<=NN; n2++)
       { 
         if ( (abs(n1)<NN) && (abs(n2)<NN) )
          continue;
         nCells++;
         AddGShort(R, k, kBloch,
                   n1*LBV[0][0] + n2*LBV[1][0],
                   n1*LBV[0][1] + n2*LBV[1][1], 
                   E, Sum);
       };

     /*--------------------------------------------------------------*/
     /* convergence analysis ----------------------------------------*/
     /*--------------------------------------------------------------*/
     MaxAbsDelta=MaxRelDelta=0.0;
     for(i=0; i<NSUM; i++)
      { Delta=abs(Sum[i]-LastSum[i]);
        if ( Delta>MaxAbsDelta )
         MaxAbsDelta=Delta;
        AbsSum=abs(Sum[i]);
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
                       double Lx, double Ly, double E, cdouble *Sum)
{
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  double RmL[3];
  RmL[0] = (R[0]-Lx);
  RmL[1] = (R[1]-Ly);
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

     cdouble PhaseFactor=exp( II * (kBloch[0]*Lx + kBloch[1]*Ly) );

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
     AddGFull(R, k, kBloch, Lx, Ly, GFull);

     cdouble GShort[NSUM];
     memset(GShort,0,NSUM*sizeof(cdouble));
     AddGShort(R, k, kBloch, Lx, Ly, E, GShort);

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
void GBarVDEwald(double *R, cdouble k, int LDim, double *kBloch,
                 double *LBV[2], double E, bool ExcludeInnerCells, 
                 cdouble *GBarVD)
{ 
  /*--------------------------------------------------------------*/
  /* the periodic green's function is well-defined at k==0, but   */
  /* in that case the method of this file doesn't work; moreover, */
  /* in practice the only situation in which this case arises in  */
  /* SCUFF-EM is when we have 'zeroed out' the material properties*/
  /* of a region in order to neglect its contributions to the BEM */
  /* matrix, so we want to return zeros in this case.             */
  /*--------------------------------------------------------------*/
   if (k==0.0)
    { memset(GBarVD, 0, 8*sizeof(cdouble));
      return;
    };

  if (LDim==1)
   { //GBarVD1D(R, k, kBloch, LBV[0], E, ExcludeInnerCells, GBarVD);
     ErrExit("structures with 1D periodicity not yet implemented");
     return;
   };

  double LBVinv[2][2];
  if (!Matrix2x2_Inverse(LBV, LBVinv)) 
   ErrExit("lattice has empty unit cell");

  /***************************************************************/
  /* E is the separation parameter, which we set to its optimal  */
  /* value if the user didn't specify it already                 */
  /***************************************************************/
  if (E==-1.0)
   { double EOpt1=sqrt( M_PI 
			/ fabs(LBV[0][0]*LBV[1][1] - LBV[0][1]*LBV[1][0]) );

     double Gamma1[2], Gamma2[2];
     GetReciprocalBasis(LBVinv, Gamma1, Gamma2);
     double G12 = Gamma1[0]*Gamma1[0] + Gamma1[1]*Gamma1[1];
     double G22 = Gamma2[0]*Gamma2[0] + Gamma2[1]*Gamma2[1];
     double EOpt2 = sqrt( norm(k) + G12 + G22 ) / 10.0; // H=10

     E=fmax(EOpt1, EOpt2);
   };

  /***************************************************************/
  /* evaluate 'nearby' and 'distant' sums                        */
  /***************************************************************/
  cdouble GBarNearby[NSUM], GBarDistant[NSUM];
  GetGBarNearby(R, k, kBloch, LBV, E, ExcludeInnerCells, 0, GBarNearby);
  GetGBarDistant(R, k, kBloch, LBVinv, E, 0, GBarDistant);
   
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
     for(int n1=-1; n1<=1; n1++)
      for(int n2=-1; n2<=1; n2++)
       AddGLongRealSpace(R, k, kBloch,
                         n1*LBV[0][0] + n2*LBV[1][0],
                         n1*LBV[0][1] + n2*LBV[1][1],
                         E, GLongInner);

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
  // convert R to lattice basis:
  double RL[2];
  for (int i = 0; i < 2; ++i) // RL = inv(LBV') * R; note LBV is transposed
   RL[i] = LBVinv[0][i] * R[0] + LBVinv[1][i] * R[1];
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
  GBarVDEwald(R, GBD->k, GBD->LDim, GBD->kBloch, GBD->LBV, GBD->E, 
              GBD->ExcludeInnerCells, GBarVD);
 
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
