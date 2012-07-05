/*
 * GBarVD[0] = GBar
 * GBarVD[1] = dGBar/dX
 * GBarVD[2] = dGBar/dY
 * GBarVD[3] = dGBar/dZ
 * GBarVD[4] = d^2GBar/dXdY
 * GBarVD[5] = d^2GBar/dXdZ
 * GBarVD[6] = d^2GBar/dYdZ
 * GBarVD[7] = d^3GBar/dXdYdZ
 */
#include <stdlib.h>
#include <math.h>

#include "libhrutil.h"
#include "libMDInterp.h"

#define ABSTOL 0.0
#define RELTOL 1.0e-8

#define II cdouble(0,1)

#define NSUM 8
#define NFIRSTROUND 1
#define NMAX 10000

int RetainFirst9=0;

/***************************************************************/
/***************************************************************/
/***************************************************************/
cdouble cerfc(cdouble z);

void GBarVDBF(cdouble k, double *P, double *L1, double *L2, double *R,
              double AbsTol, double RelTol, int *pnCells, cdouble *GBarVD);

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

  // PlusTerm  = exp(  kz*R[2] ) * cerfc( 0.5*kz/E  + R[2]*E );
  ExpFac    = exp( Q*z );
  dExpFac   = Q*ExpFac;
  Arg       = 0.5*Q/E + z*E;
  ErfcFac   = cerfc( Arg );
  dErfcFac  = -E*exp( -Arg*Arg );
  PlusTerm  = ExpFac*ErfcFac;
  dPlusTerm = dExpFac*ErfcFac + ExpFac*dErfcFac;

  // MinusTerm  = exp( -kz*R[2] ) * cerfc( 0.5*kz/E  - R[2]*E );
  ExpFac     = exp( -Q*z );
  dExpFac    = -Q*ExpFac;
  Arg        = 0.5*Q/E - z*E;
  ErfcFac    = cerfc( Arg );
  dErfcFac   = +E*exp( -Arg*Arg );
  MinusTerm  = ExpFac*ErfcFac;
  dMinusTerm = dExpFac*ErfcFac + ExpFac*dErfcFac;
  
  *EEF      = PlusTerm + MinusTerm;
  *EEFPrime = dPlusTerm + dMinusTerm;

  if ( !isfinite(*EEF) )
   *EEF=0.0;
  if ( !isfinite(*EEFPrime) )
   *EEFPrime=0.0;

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void AddG1Contribution(cdouble k, double *P, double *R, 
                       double GammaX, double GammaY,
                       double E, cdouble *GBarVD)
{ 
  double GmP[2];
  cdouble PreFactor, Q, EEF, EEFPrime;

  GmP[0] = GammaX - P[0];
  GmP[1] = GammaY - P[1];

  Q = sqrt ( GmP[0]*GmP[0] + GmP[1]*GmP[1] - k*k );

  PreFactor = exp( II * (GmP[0]*R[0] + GmP[1]*R[1]) ) / Q;

  GetEEF(R[2], E, Q, &EEF, &EEFPrime);

  GBarVD[0] += PreFactor * EEF;
  GBarVD[1] += II*GmP[0]*PreFactor*EEF;
  GBarVD[2] += II*GmP[1]*PreFactor*EEF;
  GBarVD[3] += PreFactor*EEFPrime;
  GBarVD[4] += -GmP[0]*GmP[1]*PreFactor*EEF;
  GBarVD[5] += II*GmP[0]*PreFactor*EEFPrime;
  GBarVD[6] += II*GmP[1]*PreFactor*EEFPrime;
  GBarVD[7] += -GmP[0]*GmP[1]*PreFactor*EEFPrime;

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void ComputeG1(cdouble k, double *P, double *L1, double *L2, 
               double *R, double E, int *pnCells, cdouble *Sum)
{ 
  int n1, n2;
  double RelDelta, AbsDelta;
  cdouble LastSum[NSUM];
  double MaxRelDelta, MaxAbsDelta;
  double Delta, AbsSum;
  int i, NN, ConvergedIters;
  int nCells=0;

  double Gamma1[2], Gamma2[2], AGamma;

  if ( !(L1[1]==0.0 && L2[0]==0.0) )
   ErrExit("non-square lattices not supported");
  Gamma1[0] = 2.0*M_PI / L1[0];
  Gamma1[1] = 0.0;
  Gamma2[0] = 0.0;
  Gamma2[1] = 2.0*M_PI / L2[1];
  AGamma    = 4.0*M_PI*M_PI/(L1[0]*L2[1]);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  memset(Sum,0,NSUM*sizeof(cdouble));

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  for (n1=-NFIRSTROUND; n1<=NFIRSTROUND; n1++)
   for (n2=-NFIRSTROUND; n2<=NFIRSTROUND; n2++, nCells++)
    AddG1Contribution(k, P, R,
                      n1*Gamma1[0] + n2*Gamma2[0],
                      n1*Gamma1[1] + n2*Gamma2[1],
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

          AddG1Contribution(k, P, R, 
                            n1*Gamma1[0] + n2*Gamma2[0],
                            n1*Gamma1[1] + n2*Gamma2[1], 
                            E, Sum);

          nCells++;

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
/*  exp(i P\dot L) * g_1 *( g_{2+} g_{3+} + g_{2-}g_{3-})      */  
/*                                                             */
/* where:                                                      */
/*  g_1       = 1/(8*pi*|R|)                                   */
/*  g_{2+}    = exp(i*k*|R|)                                   */
/*  g_{3+}    = erfc( E*R + i*k/(2*E) )                        */
/*  g_{2-}    = exp(-i*k*|R|)                                  */
/*  g_{3-}    = erfc( E*R - i*k/(2*E) )                        */
/* (here |R| = |r+L|.)                                         */
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
/* where g4 = (-4E/sqrt(pi)) * exp( -E^2R^2 + k^2/(4E^2).      */
/*                                                             */
/***************************************************************/
void AddG2Contribution(cdouble k, double *P, double *R, 
                       double Lx, double Ly, double E, cdouble *Sum)
{ 
  cdouble PhaseFactor; 
  double RpL[3], rpl2, rpl, rpl3, rpl4, rpl5, rpl6, rpl7;
  cdouble g2p, g3p, g2m, g3m, g4, ggPgg, ggMgg, Term;

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  PhaseFactor=exp( II * (P[0]*Lx + P[1]*Ly) ) / (8.0*M_PI);

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  RpL[0] = (R[0]+Lx);
  RpL[1] = (R[1]+Ly);
  RpL[2] =  R[2];

  rpl2=RpL[0]*RpL[0] + RpL[1]*RpL[1] + RpL[2]*RpL[2];
  rpl=sqrt(rpl2);
  if ( rpl < 1.0e-7 ) 
   return;
  rpl3=rpl2*rpl;
  rpl4=rpl3*rpl;
  rpl5=rpl4*rpl;
  rpl6=rpl5*rpl;
  rpl7=rpl6*rpl;

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  g2p = exp( II*k*rpl );
  g3p = cerfc( E*rpl + II*k/(2.0*E) );

  g2m = exp( -II*k*rpl );
  g3m = cerfc( E*rpl - II*k/(2.0*E) );

  g4 = -2.0*M_2_SQRTPI*E*exp(-E*E*rpl2 + k*k/(4.0*E*E));

  ggPgg = g2p*g3p + g2m*g3m;
  ggMgg = g2p*g3p - g2m*g3m;

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  Sum[0] += PhaseFactor * ggPgg / rpl;

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  Term = -ggPgg/rpl3 + (g4 + II*k*ggMgg)/rpl2;

  Sum[1] += PhaseFactor * RpL[0] * Term;
  Sum[2] += PhaseFactor * RpL[1] * Term;
  Sum[3] += PhaseFactor * RpL[2] * Term;

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  Term = 3.0*ggPgg/rpl5 - 3.0*(g4+II*k*ggMgg)/rpl4 - k*k*ggPgg/rpl3 - 2.0*E*E*g4/rpl2;

  Sum[4] += PhaseFactor * RpL[0] * RpL[1] * Term;
  Sum[5] += PhaseFactor * RpL[0] * RpL[2] * Term;
  Sum[6] += PhaseFactor * RpL[1] * RpL[2] * Term;

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  Term = -15.0*ggPgg/rpl7 + 15.0*(g4+II*k*ggMgg)/rpl6 
         + 6.0*k*k*ggPgg/rpl5 + 10.0*E*E*g4/rpl4
         -k*k*(II*k*ggMgg + g4)/rpl4 + 4.0*E*E*E*E*g4/rpl2;

  Sum[7] += PhaseFactor * RpL[0] * RpL[1] * RpL[2] * Term;

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void ComputeG2(cdouble k, double *P, double *L1, double *L2,
               double *R, double E, int *pnCells, cdouble *Sum)
{ 
  int n1, n2;
  double RelDelta, AbsDelta;
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
    { 
#if 0
      if ( (abs(n1)<=1) && (abs(n2)<=1) )
       continue; // skip the innermost 9 grid cells 
#endif

      AddG2Contribution(k, P, R, 
                        n1*L1[0] + n2*L2[0],
                        n1*L1[1] + n2*L2[1], 
                        E, Sum);
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
//printf("%i %i %s \n",n1,n2,CD2S(Sum[0]));
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
    };
         
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

          AddG2Contribution(k, P, R, 
                            n1*L1[0] + n2*L2[0],
                            n1*L1[1] + n2*L2[1], 
                            E, Sum);
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
//printf("%i %i %s \n",n1,n2,CD2S(Sum[0]));
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/

          nCells++;
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
/***************************************************************/
/* compute the contribution of the innermost 9 real-space      */
/* cells to the fourier-space sum                              */
/***************************************************************/
/***************************************************************/
#if 0 // this seems not to have worked
void ComputeG1First9(cdouble k, double *P, double *L1, double *L2, 
                     double *R, double E, cdouble *Sum)
{
  int n1, n2;
  double L[2], RpL[3], rpl2, rpl;
  cdouble PhaseFactor, IKR, Phi, Psi, Zeta, Upsilon;
  cdouble Part1[8], Part2[8];

  memset(Part1,0,NSUM*sizeof(cdouble));
  memset(Part2,0,NSUM*sizeof(cdouble));
  for(n1=-1; n1<=1; n1++)
   for(n2=-1; n2<=1; n2++)
    { 
      L[0] = n1*L1[0] + n2*L2[0];
      L[1] = n1*L1[1] + n2*L2[1];

      PhaseFactor = exp( II*(P[0]*L[0] + P[1]*L[1]) ) / (4.0*M_PI);

      RpL[0] = R[0]-L[0];
      RpL[1] = R[1]-L[1];
      RpL[2] = R[2];

      rpl2=RpL[0]*RpL[0] + RpL[1]*RpL[1] + RpL[2]*RpL[2];
      rpl=sqrt(rpl2);
      if (fabs(rpl)<1.0e-7) 
       continue;

      IKR=II*k*rpl;
      Phi=exp(IKR) / (4.0*M_PI*rpl);
      Psi=(-1.0 + IKR) * Phi / rpl2;
      Zeta=(3.0 + IKR*(-3.0 + IKR))*Phi/(rpl2*rpl2);
      Upsilon=(-15.0 + IKR*(15.0 + IKR*(-6.0 + IKR)))*Phi/(rpl2*rpl2*rpl2);
 
      Part1[0] += PhaseFactor * Phi;
      Part1[1] += PhaseFactor * RpL[0] * Psi;
      Part1[2] += PhaseFactor * RpL[1] * Psi;
      Part1[3] += PhaseFactor * RpL[2] * Psi;
      Part1[4] += PhaseFactor * RpL[0] * RpL[1] * Zeta;
      Part1[5] += PhaseFactor * RpL[0] * RpL[2] * Zeta;
      Part1[6] += PhaseFactor * RpL[1] * RpL[2] * Zeta;
      Part1[7] += PhaseFactor * RpL[1] * RpL[2] * RpL[3] * Zeta;

      AddG2Contribution(k, P, R, L[0], L[1], E, Part2);

    };

  int ns;
  for(ns=0; ns<NSUM; ns++)
   Sum[ns] = Part1[ns] - 0.5*Part2[ns];

};
#endif

/***************************************************************/
/* get the contributions ***************************************/
/***************************************************************/
void AddGBFContribution(cdouble k, double P[2], double R[3],
                        double Lx, double Ly, cdouble *Sum)
{ 

  double RpL[3];
  double r2, r, kr; 
  cdouble PhaseFactor, IKR, Phi, Psi, Zeta, Upsilon;
   
  PhaseFactor=exp( II*(Lx*P[0] + Ly*P[1]) );

  RpL[0]=R[0] + Lx;
  RpL[1]=R[1] + Ly;
  RpL[2]=R[2];

  r2=RpL[0]*RpL[0] + RpL[1]*RpL[1] + RpL[2]*RpL[2];
  r=sqrt(r2);
  if (r<1.0e-7) 
   return;
  IKR=II*k*r;
  Phi=exp(IKR)/(4.0*M_PI*r);
  Psi=(IKR-1.0)*Phi/r2;
  Zeta=(3.0 + IKR*(-3.0 + IKR))*Phi/(r2*r2);
  Upsilon=(-15.0 + IKR*(15.0 + IKR*(-6.0 + IKR)))*Phi/(r2*r2*r2);

  Sum[0] += PhaseFactor * Phi;
  Sum[1] += PhaseFactor * RpL[0] * Psi;
  Sum[2] += PhaseFactor * RpL[1] * Psi;
  Sum[3] += PhaseFactor * RpL[2] * Psi;
  Sum[4] += PhaseFactor * RpL[0] * RpL[1] * Zeta;
  Sum[5] += PhaseFactor * RpL[0] * RpL[2] * Zeta;
  Sum[6] += PhaseFactor * RpL[1] * RpL[2] * Zeta;
  Sum[7] += PhaseFactor * RpL[0] * RpL[1] * RpL[2] * Upsilon;
 
}

/***************************************************************/
/* sum the contributions to the innermost 9 real-space cells   */
/***************************************************************/
void ComputeGBFFirst9(cdouble k, double *P, double *L1, double *L2,
                      double *R, cdouble *Sum)
{ 
  int n1, n2;

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  memset(Sum,0,NSUM*sizeof(cdouble));
  for (n1=-1; n1<=1; n1++)
   for (n2=-1; n2<=1; n2++)
    AddGBFContribution(k, P, R, 
                       n1*L1[0] + n2*L2[0],
                       n1*L1[1] + n2*L2[1], 
                       Sum);
}

/***************************************************************/
/* 'GBar values and derivatives,' computed via Ewald's method  */
/*                                                             */
/* inputs:                                                     */
/*                                                             */
/*  k = photon wavenumber                                      */
/*  P = 2D Bloch wavevector                                    */
/*  L1, L2 = 2D lattice vectors                                */
/*  R = 3D coordinate vector                                   */
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
void GBarVDEwald(cdouble k, double *P, double *L1, double *L2, 
                 double *R, cdouble *GBarVD, double E)
{ 
  
  if ( (L1[1]!=0.0) || (L2[0]!=0.0) ) 
   ErrExit("%s:%i: internal error",__FILE__,__LINE__);

  /* E is the separation parameter, which we set to its  */
  /* optimal value if the user didn't specify it already */
  if (E==-1.0)
   E=sqrt( M_PI / (L1[0]*L2[1]) );

  /***************************************************************/
  /* FIXME *******************************************************/
  /***************************************************************/
  int OnLatticePoint=0.0;
  if ( fabs(R[2]) < 1.0e-7
      && ( (fabs(R[0]) < 1.0e-7) || (fabs(fabs(R[0])-1.0) < 1.0e-7) )
      && ( (fabs(R[1]) < 1.0e-7) || (fabs(fabs(R[1])-1.0) < 1.0e-7) )
     ) E=0.1;

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  cdouble G1[NSUM], G2[NSUM], GBFFirst9[NSUM];

  ComputeG1(k, P, L1, L2, R, E, 0, G1);
  ComputeG2(k, P, L1, L2, R, E, 0, G2);
if (RetainFirst9) 
  memset(GBFFirst9,0,NSUM*sizeof(double));
else
  ComputeGBFFirst9(k, P, L1, L2, R, GBFFirst9);

  int ns;
  for(ns=0; ns<NSUM; ns++)
   GBarVD[ns] = G1[ns] + G2[ns] - GBFFirst9[ns];

} 
