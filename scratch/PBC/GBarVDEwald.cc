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
 * GBarVDEwald.cc -- routines for computing the periodic Green's function
 *                   and its derivatives using Ewald summation
 *
 * homer reid     -- 7/2012
 *
 */

#include <stdlib.h>
#include <math.h>

#include <libhrutil.h>
#include <libMDInterp.h>
#include "PBCGeometry.h"

#define ABSTOL 0.0
#define RELTOL 1.0e-8

#define II cdouble(0,1)

#define NSUM 8
#define NFIRSTROUND 1
#define NMAX 10000

namespace scuff{

/***************************************************************/
/* complex error function **************************************/
/***************************************************************/
extern "C" { 
 void wofz_(double *XI, double *YI, double *U, double *V, int *FLAG); 
};

cdouble cerfc(cdouble z)
{ 
  double zr, zi, wr, wi;
  int flag;

  /* handle some special cases that can cause the wofz routine to barf */
  if ( real(z)<-4.0 && fabs(imag(z))<1.0 )
   return cdouble(2.0,0.0);

  /* note the wofz_ function interprets its input as i*Z */
  zr = -imag(z);
  zi = real(z);

  wofz_(&zr, &zi, &wr, &wi, &flag);

  if (flag)
   { fprintf(stderr,"** warning: wofz(%e,%e) barfed\n",zr,zi);
     return 0.0;
   };

  return exp(-z*z) * cdouble(wr, wi);
} 

cdouble cerf(cdouble z)
{ 
  return cdouble(1.0,0.0) - cerfc(z);
} 

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
void AddG1Contribution(double *R, cdouble k, double *P,
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
void ComputeG1(double *R, cdouble k, double *P, double **LBV,
               double E, int *pnCells, cdouble *Sum)
{ 
  int n1, n2;
  double RelDelta, AbsDelta;
  cdouble LastSum[NSUM];
  double MaxRelDelta, MaxAbsDelta;
  double Delta, AbsSum;
  int i, NN, ConvergedIters;
  int nCells=0;

  double Gamma1[2], Gamma2[2], AGamma;

  /*--------------------------------------------------------------*/  
  /*--------------------------------------------------------------*/  
  /*--------------------------------------------------------------*/  
  if ( !(LBV[0][1]==0.0 && LBV[1][0]==0.0) )
   ErrExit("non-square lattices not yet supported");
  Gamma1[0] = 2.0*M_PI / LBV[0][0];
  Gamma1[1] = 0.0;
  Gamma2[0] = 0.0;
  Gamma2[1] = 2.0*M_PI / LBV[1][1];
  AGamma    = 4.0*M_PI*M_PI/(LBV[0][0]*LBV[1][1]);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  memset(Sum,0,NSUM*sizeof(cdouble));

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  for (n1=-NFIRSTROUND; n1<=NFIRSTROUND; n1++)
   for (n2=-NFIRSTROUND; n2<=NFIRSTROUND; n2++, nCells++)
    AddG1Contribution(R, k, P,
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

          AddG1Contribution(R, k, P,
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
void AddG2Contribution(double *R, cdouble k, double *P,
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
void ComputeG2(double *R, cdouble k, double *P, double **LBV,
               double E, int *pnCells, cdouble *Sum)
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
    AddG2Contribution(R, k, P,
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

          AddG2Contribution(R, k, P,
                            n1*LBV[0][0] + n2*LBV[1][0],
                            n1*LBV[0][1] + n2*LBV[1][1], 
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

  if (pnCells) *pnCells=nCells;

}

/***************************************************************/
/* get the contributions of a single real-space lattice cell to*/
/* the full periodic green's function (no ewald decomposition) */
/***************************************************************/
void AddGBFContribution(double R[3], cdouble k, double P[2],
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
void ComputeGBFFirst9(double *R, cdouble k, double *P, double *LBV[2], cdouble *Sum)
{ 
  memset(Sum,0,NSUM*sizeof(cdouble));
  for (int n1=-1; n1<=1; n1++)
   for (int n2=-1; n2<=1; n2++)
    AddGBFContribution(R, k, P,
                       n1*LBV[0][0] + n2*LBV[1][0],
                       n1*LBV[0][1] + n2*LBV[1][1],
                       Sum);
}

/***************************************************************/
/* 'GBar values and derivatives,' computed via Ewald's method  */
/*                                                             */
/* inputs:                                                     */
/*                                                             */
/*  R = 3D coordinate vector                                   */
/*  k = photon wavenumber                                      */
/*  P = 2D Bloch wavevector                                    */
/*  LBV[0][0..1] = x,y coordinates of lattice basis vector 1   */
/*  LBV[1][0..1] = x,y coordinates of lattice basis vector 2   */
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
void GBarVDEwald(double *R, cdouble k, double *P, double *LBV[2],
                 double E, int ExcludeFirst9, cdouble *GBarVD)
{ 

/***************************************************************/
/***************************************************************/
/***************************************************************/
#if 0
memset(GBarVD,0,NSUM*sizeof(cdouble));
//AddGBFContribution(R, k, P, -2.0, -2.0, GBarVD);
//return;
//AddGBFContribution(R, k, P, +2.0, +2.0, GBarVD);
{ double Lx, Ly;
  for(Lx=-2.0; Lx<=2.1; Lx+=1.0)
   { AddGBFContribution(R, k, P, Lx, +2.0, GBarVD);
     AddGBFContribution(R, k, P, Lx, -2.0, GBarVD);
   };
  for(Ly=-1.0; Ly<=1.1; Ly+=1.0)
   { AddGBFContribution(R, k, P, +2.0, Ly, GBarVD);
     AddGBFContribution(R, k, P, -2.0, Ly, GBarVD);
   };
};
return;
#endif
/***************************************************************/
/***************************************************************/
/***************************************************************/

  if ( (LBV[0][1]!=0.0) || (LBV[1][0]!=0.0) )
   ErrExit("non-square lattices not yet supported");

  /* E is the separation parameter, which we set to its  */
  /* optimal value if the user didn't specify it already */
  if (E==-1.0)
   E=sqrt( M_PI / (LBV[0][0]*LBV[1][1]) );

  /***************************************************************/
  /* FIXME *******************************************************/
  /***************************************************************/
  int OnLatticePoint=0.0;
  double MyR[3];
  memcpy(MyR,R,3*sizeof(double));
  if (    (EqualFloat( R[0], 0.0 ) || EqualFloat( R[0], LBV[0][0]))
       && (EqualFloat( R[1], 0.0 ) || EqualFloat( R[1], LBV[1][1]))
       &&  EqualFloat( R[2], 0.0 ) 
     )
   { E=0.1;
     MyR[0] = 1.0e-6;
     MyR[1] = 1.0e-6;
     MyR[2] = 1.0e-6;
   };

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  cdouble G1[NSUM], G2[NSUM], GBFFirst9[NSUM];

  ComputeG1(MyR, k, P, LBV, E, 0, G1);
  ComputeG2(MyR, k, P, LBV, E, 0, G2);

  if (ExcludeFirst9)
   ComputeGBFFirst9(MyR, k, P, LBV, GBFFirst9);
  else
   memset(GBFFirst9,0,NSUM*sizeof(cdouble));

  for(int ns=0; ns<NSUM; ns++)
   GBarVD[ns] = G1[ns] + G2[ns] - GBFFirst9[ns];

} 

/***************************************************************/
/* this is an entry point for GBarVD that has the proper       */
/* prototype for passage to the Interp3D() initialization      */
/* routine; the structure GBarData is defined in PBCGeometry.h.*/
/***************************************************************/
void GBarVDPhi3D(double X1, double X2, double X3, void *UserData, double *PhiVD)
{

  GBarData *GBD = (GBarData *)UserData;

  double R[3];
  R[0]=X1;
  R[1]=X2;
  R[2]=X3;

  cdouble GBarVD[8];
  GBarVDEwald(R, GBD->k, GBD->BlochP, GBD->LBV, GBD->E, 
              GBD->ExcludeInner9, GBarVD);
 
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
