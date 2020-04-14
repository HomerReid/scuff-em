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
#include <stdlib.h>
#include <stdlib.h>
#include <math.h>

#include <libhrutil.h>
#include <libMDInterp.h>
#include <libSpherical.h>
#include <libscuff.h>
#include <libscuffInternals.h>

#define ABSTOL 0.0
#define RELTOL 1.0e-8

#define II cdouble(0,1)

#define NSUM 8
#define NFIRSTROUND 1
#define NMAX 10000

#define EULERGAMMA 0.57721566490153286

#include "Faddeeva.hh"

using namespace Faddeeva;
namespace scuff{

/*****************************************************/
/* Continued fraction for ExpIntegral[1,z] using     */
/* modified Lentz's method.                          */
/*****************************************************/
#define TINY 1.0e-30
cdouble ExpInt_CF(cdouble z, double RelTol, int *nIters)
{
  // initialization
  cdouble aj=0.0, bj=0.0;
  cdouble fjm1, fj = TINY;
  cdouble Cjm1, Cj = fj;
  cdouble Djm1, Dj = 0.0;
  int ConvergedIters=0;
  cdouble Delta;
  int j;
  for(j=1; j<1000; j++)
    {   
     if (j==1)
      { aj=exp(-z); bj=z; }
     else
      { aj = j/2; bj = (j%2) ? z : 1.0; }

     Djm1 = Dj;
     Dj = bj + aj*Djm1; 
     if (Dj==0.0) Dj = TINY;

     Cjm1 = Cj;
     Cj = bj + aj/Cjm1;
     if (Cj==0.0) Cj = TINY;

     Dj = 1.0/Dj;
     Delta = Cj*Dj;

     fjm1 = fj;
     fj = fjm1*Delta;

     if ( abs(Delta-1.0) < RelTol )
      ConvergedIters++;
     else
      ConvergedIters=0;
     if (ConvergedIters==5)
      break;
   };
  if (ConvergedIters!=5)
   Warn("potentially large error in ExpInt_CF(%s) [%.1e %%]",
         z2s(z), 100.0*abs(Delta-1.0));

  if (nIters) *nIters=j;
  return fj;

}

/*****************************************************/
/* power series for ExpIntegralE[1,z] ****************/
/*****************************************************/
cdouble ExpInt_PS(cdouble z, double RelTol=1.0e-8, int *pnTerms=0)
{ 
  cdouble Result = -EULERGAMMA - log(z);

  double nFact=1.0;
  cdouble zPower=1.0;
  cdouble LastResult;
  int n;
  for(n=1; n<=100; n++)
   { 
     nFact *= n;
     zPower *= -z;
     LastResult = Result;
     Result -= zPower/(nFact*n);
     if ( abs(Result-LastResult) < RelTol*abs(Result) )
      break;
   };
  if (n==100)
   Warn("potentially large error in ExpInt_PS(%s) [%.1e %%]",
         CD2S(z), 100.0*abs(Result-LastResult)/abs(Result));

  if (pnTerms) *pnTerms=n;

  return Result;
  
}

/*****************************************************/
/* large-z power series for ExpIntegralE[1,z] ********/
/*****************************************************/
cdouble ExpInt_Asymptotic(cdouble z, double RelTol=1.0e-8, int *pnTerms=0)
{ 
  if ( abs(z)>100.0 ) 
   return 0.0;

  cdouble Sum  = 1.0;
  cdouble Term = 1.0;
  int n;
  for(n=1; n<=100; n++)
   { 
     Term *= -((double)n) / z;
     Sum += Term;
     if ( abs(Term) < RelTol*abs(Sum) )
      break;
   };
  if (n==100)
   Warn("potentially large error in ExpInt_PS(%s) [%.1e %%]",
         CD2S(z), 100.0*abs(Term)/abs(Sum));

  if (pnTerms) *pnTerms=n;

  return Sum * exp(-z) / z;
  
}

/***************************************************************/
/* exponential integral E[1,x]                                 */
/***************************************************************/
cdouble ExpInt(cdouble z)
{
  double absz = abs(z); 

  if ( absz < 5.0 )
   return ExpInt_PS(z, 1.0e-8, 0); 
  else if ( absz < 30.0 )
   { if ( fabs( fabs(arg(z))-M_PI ) < 0.2 )
      return ExpInt_PS(z, 1.0e-8, 0);
     else
      return ExpInt_CF(z, 1.0e-8, 0);
   }
  else 
   return ExpInt_Asymptotic(z, 1.0e-8, 0);
}

/***************************************************************/
/* Given a 1D or 2D basis for a direct lattice, compute a      */
/* basis for the reciprocal lattice and the 1D or 2D "volume"  */
/* (i.e. length or area) of the Brillouin zone.                */
/* If EOpt is non-null, then the optimal value of the Ewald    */
/* separation parameter Eta is also computed and returned in   */
/* EOpt. (The photon wavenumber k is used in this computation.)*/
/*                                                             */
/* In the 1D case, if pRho is nonzero then on return it        */
/* contains the quantity Rho, the perpendicular distance from  */
/* R to the line defined by the lattice vector.                */
/***************************************************************/
void GetRLBasis(int LDim, double (*L)[3], double (*Gamma)[3],
                cdouble k, double *EOpt, double R[3], double *pRho)
{
  /*--------------------------------------------------------------*/
  /*- get a basis for the reciprocal lattice and choose the      -*/
  /*- separation parameter eta if the user didn't choose it      -*/
  /*--------------------------------------------------------------*/
  if (LDim==1)
   { 
     double L2 = L[0][0]*L[0][0] + L[0][1]*L[0][1];
     Gamma[0][0] = (2.0*M_PI/L2) * L[0][0];
     Gamma[0][1] = (2.0*M_PI/L2) * L[0][1];
     Gamma[1][0] = Gamma[1][1] = 0.0;
     
     double Factor = (R[0]*L[0][0] + R[1]*L[0][1]) / L2;
     double Rho2D[2];
     Rho2D[0] = R[0] - Factor*L[0][0];
     Rho2D[1] = R[1] - Factor*L[0][1];
     double Rho2 = Rho2D[0]*Rho2D[0] + Rho2D[1]*Rho2D[1] + R[2]*R[2];
     double Rho=sqrt(Rho2);
     if (pRho) *pRho=Rho;

     /***************************************************************/
     /* choose optimal E parameter following the discussion of      */
     /* Valerio et al, IEEE TAP *55* 1630 (2007)                    */
     /***************************************************************/
     double E  = sqrt(M_PI / L2);
     double E1 = abs(k) / 20.0; // H=10
     double E2 = (Rho2==0.0) ? 1.0e100 : 1.2 / Rho;
     if ( E < E1 ) 
      E=E1;
     else if ( E > E2 )
      E=E2;
     *EOpt = E;
   }
  else if (LDim==2)
   { 
     double Area= L[0][0]*L[1][1] - L[0][1]*L[1][0];
     if (Area==0.0)
      ErrExit("%s:%i: lattice has empty unit cell",__FILE__,__LINE__);
     Gamma[0][0] =  2.0*M_PI*L[1][1] / Area;
     Gamma[0][1] = -2.0*M_PI*L[0][1] / Area;
     Gamma[1][0] = -2.0*M_PI*L[1][0] / Area;
     Gamma[1][1] =  2.0*M_PI*L[0][0] / Area;

     if (EOpt)
      { double EOpt1 = sqrt(M_PI / Area);
        double G12 = Gamma[0][0]*Gamma[0][0] + Gamma[0][1]*Gamma[0][1];
        double G22 = Gamma[1][0]*Gamma[1][0] + Gamma[1][1]*Gamma[1][1];
        double EOpt2 = sqrt( norm(k) + G12 + G22 ) / 10.0; // H=10
        *EOpt = fmax(EOpt1, EOpt2);
      };
   }
  else
   ErrExit("only 1D or 2D periodicity implemented in GBarVDEwald");
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

  if ( !IsFinite(*EEF) || cisnan(*EEF) )
   *EEF=0.0;
  if ( !IsFinite(*EEFPrime) || cisnan(*EEFPrime) )
   *EEFPrime=0.0;

}

/***************************************************************/
/* add the contribution of a single reciprocal-lattice vector  */
/* G = n1*Gamma1 + n2*Gamma2 to the reciprocal-lattice sum     */
/* that defines GBarDistant in the 2D case.                    */
/***************************************************************/
void AddGLong2D(double R[3], cdouble k, double P[2],
                int n1, int n2, double Gamma[3][3],
                double E, cdouble *GBarVD, bool &Singular)
{ 
  if (Singular) return;

  double PmG[2];
  cdouble PreFactor, Q, EEF, EEFPrime;
   
  PmG[0] = P[0] - n1*Gamma[0][0] - n2*Gamma[1][0];
  PmG[1] = P[1] - n1*Gamma[0][1] - n2*Gamma[1][1];

  Q = sqrt ( PmG[0]*PmG[0] + PmG[1]*PmG[1] - k*k );

  if ( abs(Q) < 1.0e-4*abs(k) )
   { Singular=true;
     return;
   };

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

/*****************************************************/
/* Compute the 1D Fourier transform of GFull (the    */
/* full Helmholtz Green's function).                 */
/*****************************************************/
cdouble GFullTwiddle1D(double kx, double Rho, cdouble k,
                       cdouble *dGdRho=0)
{
  cdouble kt2 = kx*kx - k*k;
  cdouble kt = sqrt(kt2);
  cdouble K[3];
  AmosBessel('K',kt*Rho,0.0,3,false,K,0);
  double Denom = 4.0*M_PI*M_PI;

  if (dGdRho)
   { dGdRho[0] = -kt*K[1] / Denom;
     dGdRho[1] = kt2*(K[0] + K[2])/(2.0*Denom);
  };
  
  return K[0] / Denom;
}

/*****************************************************/
/* Compute the 1D Fourier transform of GLong.        */
/* If dGdRho is non-null, then on return we have     */
/* dGdRho[0] = dG   / dRho                           */
/* dGdRho[1] = dG^2 / dRho^2                         */
/*****************************************************/
cdouble GetGLongTwiddle1D(double kx, double Rho, cdouble k, double E,
                          cdouble *dGdRho, bool &Singular)
{
  if (Singular)
   return 0.0;

  if ( Rho*E > 4.5 )
   return GFullTwiddle1D(kx, Rho, k, dGdRho);

  if (dGdRho) 
   dGdRho[0]=dGdRho[1]=0.0;

  cdouble kt2      = kx*kx - k*k;
  if ( abs(kt2) < 1.0e-8 * abs(k*k) )
   { Singular=true;
     return 0.0;
   };

  double E2        = E*E;
  cdouble Arg      = kt2 / (4.0*E2);
  cdouble Eqp1     = ExpInt(Arg);
  double NormFac   = 8.0*M_PI*M_PI;

  if (Rho==0.0)
   return Eqp1 / NormFac;

  cdouble ExpFac   = exp(-Arg);
  double RhoE2     = Rho*Rho*E2;
  double PreFactor = 1.0;
  cdouble Sum      = Eqp1;
  int ConvergedIters=0;
  for(int q=1; q<1000; q++)
   {
     // use recurrence to get ExpIntegral_{q+1}
     Eqp1 = (ExpFac - Arg*Eqp1) / ((double)q);

     PreFactor *= -1.0*RhoE2 / ((double)q);

     cdouble Summand = PreFactor*Eqp1;
     Sum += Summand;

     if (dGdRho)
      { double Factor = 2.0*((double)q)/Rho;
        dGdRho[0] += Factor*Summand;
        dGdRho[1] += Factor*Summand*(2.0*((double)q)-1.0)/Rho;
      };

     /* monitor convergence of value sum only; assume */
     /* derivative sums converge at same rate         */
     if ( abs(Summand/Sum) < 1.0e-8 )
      ConvergedIters++;
     else
      ConvergedIters=0;
     if (ConvergedIters==2) break;

   };
  if (ConvergedIters!=2 && abs(Sum)>1.0e-8 )
   Warn("potential nonconvergence in GLongTwiddle1D(%g,%g,%s,%e)",
         kx, Rho, z2s(k), E);

  if (dGdRho)
   { dGdRho[0] /= NormFac;
     dGdRho[1] /= NormFac;
   };
  return Sum / NormFac;
  
}

// 
// 
void AddGLong1D(double R[3], double Rho, cdouble k, double P[2],
                int m, double Gamma[3][3], double E, cdouble *GBarVD,
                bool &Singular)
{
  if (Singular) return;

  double PmG[2], PmGMag;
  PmG[0] = P[0] - m*Gamma[0][0];
  PmG[1] = P[1] - m*Gamma[0][1];
  if (PmG[1]!=0.0)
   ErrExit("1D lattice vectors must point in the x direction");
  PmGMag = PmG[0];

  cdouble ExpFac = exp( II * ( PmG[0]*R[0] + PmG[1]*R[1]) );

  cdouble dGdRho[2];
  cdouble GT=GetGLongTwiddle1D(PmGMag, Rho, k, E, dGdRho, Singular);
  cdouble dGTdRho = dGdRho[0];
  cdouble dGT2dRho2 = Rho==0.0 ? 0.0 : (dGdRho[1] - dGdRho[0]/Rho);
  double YOverRho = (Rho==0.0) ? 0.0 : R[1]/Rho;
  double ZOverRho = (Rho==0.0) ? 0.0 : R[2]/Rho;

  GBarVD[0] += ExpFac * GT;
  GBarVD[1] += II*PmG[0] * ExpFac * GT;
  GBarVD[2] += YOverRho * ExpFac * dGTdRho;
  GBarVD[3] += ZOverRho * ExpFac * dGTdRho;
  GBarVD[4] += II*PmG[0] * YOverRho * ExpFac * dGTdRho;
  GBarVD[5] += II*PmG[0] * ZOverRho * ExpFac * dGTdRho;
  GBarVD[6] += YOverRho * ZOverRho * ExpFac * dGT2dRho2;
  GBarVD[7] += II*PmG[0] * YOverRho * ZOverRho * ExpFac * dGT2dRho2;

}

/***************************************************************/
/* set Singular=true if singularity encountered ****************/
/***************************************************************/
bool GetGBarDistant(double *R, double Rho, cdouble k, double *kBloch,
                    double Gamma[3][3], int LDim,
                    double E, int *pnCells, cdouble *Sum)
{
  memset(Sum,0,NSUM*sizeof(cdouble));
  if (E==0.0) return false;

  bool Singular=false;

  /***************************************************************/
  /* start by summing the contributions of a ``first round''     */
  /* of cells near the origin                                    */
  /***************************************************************/
  int nCells=0;
  if (LDim==1)
   { for (int m=-NFIRSTROUND; m<=NFIRSTROUND; m++)
      AddGLong1D(R, Rho, k, kBloch, m, Gamma, E, Sum, Singular);
   }
  else // LDim==2
   { for (int m1=-NFIRSTROUND; m1<=NFIRSTROUND; m1++)
      for (int m2=-NFIRSTROUND; m2<=NFIRSTROUND; m2++, nCells++)
       AddGLong2D(R, k, kBloch, m1, m2, Gamma, E, Sum, Singular);
   };

  if (Singular) return true;

  /***************************************************************/
  /* continue to add contributions of outer cells until converged*/
  /***************************************************************/
  cdouble LastSum[NSUM];
  memcpy(LastSum,Sum,NSUM*sizeof(cdouble));
  int ConvergedIters=0;
  for(int NN=NFIRSTROUND+1; ConvergedIters<3 && NN<=NMAX; NN++)
   {  
     if (LDim==1)
      { AddGLong1D(R, Rho, k, kBloch,  NN, Gamma, E, Sum, Singular);
        AddGLong1D(R, Rho, k, kBloch, -NN, Gamma, E, Sum, Singular);
        nCells+=2;
      }
     else // LDim==2
      { 
        /*--------------------------------------------------------------*/
        /* sum the contributions of the outer perimeter of the innermost*/
        /* NNxNN square of grid cells                                   */
        /*--------------------------------------------------------------*/
        for(int m=-NN; m<NN; m++)
         { AddGLong2D(R, k, kBloch,   m,  NN, Gamma, E, Sum, Singular);
           AddGLong2D(R, k, kBloch,  NN,  -m, Gamma, E, Sum, Singular);
           AddGLong2D(R, k, kBloch,  -m, -NN, Gamma, E, Sum, Singular);
           AddGLong2D(R, k, kBloch, -NN,   m, Gamma, E, Sum, Singular);
           nCells+=4;
         };
      };

     if (Singular) return true;

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
  double PreFactor;
  if (LDim==1)
   PreFactor = sqrt(Gamma[0][0]*Gamma[0][0] + Gamma[0][1]*Gamma[0][1]);
  else
   PreFactor = (Gamma[0][0]*Gamma[1][1] - Gamma[0][1]*Gamma[1][0])/(16.0*M_PI*M_PI);

  for(int ns=0; ns<NSUM; ns++)
   Sum[ns] *= PreFactor;

  if (pnCells) 
   *pnCells=nCells;

  return false;

}

/***************************************************************/
/* get the contributions of a single real-space lattice cell to*/
/* the FULL periodic green's function (no ewald decomposition).*/
/* If ValueOnly==true, derivative calculation is skipped.      */
/***************************************************************/
void AddGFull(double R[3], cdouble k, double kBloch[2],
              double Lx, double Ly, cdouble *Sum, 
              bool ValueOnly=false)
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
  Sum[0] += PhaseFactor * Phi;

  if (ValueOnly) return;

  Psi=(IKR-1.0)*Phi/r2;
  Zeta=(3.0 + IKR*(-3.0 + IKR))*Phi/(r2*r2);
  Upsilon=(-15.0 + IKR*(15.0 + IKR*(-6.0 + IKR)))*Phi/(r2*r2*r2);

  Sum[1] += PhaseFactor * RmL[0] * Psi;
  Sum[2] += PhaseFactor * RmL[1] * Psi;
  Sum[3] += PhaseFactor * RmL[2] * Psi;
  Sum[4] += PhaseFactor * RmL[0] * RmL[1] * Zeta;
  Sum[5] += PhaseFactor * RmL[0] * RmL[2] * Zeta;
  Sum[6] += PhaseFactor * RmL[1] * RmL[2] * Zeta;
  Sum[7] += PhaseFactor * RmL[0] * RmL[1] * RmL[2] * Upsilon;
 
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
               int n1, int n2, double (*LBV)[3], int LDim,
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
  else // (LDim==2)
   { L[0] = n1*LBV[0][0] + n2*LBV[1][0];
     L[1] = n1*LBV[0][1] + n2*LBV[1][1];
   };

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  if (E==0.0)
   { AddGFull(R, k, kBloch, L[0], L[1], Sum);
     return;
   };

  cdouble PhaseFactor=exp( II * (kBloch[0]*L[0] + kBloch[1]*L[1]) ) / (8.0*M_PI);

  double RmL[3], rml2, rml, rml3, rml4, rml5, rml6, rml7;
  cdouble g2p, g3p, g2m, g3m, g4, ggPgg, ggMgg, Term;

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
                   double (*LBV)[3], int LDim,
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
        /* sum the contributions of the outer perimeter of the innermost*/
        /* NNxNN square of grid cells.                                  */
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
/***************************************************************/
/***************************************************************/
#define PI32 5.5683279968317078453 // pi^{3/2}
void AddGLongRealSpace(double *R, cdouble k, double *kBloch,
                       int n1, int n2, double (*LBV)[3], int LDim, 
                       double E, cdouble *Sum)
{  
  if (E==0.0) return;

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
  bool Smallr = ( rml*E < 0.1 ) && ( rml*abs(k) < 1.0 );

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

        if (    (ns==0)
             && (abs(GLong) < 1.0e-6*(abs(GFull[ns]) + abs(GShort[ns])))
             && (abs(GLong) > 1.0e-8)
           )
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
/***************************************************************/
void GBarVDEwald(double *R, cdouble k, double *kBloch0,
                 double (*LBV)[3], int LDim,
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

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  double Gamma[3][3], EOpt, Rho;
  GetRLBasis(LDim, LBV, Gamma, k, &EOpt, R, &Rho);
  if (E==-1.0) E=EOpt;

  /***************************************************************/
  /* evaluate 'nearby' and 'distant' sums                        */
  /***************************************************************/
  double kBloch[3]={0.0, 0.0, 0.0};
  memcpy(kBloch, kBloch0, LDim*sizeof(double));

  cdouble GBarNearby[NSUM], GBarDistant[NSUM];
  GetGBarNearby(R, k, kBloch, LBV, LDim,
                E, ExcludeInnerCells, 0, GBarNearby);
  bool Singular
   = GetGBarDistant(R, Rho, k, kBloch, Gamma, LDim, E, 0, GBarDistant);

  if (Singular)
   { Log("Ewald spectral sum is singular (recomputing at displaced kBloch)");
     kBloch[0] += (1.0e-2)*abs(k);
     GetGBarNearby(R, k, kBloch, LBV, LDim,
                   E, ExcludeInnerCells, 0, GBarNearby);
     GetGBarDistant(R, Rho, k, kBloch, Gamma, LDim, E, 0, GBarDistant);
   };
   
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

} 

} // namespace scuff
